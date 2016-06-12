#include "feature.h"
#include "../gauss_iir/gauss_iir.h"

#include "sift_base.h"
#include "sift_pca_base.h"
#include "esift_base.h"
#include "esift_pca_base.h"
#include "sc_base.h"
#include "sc_pca_base.h"
#include "pca_base.h"
#include "pca_eig.h"

const int SegLocSize=4;
const int SegOriSize=6;
const int SegSiftSize=SegLocSize*SegLocSize*SegOriSize;
const int LocSize=4;
const int OriSize=8;
const int SiftSize=LocSize*LocSize*OriSize;
const int MagFactor=2;
const float MaxIndexVal = 0.2;  /* Good value is 0.2 */

const int ELocSize=4; 
const int EOriSize=6;
const int ESiftSize=ELocSize*ELocSize*EOriSize;

const int SrSize=3;
const int ScSize=4;
const int SLocSize=4;
const int SOriSize=8;
//const int ShapeSize=SLocSize*SLocSize*SOriSize;
const int ShapeSize=SrSize*ScSize*SOriSize;

const int GPLEN=3042;
const int PCALEN=36;

void MakeKeypointSample(DARY *grad, DARY *ori, float *vec, int size, int oriSize, int locSize);
void AddSample(float *index, DARY *grad, DARY *orim, float angle, int r, int c, float rpos, 
	       float cpos,float rx, float cx, int oriSize, int locSize, float max_angle);
void PlaceInIndex(float *index,
		  float mag, float ori, float rx, float cx, int oriSize, int locSize, float max_angle);
void KeySample(float *index, DARY *grad, DARY *ori, float angle, int oriSize, int locSize, float max_angle);
void NormalizeVect(float *vec, int len);
void NormalizeEVect(float *vec, int len, float norm);

void normalizeFeature(float *vec, int size){

  NormalizeVect(vec, size); 
  int intval, changed = FALSE;
   for (int i = 0; i < size; i++)
    if (vec[i] > MaxIndexVal) { 
    vec[i] = MaxIndexVal;
    changed = TRUE;
    }
    if (changed)
      NormalizeVect(vec, size); 
  
    
  /* Convert float vector to integer. Assume largest value in normalized
    vector is likely to be less than 0.5. */
    for (int i = 0; i < size; i++) {
      intval = (int) (512.0 * vec[i]);
      vec[i] = (255 < intval) ? 255 : intval;
    }
}


/* David Lowe's code*/

/* Increment the appropriate locations in the index to incorporate
   this image sample.  The location of the sample in the index is (rx,cx). 
*/

void PlaceInLogPolIndex(float *index,
			float mag, float ori, float rx, float cx, int oriSize, int rSize, int cSize)
{
   int r, c, ort, ri, ci, oi, rindex, cindex, oindex, rcindex;
   float oval, rfrac, cfrac, ofrac, rweight, cweight, oweight;
   
   oval = oriSize * ori / (M_2PI);
   
   ri = (int)((rx >= 0.0) ? rx : rx - 1.0);  /* Round down to next integer. */
   ci = (int)((cx >= 0.0) ? cx : cx - 1.0);
   oi = (int)((oval >= 0.0) ? oval : oval - 1.0);
   rfrac = rx - ri;         /* Fractional part of location. */
   cfrac = cx - ci;
   ofrac = oval - oi; 
   /*   assert(ri >= -1  &&  ri < IndexSize  &&  oi >= 0  &&  oi <= OriSize  && rfrac >= 0.0  &&  rfrac <= 1.0);*/
   //cout << ri << " " << ci << " " << oi << endl;
   /* Put appropriate fraction in each of 8 buckets around this point
      in the (row,col,ori) dimensions.  This loop is written for
      efficiency, as it is the inner loop of key sampling. */
   for (r = 0; r < 2; r++) {
      rindex = ri + r;
      if (rindex >=0 && rindex < rSize) {
         rweight = mag * ((r == 0) ? 1.0 - rfrac : rfrac);
         
         for (c = 0; c < 2; c++) {
            cindex = ci + c;
	    if(cindex >= cSize)
	      cindex=0;	    
            if (cindex >=0 && cindex < cSize) {
               cweight = rweight * ((c == 0) ? 1.0 - cfrac : cfrac);
               rcindex=(rindex*cSize+cindex)<<3;//remember when you change the orientation number
               for (ort = 0; ort < 2; ort++) {
                  oindex = oi + ort;
                  if (oindex >= oriSize)  /* Orientation wraps around at PI. */
                     oindex = 0;
                  oweight = cweight * ((ort == 0) ? 1.0 - ofrac : ofrac);
                  //cout << rcindex+oindex<< endl;
		  index[rcindex+oindex]+=oweight;
               }
            }  
         }
      }
   } 
}

/* Given a sample from the image gradient, place it in the index array.
*/
void AddLogPolSample(float *index,
		     DARY *grad, DARY *orim, float angle, int r, int c, float rpos, float cpos,
		     float rx, float cx, int oriSize, int rSize, int cSize)
{
    float mag, ori;
    
    /* Clip at image boundaries. */
    if (r < 0  ||  r >= (int)grad->y()  ||  c < 0  ||  c >= (int)grad->x())
       return;
    
    mag = patch_mask->fel[r][c] * grad->fel[r][c];
    //mag = grad->fel[r][c];
    /* Subtract keypoint orientation to give ori relative to keypoint. */
    ori = orim->fel[r][c]-angle;
    
    /* Put orientation in range [0, 2*PI].  If sign of gradient is to
       be ignored, then put in range [0, PI]. */
    
    while (ori > M_2PI)
      ori -= M_2PI;
    while (ori < 0.0)
      ori += M_2PI;     
    PlaceInLogPolIndex(index, mag, ori, rx, cx, oriSize, rSize, cSize);
} 
void KeyLogPolSample(float *index, DARY *grad, DARY *ori, float angle, int oriSize, int rSize, int cSize)
{
   int i, j, iradius;
   float rspacing,cspacing, rpos, cpos, lrpos, lcpos, rx, cx;
   
   /* Radius of index sample region must extend to diagonal feature of
      index patch plus half sample for interpolation. */
   
   iradius = grad->x()>>1;
   rspacing = (rSize) / (iradius);
   cspacing = (cSize) / (M_2PI);
   float sine = sin(-angle);
   float cosine = cos(-angle);

   //   printf("spacing %f, scale %f, radius %d\n",spacing,scale, iradius);
   // Examine all points from the gradient image that could lie within the index square. 
   for (i = -iradius; i <= iradius; i++)
     for (j = -iradius; j <= iradius; j++) {
       

       rpos = (cosine*i + sine*j);// coreect with dominant angle
       cpos = (-sine*i + cosine*j);
       lcpos=(M_PI+atan2(rpos,cpos))*cspacing;
       lrpos=log(1+sqrt((float)i*i+j*j)/iradius)*rSize;
       
       // Compute location of sample in terms of real-valued index array
       //  coordinates.  Subtract 0.5 so that rx of 1.0 means to put full
       // weight on index[1] (e.g., when rpos is 0 and IndexSize is 3. 
       rx = lrpos;// + (rSize - 1) / 2.0;
       cx = lcpos;// + (cSize - 1) / 2.0;
       //cout << rx << " "<< cx << endl;
       
       //cout <<"rx " << rx << " " << cx << endl;
       // Test whether this sample falls within boundary of index patch. 
       if (rx > -1.0 && rx < (float) rSize  &&
	   cx > -1.0 && cx < (float) cSize)
         //cout << "in" << cpos << " " << rpos << endl;
	 AddLogPolSample(index, grad, ori, angle, iradius + i, iradius + j, lrpos, lcpos,
		   rx, cx, oriSize, rSize, cSize);
   } 
}


void PlaceInIndex(float *index,
		  float mag, float ori, float rx, float cx, int oriSize, int locSize,float max_angle)
{
   int r, c, ort, ri, ci, oi, rindex, cindex, oindex, rcindex;
   float oval, rfrac, cfrac, ofrac, rweight, cweight, oweight;
   
   oval = oriSize * ori / (max_angle);
   
   ri = (int)((rx >= 0.0) ? rx : rx - 1.0);  /* Round down to next integer. */
   ci = (int)((cx >= 0.0) ? cx : cx - 1.0);
   oi = (int)((oval >= 0.0) ? oval : oval - 1.0);
   rfrac = rx - ri;         /* Fractional part of location. */
   cfrac = cx - ci;
   ofrac = oval - oi; 

   //cout << ri<< " "<< ci << " " << oi << endl;
   /*   assert(ri >= -1  &&  ri < IndexSize  &&  oi >= 0  &&  oi <= OriSize  && rfrac >= 0.0  &&  rfrac <= 1.0);*/
   
   /* Put appropriate fraction in each of 8 buckets around this point
      in the (row,col,ori) dimensions.  This loop is written for
      efficiency, as it is the inner loop of key sampling. */
   for (r = 0; r < 2; r++) {
      rindex = ri + r;
      if (rindex >=0 && rindex < locSize) {
         rweight = mag * ((r == 0) ? 1.0 - rfrac : rfrac);         
         for (c = 0; c < 2; c++) {
            cindex = ci + c;
            if (cindex >=0 && cindex < locSize) {
               cweight = rweight * ((c == 0) ? 1.0 - cfrac : cfrac);
               rcindex=(rindex*locSize+cindex)*oriSize;//remember when you change the orientation number
               for (ort = 0; ort < 2; ort++) {
		 oindex = oi + ort;
                  if (oindex >= oriSize)  /* Orientation wraps around at PI. */
		    oindex = 0;
                  oweight = cweight * ((ort == 0) ? 1.0 - ofrac : ofrac);
                  
		  // if((rcindex+oindex)>=SegSiftSize)cout << "indexes " << oindex << " "<< rindex <<"  " << cindex<< " " << oweight<< " " << rcindex+oindex<< endl;
		  index[rcindex+oindex]+=oweight;
               }
            }  
         }
      }
   } 
}
 
  
/* Given a sample from the image gradient, place it in the index array.
*/
void AddSample(float *index,
	       DARY *grad, DARY *orim, float angle, int r, int c, float rpos, float cpos,
	       float rx, float cx, int oriSize, int locSize, float max_angle)
{
    float mag, ori;
    
    /* Clip at image boundaries. */
    if (r < 0  ||  r >= (int)grad->y()  ||  c < 0  ||  c >= (int)grad->x())
       return;
    
    /* Compute Gaussian weight for sample, as function of radial distance
       from center.  Sigma is relative to half-width of index. */
    //sigma =  0.5*(IndexSize+1)*(IndexSize+1);
    //sigma = (IndexSize+1)*(IndexSize+1);
    //weight = 0.6*exp(- (rpos * rpos + cpos * cpos) / (sigma) );
    //cout << "rpos "<< rpos << " cpos "<< cpos << " weight " << weight<< endl;
    //mag = weight * grad->fel[r][c];
    //mag = grad->fel[r][c];
    mag = patch_mask->fel[r][c] * grad->fel[r][c];
    /* Subtract keypoint orientation to give ori relative to keypoint. */
    ori = orim->fel[r][c]-angle;
    
    /* Put orientation in range [0, 2*PI].  If sign of gradient is to
       be ignored, then put in range [0, PI]. */
    
    while (ori > max_angle)
      ori -= max_angle;
    while (ori < 0.0)
      ori += max_angle;     
    PlaceInIndex(index, mag, ori, rx, cx, oriSize, locSize, max_angle);
} 

/* Add features to vec obtained from sampling the grad and ori images
   for a particular scale.  Location of key is (scale,row,col) with respect
   to images at this scale.  We examine each pixel within a circular
   region containing the keypoint, and distribute the gradient for that
   pixel into the appropriate bins of the index array.
*/
void KeySampleOrig(float *index, DARY *grad, DARY *ori, int x, int y, int iradius, float angle, int oriSize, int locSize, float max_angle)
{
  int i, j, r,c;
  float spacing, rpos, cpos, rx, cx,mag,orim;
   
   /* Radius of index sample region must extend to diagonal feature of
      index patch plus half sample for interpolation. */
   
   spacing = (locSize + 1) / (2.0*iradius);
   float sine = sin(-angle);
   float cosine = cos(-angle);
   //   printf("spacing %f, scale %f, radius %d\n",spacing,scale, iradius);
   // Examine all points from the gradient image that could lie within the index square. 
   int patch_rad=(patch_mask->x()>>1) ;
   //cout << iradius << (patch_mask->x()>>1) << endl;
   for (i = -iradius; i <= iradius; i++)
     for (j = -iradius; j <= iradius; j++) {
       
       rpos = (cosine*i + sine*j) * spacing;
       cpos = (-sine*i + cosine*j) * spacing;
        
       // Compute location of sample in terms of real-valued index array
       //  coordinates.  Subtract 0.5 so that rx of 1.0 means to put full
       // weight on index[1] (e.g., when rpos is 0 and IndexSize is 3. 
       rx = rpos + (locSize - 1) / 2.0;
       cx = cpos + (locSize - 1) / 2.0;
       
       r=y + i;
       c=x + j;
       //cout <<"rx " << x<< " " << y << " " << r << " " << c << " " << iradius<< endl;
       // Test whether this sample falls within boundary of index patch. 
       if (rx > -1.0 && rx < (float) locSize  &&
	   cx > -1.0 && cx < (float) locSize &&
	   r >= 0  &&  r < (int)grad->y()  &&   c >= 0  &&  c < (int)grad->x()){
         //cout << "in" << cpos << " " << rpos << endl;
         if(0 && iradius==patch_rad)
	   mag = patch_mask->fel[iradius+i][iradius+j] * grad->fel[r][c];
	 else mag =  grad->fel[r][c];
         if(max_angle>100)
           orim=ori->fel[r][c];
         else {
           orim=ori->fel[r][c]-angle;
           while (orim > max_angle)orim -= max_angle;
           while (orim < 0.0)orim += max_angle;
         }
	 PlaceInIndex(index, mag, orim, rx, cx, oriSize, locSize,max_angle);
       }
   } 
}

void KeySample(float *index, DARY *grad, DARY *ori, float angle, int oriSize, int locSize, float max_angle)
{
  int i, j, iradius,r,c;
  float spacing, rpos, cpos, rx, cx,mag,orim;
   
   /* Radius of index sample region must extend to diagonal feature of
      index patch plus half sample for interpolation. */
 //   grad->writePNG("grad.png");
 //    ori->writePNG("gori.png");
   
   iradius = grad->x()>>1;
   spacing = (locSize + 1) / (2.0*iradius);
   float sine = sin(-angle);
   float cosine = cos(-angle);
   //   printf("spacing %f, scale %f, radius %d\n",spacing,scale, iradius);
   // Examine all points from the gradient image that could lie within the index square. 
   for (i = -iradius; i <= iradius; i++)
     for (j = -iradius; j <= iradius; j++) {
       
       rpos = (cosine*i + sine*j) * spacing;
       cpos = (-sine*i + cosine*j) * spacing;
        
       // Compute location of sample in terms of real-valued index array
       //  coordinates.  Subtract 0.5 so that rx of 1.0 means to put full
       // weight on index[1] (e.g., when rpos is 0 and IndexSize is 3. 
       rx = rpos + (locSize - 1) / 2.0;
       cx = cpos + (locSize - 1) / 2.0;
       
       r=iradius + i;
       c=iradius + j;
       //cout <<"rx " << rx << " " << cx << endl;
       // Test whether this sample falls within boundary of index patch. 
       if (rx > -1.0 && rx < (float) locSize  &&
           cx > -1.0 && cx < (float) locSize &&
           r >= 0  &&  r < (int)grad->y()  &&   c >= 0  &&  c < (int)grad->x()){
         //cout << "in" << cpos << " " << rpos << endl;
         if(grad->size())//==patch_mask->size())
           mag =  grad->fel[r][c];// * patch_mask->fel[r][c];
         else mag =  grad->fel[r][c];
         if(mag>0.0){
           orim=ori->fel[r][c]-angle;
           while (orim > max_angle)
             orim -= max_angle;
           while (orim < 0.0)
             orim += max_angle;     
           PlaceInIndex(index, mag, orim, rx, cx, oriSize, locSize,max_angle);
         }
      }
   } 
}
void KeySample(FD *fd, vector<edge> points, float **grad, float **ori, int oriSize, int locSize, float max_angle){
  
  
  float *index = fd->getVec();
  float scale=fd->getScale();
  float angle=fd->getAngle();
  int mx=(int)fd->getX();
  int my=(int)fd->getY();
  
  uint i;
  float spacing, rpos, cpos, rx, cx, orim;
   
   /* Radius of index sample region must extend to diagonal feature of
  index patch plus half sample for interpolation. */
   
  spacing = (locSize + 1) / (2.0*scale);
  float sine = sin(-angle);
  float cosine = cos(-angle);

   //   printf("spacing %f, scale %f, radius %d\n",spacing,scale, iradius);
   // Examine all points from the gradient image that could lie within the index square. 
  for (i = 0; i <points.size(); i++){
     
    rpos = (cosine*(points[i].b-my) + sine*(points[i].a-mx)) * spacing;
    cpos = (-sine*(points[i].b-my) + cosine*(points[i].a-mx)) * spacing;
        
       // Compute location of sample in terms of real-valued index array
       //  coordinates.  Subtract 0.5 so that rx of 1.0 means to put full
       // weight on index[1] (e.g., when rpos is 0 and IndexSize is 3. 
    rx = rpos + (locSize - 1) / 2.0;
    cx = cpos + (locSize - 1) / 2.0;
       
       // Test whether this sample falls within boundary of index patch. 
    if (rx > -1.0 && rx < (float) locSize  &&
        cx > -1.0 && cx < (float) locSize){
         //cout << "in" << cpos << " " << rpos << endl;
      orim=ori[points[i].b][points[i].a]-angle;
      while (orim > max_angle)
        orim -= max_angle;
      while (orim < 0.0)
        orim += max_angle;     
	 //cout <<sy[i] << " " <<sx[i] << " rx " << rx << " " << cx << " rpos "<<rpos << " cpos " << cpos << " orim "<< orim <<  endl;
      PlaceInIndex(index, grad[points[i].b][points[i].a], orim, rx, cx, oriSize, locSize,max_angle);
        }else {
	 //cout << i << " of "<< sx.size() << " " << rx << " "<< cx << " " << locSize << endl; 
        }
  } 
  
  
  normalizeFeature(index, fd->getSize());
  
}

void KeySample(float *index, vector<int> sx, vector<int> sy, vector<float> grad, vector<float> ori, float angle, float scale, int oriSize, int locSize, float max_angle)
{
  uint i;
  float spacing, rpos, cpos, rx, cx, orim;
   
   /* Radius of index sample region must extend to diagonal feature of
      index patch plus half sample for interpolation. */
   
    spacing = (locSize + 1) / (2.0*scale);
   float sine = sin(-angle);
   float cosine = cos(-angle);

   //   printf("spacing %f, scale %f, radius %d\n",spacing,scale, iradius);
   // Examine all points from the gradient image that could lie within the index square. 
   for (i = 0; i <sx.size(); i++){
     
       rpos = (cosine*sy[i] + sine*sx[i]) * spacing;
       cpos = (-sine*sy[i] + cosine*sx[i]) * spacing;
        
       // Compute location of sample in terms of real-valued index array
       //  coordinates.  Subtract 0.5 so that rx of 1.0 means to put full
       // weight on index[1] (e.g., when rpos is 0 and IndexSize is 3. 
       rx = rpos + (locSize - 1) / 2.0;
       cx = cpos + (locSize - 1) / 2.0;
       
       // Test whether this sample falls within boundary of index patch. 
       if (rx > -1.0 && rx < (float) locSize  &&
	   cx > -1.0 && cx < (float) locSize){
         //cout << "in" << cpos << " " << rpos << endl;
	 orim=ori[i]-angle;
	 while (orim > max_angle)
	   orim -= max_angle;
	 while (orim < 0.0)
	   orim += max_angle;     
	 //cout <<sy[i] << " " <<sx[i] << " rx " << rx << " " << cx << " rpos "<<rpos << " cpos " << cpos << " orim "<< orim <<  endl;
	 PlaceInIndex(index, grad[i], orim, rx, cx, oriSize, locSize,max_angle);
       }else {
	 //cout << i << " of "<< sx.size() << " " << rx << " "<< cx << " " << locSize << endl; 
       }
   } 
}

 
/* Normalize length of vec to 1.0.
*/
void NormalizeVect(float *vec, int len)
{
   int i;
   float val, fac, sqlen = 0.0;

   for (i = 0; i < len; i++) {
     val = vec[i];
     sqlen += val * val;
   }
   if(sqlen<0.0001){
     sqlen+=1;
   }
   fac = 1.0 / sqrt(sqlen);
  for (i = 0; i < len; i++)
     vec[i] *= fac;
}

void NormalizeEVect(float *vec, int len, float norm)
{
   float fac, val, sqlen = 0.0;

   for (int i = 0; i < len; i++) {
     val = vec[i];
      sqlen += val * val;
      //sqlen += fabs(vec[i]);
   }
   if(sqlen<0.0001){
     sqlen+=1;
      
   }
      fac = norm / sqrt(sqlen);
   //fac=sqlen/len;
   for (int i = 0; i < len; i++){
     //vec[i] *= fac;
     vec[i]=(int)vec[i];
   }
}
   


void symmetricSift(FD *fd, float x){

  float *symVec=new float[fd->getSize()];
 
  float *vec=fd->getVec();
  int locSize,oriSize;
  if((fd->getType() & DSEGM)!=0){
    locSize=SegLocSize;
    oriSize=SegOriSize; 
  }else{
    locSize=LocSize;
    oriSize=OriSize;   
  }
  for (int r=0; r<locSize; r++)
    for (int c=0; c<locSize; c++)
      for (int o=0; o<oriSize; o++){
	symVec[(r*locSize+c)*oriSize+o] =
	  vec[((locSize-1-r)*locSize+c)*oriSize+((o==0)?0:oriSize-o)];
      }        
  
  float angle=M_PI-fd->getAngle();//M_PI is correct <= atan2(dy=1,dx=0)=M_PI/2;
  float mangle=M_PI-fd->getMotion();//M_PI is correct <= atan2(dy=1,dx=0)=M_PI/2;
  while(angle<-M_PI)angle+=M_2PI;
  while(angle>M_PI)angle-=M_2PI;
  while(mangle<-M_PI)mangle+=M_2PI;
  while(mangle>M_PI)mangle-=M_2PI;
  fd->setMotion(mangle);
  fd->setAngle(angle);
  fd->setSymm(1);
  
  memcpy(vec, symVec, sizeof(float)*fd->getSize());
  fd->setMi(fd->getMi11(),-fd->getMi12(),-fd->getMi21(),fd->getMi22());
  fd->setX_Y(x-fd->getX(),fd->getY());

  delete []symVec;
} 

void computeSift(DARY *img, FeatureDescriptor *ds){
  
  DARY * grad = new DARY(patch_mask->x(),patch_mask->x(),FLOAT1);
  DARY * gori = new DARY(patch_mask->x(),patch_mask->x(),FLOAT1);    
  
  gradAngle(img,grad,gori);
  ds->allocVec(SiftSize);
  float *vec = ds->getVec();
  KeySample(vec, grad, gori, 0.0, OriSize, LocSize, M_2PI);
  normalizeFeature(vec, SiftSize); 
  
  ds->setType((ADETECTOR & ds->getType())|DSIFT);
  delete grad;delete gori;
  //int sift_pca_size=128;
   //ds->pca(sift_pca_size, sift_pca_avg, sift_pca_base);	
}


void computeSift(Segment *seg, FeatureDescriptor *ds){
  
  ds->allocVec(SegSiftSize);
  float *vec = ds->getVec();
  KeySample(vec, seg->x, seg->y, seg->grad, seg->angle, seg->ang_mean, seg->scale, SegOriSize, SegLocSize, M_PI);
  
  NormalizeVect(vec, SegSiftSize); 
  int intval, changed = FALSE;
  for (int i = 0; i < SegSiftSize; i++)
    if (vec[i] > MaxIndexVal) { 
       vec[i] = MaxIndexVal;
       changed = TRUE;
     }
  if (changed) 
    NormalizeVect(vec, SegSiftSize); 
  /* Convert float vector to integer. Assume largest value in normalized
     vector is likely to be less than 0.5. */
  for (int i = 0; i < SegSiftSize; i++) {
    intval = (int) (512.0 * vec[i]);
    vec[i] = (255 < intval) ? 255 : intval;
  }
  ds->setX_Y(seg->x_mean,seg->y_mean);
  ds->setAngle(seg->ang_mean);
  ds->setScale(seg->scale);

  ds->setType(DSEGM | DSIFT);
  //int sift_pca_size=128;
   //ds->pca(sift_pca_size, sift_pca_avg, sift_pca_base);	
}



void computeSift(DARY *dx, DARY *dy ,DARY *grad, DARY *gori , FeatureDescriptor *ds){
  ds->allocVec(SiftSize);
  float *vec = ds->getVec();
  float angle=ds->getAngle();
  //if(angle<-M_PI || angle>M_PI)angle=0;
  //float t=vec[0];
 int x=grad->x()>>1;
 int y=grad->y()>>1;
 int radius =x;
  
 //KeySampleOrig(vec, grad, gori, x, y, radius, angle, OriSize, LocSize, M_2PI);
 KeySample(vec, grad, gori, angle, OriSize, LocSize, M_2PI);
  //NormalizeVect(vec, SiftSize);
  normalizeFeature(vec, SiftSize);
  ds->setType((ADETECTOR & ds->getType())|DSIFT);
 // ds->Cout(50);
  //grad->writePNG("grad1.png");
  //if(vec[0]<0){cout << t << "  ";ds->Cout(5);grad->writePNG("grad.png");getchar();}
 //int sift_pca_size=128;
   //ds->pca(sift_pca_size, sift_pca_avg, sift_pca_base);	
}

void computeSiftOrig(DARY *dx, DARY *dy ,DARY *grad, DARY *gori , FeatureDescriptor *ds){
  ds->allocVec(SiftSize);
  float *vec = ds->getVec();
  float angle=ds->getAngle();
  int x=(int)(0.5+ds->getX());
  int y=(int)(0.5+ds->getY());
  int radius=(int)ds->getScale();
  //if(angle<-M_PI || angle>M_PI)angle=0;
  //float t=vec[0];
  //grad->writePNG("grad.png");getchar();
  KeySampleOrig(vec, grad, gori, x, y, radius, angle, OriSize, LocSize, M_2PI);
  normalizeFeature(vec, SiftSize);
  //NormalizeVect(vec, SiftSize);
  //NormalizeEVect(vec, SiftSize, 1.0);

  ds->setType((ADETECTOR & ds->getType())|DSIFT);
  //if(vec[0]<0){cout << t << "  ";ds->Cout(5);grad->writePNG("grad.png");getchar();}
 //int sift_pca_size=128;
  //ds->pca(sift_pca_size, sift_pca_avg, sift_pca_base);	
}
void computeSiftOrig(DARY *dx, DARY *dy, DARY *grad, DARY *gori , FeatureDescriptor *ds, int oriSize, int locSize, float max_gori){
  const int siftSize=locSize*locSize*oriSize;
  //cout << ds->getSize() << " " << locSize << " " << oriSize << endl;

  ds->allocVec(siftSize);
  //cout << siftSize << " " << locSize << " " << oriSize << endl;
  float *vec = ds->getVec();
  float angle=ds->getAngle();
  int x=(int)(0.5+ds->getX());
  int y=(int)(0.5+ds->getY());
  int radius=(int)ds->getScale();
 

  //if(angle<-M_PI || angle>M_PI)angle=0;
  //float t=vec[0];
  //ds->Cout();
  //grad->writePNG("grad.png");getchar();
  KeySampleOrig(vec, grad, gori, x, y, radius, angle, oriSize, locSize, M_2PI);
  //normalizeFeature(vec, siftSize);
 // ds->Cout(40);
  //NormalizeVect(vec, SiftSize);
  //NormalizeEVect(vec, SiftSize, 512.0);

  //cout << ds->getType() << " ";
  ds->setType((ADETECTOR & ds->getType())|DCOLOR);
  // grad->writePNG("grad2.png");
  //cout << ds->getType() <<endl;getchar();
  //if(vec[0]<0){cout << t << "  ";ds->Cout(5);grad->writePNG("grad.png");getchar();}
 //int sift_pca_size=128;
  //ds->pca(sift_pca_size, sift_pca_avg, sift_pca_base);	
}
 
 
 
void computeESift(DARY *img, FeatureDescriptor *ds){
  
  DARY * grad = new DARY(patch_mask->x(),patch_mask->x(),FLOAT1);
  DARY * gori = new DARY(patch_mask->x(),patch_mask->x(),FLOAT1);    
  
  gradAngle(img,grad,gori);
  ds->allocVec(ESiftSize);//EIndexSize*EIndexSize*OriSize;
  float *vec = ds->getVec();

  KeySample(vec, grad, gori, 0.0, EOriSize, ELocSize, M_2PI);
  //grad->writePNG("grad.png");cout << "patch 1" << endl;getchar();
  delete grad;delete gori;
  NormalizeEVect(vec, ESiftSize, 1.0); 
  ds->setType((ADETECTOR & ds->getType())|DGLOH);
  int esift_pca_size=128; 
  ds->pca(esift_pca_size, esift_pca_avg, esift_pca_base);

}


void computeESift(DARY *dx, DARY *dy ,DARY *grad, DARY *gori , FeatureDescriptor *ds){
  ds->allocVec(ESiftSize);//EIndexSize*EIndexSize*OriSize;
  float *vec = ds->getVec();
  float angle=ds->getAngle();
  if(angle<-M_PI || angle>M_PI)angle=0;
  KeySample(vec, grad, gori, angle, EOriSize, ELocSize, M_2PI);
  //grad->writePNG("grad.png");cout << "patch 2" << endl;getchar();
  NormalizeEVect(vec, ESiftSize, 1.0); 
  ds->setType((ADETECTOR & ds->getType())|DGLOH);
  int esift_pca_size=128; 
  ds->pca(esift_pca_size, esift_pca_avg, esift_pca_base);	
}


/******************SHAPE CONTEXT******************************************/


void computeShape(DARY *img, FeatureDescriptor *ds){

  DARY * grad = new DARY(patch_mask->x(),patch_mask->x(),FLOAT1);
  DARY * ori = new DARY(patch_mask->x(),patch_mask->x(),FLOAT1);    
  DARY *dx = new DARY(img->y(),img->x(),FLOAT1);
  DARY *dy = new DARY(img->y(),img->x(),FLOAT1);
  DARY *edge = new DARY(img->y(),img->x(),FLOAT1);
  //float angle=ds->getAngle();
  //canny(img,grad,ori,5,15);
  dX2(img,dx);
  dY2(img,dy);
  for(uint j=0;j<grad->y();j++){
    for(uint i=0;i<grad->x();i++){
      grad->fel[j][i]=sqrt(dx->fel[j][i]*dx->fel[j][i]+dy->fel[j][i]*dy->fel[j][i]); 
      ori->fel[j][i]=atan2(dy->fel[j][i],dx->fel[j][i]);
    }
  } 
  cannyEdges(dx, dy, grad, edge, 5, 15);
  //grad->write("edge.pgm");cout << "OK"<< endl;getchar();
  delete dx; delete dy;delete grad;


  ds->allocVec(ShapeSize);
  float *vec = ds->getVec();
  //KeyLogPolSample(vec, edge, ori,  angle, SOriSize, SrSize, ScSize);
  //getchar();
   KeySample(vec, edge, ori,0.0, SOriSize, SLocSize, M_2PI);
  delete edge;delete ori;
  NormalizeVect(vec, ShapeSize); 
  int intval, changed = FALSE;
  for (int i = 0; i < ShapeSize; i++)
    if (vec[i] > MaxIndexVal) { 
      vec[i] = MaxIndexVal;
       changed = TRUE;
     }
  if (changed)
    NormalizeVect(vec, ShapeSize); 
  /* Convert float vector to integer. 
     Assume largest value in normalized
     vector is likely to be less than 0.5. */
  for (int i = 0; i < ShapeSize; i++) {
    intval = (int) (512.0 * vec[i]);
    vec[i] = (255 < intval) ? 255 : intval;
  }
   
  ds->setType((ADETECTOR & ds->getType())|DSHAPE);
  int sc_pca_size=36;
  ds->pca(sc_pca_size, sc_pca_avg, sc_pca_base);	

} 
  
void computeShape(DARY *dx, DARY *dy ,DARY *grad, DARY *ori, FeatureDescriptor *ds){

  DARY *edge = new DARY(dx->y(),dx->x(),FLOAT1);

  cannyEdges(dx, dy, grad, edge,   5, 15);
  //  edge->write("edge.pgm");cout << "OK"<< endl;getchar();

  float angle=ds->getAngle();
  if(angle<-M_PI || angle>M_PI)angle=0;

  ds->allocVec(ShapeSize);
  float *vec = ds->getVec();
  KeyLogPolSample(vec, edge, ori,  angle, SOriSize, SrSize, ScSize);
 //getchar();

  //KeySample(vec, edge, ori, angle, SOriSize, SLocSize);
  delete edge;
  NormalizeVect(vec, ShapeSize); 
  int intval, changed = FALSE;
  for (int i = 0; i < ShapeSize; i++)
    if (vec[i] > MaxIndexVal) { 
      vec[i] = MaxIndexVal;
       changed = TRUE;
     }
  if (changed)
    NormalizeVect(vec, ShapeSize); 
  /* Convert float vector to integer. 
     Assume largest value in normalized
     vector is likely to be less than 0.5. */
  for (int i = 0; i < ShapeSize; i++) {
    intval = (int) (512.0 * vec[i]);
    vec[i] = (255 < intval) ? 255 : intval;
  }
 
  ds->setType((ADETECTOR & ds->getType())| DSHAPE);
  //int sc_pca_size=36;
  //ds->pca(sc_pca_size, sc_pca_avg, sc_pca_base);	

} 


void computePcaSift(DARY *img, FeatureDescriptor *ds){
  

  float *tvec = new float[GPLEN];
  uint count=0;
  for(int j=1;j<patch_mask->x()-1;j++){
    for(int i=1;i<patch_mask->x()-1;i++){      
      tvec[count++]=img->fel[j][i+1]-img->fel[j][i-1];
      tvec[count++]=img->fel[j+1][i]-img->fel[j-1][i];
    }
  }

  NormalizeEVect(tvec, GPLEN, 1.0); 


  for (int i = 0; i < GPLEN; i++) {
    tvec[i] -= avgs[i];
  }

  ds->allocVec(PCALEN);//EIndexSize*EIndexSize*OriSize;
  float *vec = ds->getVec();

  for (int ldi = 0; ldi < PCALEN; ldi++) {
    for (int x = 0; x < GPLEN; x++)
      vec[ldi] += eigs[ldi][x] * tvec[x];
  }
  ds->setType((ADETECTOR & ds->getType())| DPCA);
  delete tvec;
}

void computeNON(DARY *img, FeatureDescriptor *ds){}
