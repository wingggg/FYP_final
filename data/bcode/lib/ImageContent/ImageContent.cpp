// Implementation de la classe image
#include "imageread.h"
#include "imageContent.h"
using namespace std;
using namespace imgcon;

float ImageContent::getValue(float x, float y){

  uint ix=(int)floor(x);
  uint iy=(int)floor(y);
  float dx = x-ix;
  float dy = y-iy; 
  float val=0;
  if(ix>=0 && iy>=0 && ix<x_size-1 && iy<y_size-1){
    val =  (1.0 - dy)*((1.0 - dx)*fel[iy][ix] +
		       dx*fel[iy][ix+1]) +
      dy*((1.0 - dx)*fel[iy+1][ix] +
	  dx *fel[iy+1][ix+1]);    
  } else if(ix==(x_size-1) && iy<y_size-1){
    val= (1.0 - dy)*fel[iy][ix] +
      dy*(fel[iy+1][ix]);
  }else if(iy==(y_size-1) && ix<x_size-1){
    val= (1.0 - dx)*fel[iy][ix] +
      dx*(fel[iy][ix+1]);
  }else if(iy==(y_size-1) && ix==(x_size-1)){
    val = fel[iy][ix];
  }else val=0;
  return val;
}

  

/****************************************************************/
int verif_format(const char* name)
{


  FILE* ident;
  //char *buf1 = new char[2];
  unsigned char buf1[PNG_BYTES_TO_CHECK];

  if((  ident = fopen(name,"r") ) == NULL) {cout  << endl<<"\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\nUnable to open " << name << endl;   return(0); }
 
 fread(buf1,1,PNG_BYTES_TO_CHECK,ident);
 fclose(ident);
 
 if(!png_sig_cmp(&buf1[0], (png_size_t)0, PNG_BYTES_TO_CHECK))return(4);
 else if((buf1[0]) != 'P')  return(0);
 else if(buf1[1] == '5') return(1);
 else if( buf1[1] == 'g' ) return(2);
 else if(buf1[1] == '6' ) return(3);
 
  return 0; 
}

/****************************************************************/
ImageContent::ImageContent(uint y_size_in,uint x_size_in, int buftype_in)
{
  y_size=0;x_size=0;buftype=0;
  init(y_size_in,x_size_in, buftype_in);
  buftype=buftype_in;
}
/****************************************************************/
ImageContent::ImageContent(ImageContent *im){
  y_size=0;x_size=0;buftype=0;
  init(im->y(),im->x(),im->getType());    
  buftype=im->getType();
  set(im);
}	   


/****************************************************************/
void ImageContent::init(uint y_size_in,uint x_size_in, int type_in)
{
  //cout << buftype << " " << type_in << endl;
    
  if(buftype !=  FLOAT1 && buftype !=  UCHAR1 && buftype !=  UCHAR1FLOAT1 && buftype !=  FLOAT3 && 
     buftype !=  UCHAR3 && buftype !=  UCHAR3FLOAT3 && buftype !=  UCHAR3FLOAT1)filename = new char[100];
  
  if(type_in==FLOAT1 || type_in==UCHAR1FLOAT1 || type_in==UCHAR3FLOAT1){
    if((buftype==FLOAT1 || buftype==UCHAR1FLOAT1 || buftype==UCHAR3FLOAT1) && x_size==x_size_in && y_size==y_size_in);
    else{
      if(buftype==FLOAT1 || buftype==UCHAR1FLOAT1 || buftype==UCHAR3FLOAT1){
        delete [] fel[0];delete [] fel;fel=NULL; 
      }
      x_size=x_size_in;
      y_size=y_size_in;
      tsize=y_size*x_size;
      fel = new float*[y_size];
      fel[0]=new float[tsize];
      for(uint i=1;i<y_size;i++)fel[i]=fel[0]+i*x_size;  
    }
  }

  if(type_in==UCHAR1 || type_in==UCHAR1FLOAT1){
    if((buftype==UCHAR1 || buftype==UCHAR1FLOAT1) && x_size==x_size_in && y_size==y_size_in);
    else{
      if(buftype==UCHAR1 || buftype==UCHAR1FLOAT1){
        delete [] bel[0];delete [] bel;bel=NULL; 
      }
	x_size=x_size_in;
	y_size=y_size_in;
	tsize=y_size*x_size;
	bel = new unsigned char*[y_size];
	bel[0]=new unsigned char[tsize];
	for(uint i=1;i<y_size;i++)bel[i]=bel[0]+i*x_size;  
    }
  }

  if(type_in==UCHAR3 || type_in==UCHAR3FLOAT1){
    if((buftype==UCHAR3 || buftype==UCHAR3FLOAT1) && x_size==x_size_in && y_size==y_size_in);
    else{
      if(buftype==UCHAR3 || buftype==UCHAR3FLOAT1){
        delete [] belr[0];delete [] belr;belr=NULL; 
        delete [] belg[0];delete [] belg;belg=NULL; 
        delete [] belb[0];delete [] belb;belb=NULL; 
      }
        filename = new char[100];
	x_size=x_size_in;
	y_size=y_size_in;
	tsize=y_size*x_size;

	belr = new unsigned char*[y_size];
	belr[0]=new unsigned char[tsize];
	for(uint i=1;i<y_size;i++)belr[i]=belr[0]+i*x_size;  

	belg = new unsigned char*[y_size];
	belg[0]=new unsigned char[tsize];
	for(uint i=1;i<y_size;i++)belg[i]=belg[0]+i*x_size;  

	belb = new unsigned char*[y_size];
	belb[0]=new unsigned char[tsize];
	for(uint i=1;i<y_size;i++)belb[i]=belb[0]+i*x_size;  
    }
  }

  if(type_in==FLOAT3){
    if((buftype==FLOAT3) && x_size==x_size_in && y_size==y_size_in)return;
    else{
      if(buftype==FLOAT3){
        delete [] felr[0];delete [] felr;felr=NULL; 
        delete [] felg[0];delete [] felg;felg=NULL; 
        delete [] felb[0];delete [] felb;felb=NULL; 
      }
	x_size=x_size_in;
	y_size=y_size_in;
	tsize=y_size*x_size;

	felr = new float*[y_size];
	felr[0]=new float[tsize];
	for(uint i=1;i<y_size;i++)felr[i]=felr[0]+i*x_size;  

	felg = new float*[y_size];
	felg[0]=new float[tsize];
	for(uint i=1;i<y_size;i++)felg[i]=felg[0]+i*x_size;  

	felb = new float*[y_size];
	felb[0]=new float[tsize];
	for(uint i=1;i<y_size;i++)felb[i]=felb[0]+i*x_size;  
    }
  }
  //buftype=type_in;
}


/****************************************************************/
void ImageContent::set(float val)
{
  if(buftype==FLOAT1 || buftype==UCHAR1FLOAT1 || buftype==UCHAR3FLOAT1){
     if(val==0.0){
      bzero(fel[0],tsize*sizeof(float));
    }else
      for (uint col = 0; col < tsize; col++){
	fel[0][col] = val;
      }    
  }
  if(buftype==UCHAR1 || buftype==UCHAR1FLOAT1 ){
    unsigned char value = (unsigned char)val;
    memset(bel[0],(int)value,tsize*sizeof(unsigned char));
  }
  if(buftype==UCHAR3 || buftype==UCHAR3FLOAT1 ){
    unsigned char value = (unsigned char)val;
    memset(belr[0],(int)value,tsize*sizeof(unsigned char));
    memset(belg[0],(int)value,tsize*sizeof(unsigned char));
    memset(belb[0],(int)value,tsize*sizeof(unsigned char));
  }
}

/****************************************************************/
void ImageContent::flipH()
{
  if(buftype==FLOAT1){
    float tmp;
    uint middle=(x_size>>1);
   for (uint row = 0; row < y_size; row++){
      for (uint col = 0; col < middle; col++){
	tmp=fel[row][col];
	fel[row][col]=fel[row][x_size-1-col];
	fel[row][x_size-1-col]=tmp;
      }    
    }   
   
  }else if(buftype==UCHAR1){
    unsigned char btmp;
    uint middle=(x_size>>1);
    for (uint row = 0; row < y_size; row++){
      for (uint col = 0; col < middle; col++){
	btmp=bel[row][col];
	bel[row][col]=bel[row][x_size-1-col];
	bel[row][x_size-1-col]=btmp;
      }    
    }       
  }
}
/****************************************************************/
void ImageContent::flipV()
{
  if(buftype==FLOAT1){
    float tmp;
    uint middle=(y_size>>1);
   for (uint col = 0; col < x_size; col++){
     for (uint row = 0; row <= middle; row++){
       tmp=fel[row][col];
       fel[row][col]=fel[y_size-1-row][col];
       fel[y_size-1-row][col]=tmp;
     }    
   }   
   
  }else if(buftype==UCHAR1){
    unsigned char btmp;
    uint middle=(y_size>>1);
   for (uint col = 0; col < x_size; col++){
     for (uint row = 0; row <= middle; row++){
       btmp=bel[row][col];
       bel[row][col]=bel[y_size-1-row][col];
       bel[y_size-1-row][col]=btmp;
     }    
   }   
  }
}

/****************************************************************/
void ImageContent::set(ImageContent *im)
{
// cout << FLOAT1 << " " << im->getType() << " " << buftype << " " << UCHAR3 << endl;
  if(x_size!=im->x() || y_size!=im->y()){
    cout << "ERROR   ImageContent::set, ingae size problem " << endl;exit(1);
  }
  if((im->getType()==UCHAR1 || im->getType()==UCHAR1FLOAT1) && (buftype==UCHAR1 || buftype==UCHAR1FLOAT1)){
    memcpy(bel[0],im->bel[0],tsize*sizeof(unsigned char)); 
  } 
  if((im->getType()==FLOAT1 || im->getType()==UCHAR3FLOAT1)&& buftype==FLOAT1){
     memcpy(fel[0],im->fel[0],tsize*sizeof(float));
  }
  if(im->getType()==UCHAR1 && (buftype==FLOAT1 || buftype==UCHAR1FLOAT1)){
    for (uint col = 0; col < tsize; col++){
      fel[0][col] = (float)im->bel[0][col];
    }  
  }
  if(im->getType()==FLOAT1 && buftype==UCHAR1){
      for (uint col = 0; col < tsize; col++){
      bel[0][col] = (unsigned char)im->fel[0][col];
      }    
  }else if(im->getType()==UCHAR3 && buftype==FLOAT1){
        for (uint col = 0; col < tsize; col++){
         // cout << "copy " << endl;
          fel[0][col] = (im->belr[0][col]+im->belg[0][col]+im->belb[0][col])/3.0;
        }  
  }else if(im->getType()==UCHAR3 && buftype==UCHAR3){
     memcpy(belr[0],im->belr[0],tsize*sizeof(unsigned char)); 
     memcpy(belg[0],im->belg[0],tsize*sizeof(unsigned char)); 
     memcpy(belb[0],im->belb[0],tsize*sizeof(unsigned char)); 
  }
  
} 

void ImageContent::set(unsigned char*bel_in, uint y, uint x){
  if((buftype==UCHAR1 || buftype==UCHAR1FLOAT1) && x==x_size && y==y_size){
    memcpy(bel[0],bel_in,tsize*sizeof(unsigned char)); 
    if(buftype==UCHAR1FLOAT1){
      for (uint col = 0; col < tsize; col++){
        fel[0][col] = (float)bel_in[col];
      }  
    }
  }
  
}
void ImageContent::set(float*fel_in, uint y, uint x){
  if((buftype==FLOAT1 || buftype==UCHAR1FLOAT1  || buftype==UCHAR3FLOAT1)  && x==x_size && y==y_size){
    memcpy(fel[0],fel_in,tsize*sizeof(float)); 
        
  }      
}
void ImageContent::median2d(ImageContent *output, int median_s)
{
    if(getType()!=UCHAR1FLOAT1 && getType()!=UCHAR3FLOAT1){
      cout<< "ImageContent::median2d: wrong input type "<< endl; 
      return;
    }
  //  cout << "ImageContent::median2d " << median_s << " " << x_size<<" " << y_size<<  endl;
    int med_size=median_s*median_s;
    int median_s2=median_s>>1;
    float *med_vec=new float[med_size];
    bzero(med_vec,med_size*sizeof(float));
    float *med_vec_tmp=new float[med_size];
    for(int i = 0 ;i<x_size-median_s;i++){
      for(int m=0;m<median_s;m++){
	  bcopy(fel[m]+i,med_vec+m*median_s,median_s*sizeof(float));
      }
      output->fel[median_s2][i+median_s2]=median(med_vec,med_size);
      for(int j = median_s ;j<y_size;j++){
	
 	  bcopy(fel[j]+i,(med_vec+(j%median_s)*median_s),median_s*sizeof(float));
 	  bcopy(med_vec,med_vec_tmp,med_size*sizeof(float));
          output->fel[j-median_s2][i+median_s2]=median(med_vec_tmp,med_size);	 	
         // if(!(i%300))cout << j << " " << i << " " << output->fel[j-median_s2][i+median_s2] << endl;
      }
    }
    delete []med_vec;
   delete []med_vec_tmp;
 
}

void ImageContent::dilation(ImageContent *output, int elem_rad, float val)
{
    if(getType()!=UCHAR1FLOAT1 && getType()!=UCHAR3FLOAT1){
      cout<< "ImageContent::median2d: wrong input type "<< endl; 
      return;
    }
//    cout << "ImageContent::erosion " << elem_rad << " " << x_size<<" " << y_size<<  endl;
  
   float refsum = val*(2*elem_rad+1)*(2*elem_rad+1);
   //cout << "ImageContent::erosion " << elem_rad << " " << x_size<<" " << y_size<<  " " << refsum<< endl;
   //float rcount=(2*elem_rad+1)*(2*elem_rad+1);
   float sum,count;
    for(uint j=elem_rad;j<y_size-elem_rad;j++){
      for(uint i=elem_rad;i<x_size-elem_rad;i++){
	if(fel[j][i]==val){
        for(uint m=j-elem_rad;m<=j+elem_rad;m++){
            for(uint n=i-elem_rad;n<=i+elem_rad;n++){
	      output->fel[m][n]=val;
	    }
	 }
	  
	}	
      }      
    } 
}

void ImageContent::erosion(ImageContent *output, int elem_rad, float val)
{
    if(getType()!=UCHAR1FLOAT1 && getType()!=UCHAR3FLOAT1){
      cout<< "ImageContent::median2d: wrong input type "<< endl; 
      return;
    }
//    cout << "ImageContent::erosion " << elem_rad << " " << x_size<<" " << y_size<<  endl;
  
   float refsum = val*(2*elem_rad+1)*(2*elem_rad+1);
   //cout << "ImageContent::erosion " << elem_rad << " " << x_size<<" " << y_size<<  " " << refsum<< endl;
   //float rcount=(2*elem_rad+1)*(2*elem_rad+1);
   float sum,count;
    for(uint j=elem_rad;j<y_size-elem_rad;j++){
      for(uint i=elem_rad;i<x_size-elem_rad;i++){
	sum=0;//count=0;
         for(uint m=j-elem_rad;m<=j+elem_rad;m++){
            for(uint n=i-elem_rad;n<=i+elem_rad;n++){
	      sum+=fel[m][n];
	      //count++;
	    }
	 }
	 //cout << j << " " << i << " sum " << sum << " ref " << refsum<< " " << count << " "<< rcount << endl;
	 if(sum<refsum)output->fel[j][i]=0;
	 else {
	 //  cout << j << " " << i << " sum " << sum << endl;
	   output->fel[j][i]=val;
	 }
	
      }      
    }
   
}

void ImageContent::set(ImageContent *r, ImageContent *g, ImageContent *b)
{
  if(r->size()!=tsize || g->size()!=tsize || b->size()!=tsize || buftype!=UCHAR3){
    cout << tsize << " " << r->size()<< endl;
    cout << "ERROR ImageContent::set different image size"<< endl;exit(1);  
  }
  if(r->getType()==UCHAR1 && g->getType()==UCHAR1 &&  b->getType()==UCHAR1){
    memcpy(belr[0],r->bel[0],tsize*sizeof(unsigned char)); 
    memcpy(belg[0],g->bel[0],tsize*sizeof(unsigned char)); 
    memcpy(belb[0],b->bel[0],tsize*sizeof(unsigned char)); 
  }else if((r->getType()==FLOAT1 && g->getType()==FLOAT1 && b->getType()==FLOAT1) ||
            (r->getType()==UCHAR1FLOAT1 && g->getType()==UCHAR1FLOAT1 && b->getType()==UCHAR1FLOAT1) ||
            (r->getType()==UCHAR3FLOAT1 && g->getType()==UCHAR3FLOAT1 && b->getType()==UCHAR3FLOAT1)){
    for (uint col = 0; col < tsize; col++){
      belr[0][col] = (unsigned char)r->fel[0][col];
      belg[0][col] = (unsigned char)g->fel[0][col];
      belb[0][col] = (unsigned char)b->fel[0][col];
    }  
  }else{
    
    cout << "ERROR ImageContent::set format problem "<< endl;exit(0); 
  }
}

void ImageContent::get(ImageContent *r, ImageContent *g, ImageContent *b)
{ 
  if(r->size()!=tsize || g->size()!=tsize || b->size()!=tsize || r->getType()!=FLOAT1 || g->getType()!=FLOAT1 || b->getType()!=FLOAT1){
    cout << "error:put different image size or wrong type "<< tsize << "  "  << g->size() << endl;return;  
  }
  if(buftype==UCHAR3 )
    for (uint col = 0; col < tsize; col++){
      r->fel[0][col] = (float)belr[0][col];
      g->fel[0][col] = (float)belg[0][col];
      b->fel[0][col] = (float)belb[0][col];
    }
    else  if(buftype==UCHAR1 ){
      for (uint col = 0; col < tsize; col++){
        r->fel[0][col] = (float)bel[0][col];
        g->fel[0][col] = (float)bel[0][col];
        b->fel[0][col] = (float)bel[0][col];
      } 
    }
}

/****************************************************************/
ImageContent::ImageContent(const char *name)
{
  int type;

   ConImage im;
   size_t w, h, c; 
   int ret = read_image(name, im.data, w, h, c);
   if (ret)
   {
      printf("Error reading image %s, ", name);
      switch(ret)
      {
      case -1:
         printf("cannot open file.\n"); break;
      case -2:
         printf("invalid file or unknown format.\n"); break;
      case -3:
         printf("unsupported format.\n"); break;
      case -4:
         printf("corrupted file?\n"); break;
      }
      exit(1);
   }
    y_size=h;x_size=w;
    if(c==1){
      init(y_size,x_size, UCHAR1);buftype=UCHAR1;
      uint k=0;
      for(uint row = 0; row < y_size; row++) {
	for(uint col = 0; col < x_size; col++) {
	  bel[row][col] = im.data[k++];
	}
      }
       
    coltype=GRAY;
    }else if(c==3){
      init(y_size,x_size,UCHAR3);buftype=UCHAR3;
      uint k=0;
       for(uint row = 0; row < y_size; row++) {
	 for(uint col = 0;  col < x_size; col++) {
	  belr[row][col] = im.data[k++] ;
	  belg[row][col] = im.data[k++] ;
	  belb[row][col] = im.data[k++] ;	    
	}  
       }
       coltype=RGB;
    
    }
//cout << "IMAGE READ "<< name << endl;
delete [] im.data;
return;

  type =  verif_format(name);
  //  cout << type << endl;
  if(type ==0){
      cerr  << "unrecognize image format: " << name << endl;
  }else if(type==4){

    png_structp png_ptr;
    png_infop info_ptr;
    
    FILE *fp = fopen(name, "rb");
    if(fp == NULL)return;

    unsigned char buf[PNG_BYTES_TO_CHECK];
    if ((int)fread(buf, 1, PNG_BYTES_TO_CHECK, fp) != PNG_BYTES_TO_CHECK)
      return;
    if(png_sig_cmp(&buf[0], (png_size_t)0, PNG_BYTES_TO_CHECK))return;

    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (png_ptr == NULL) {
      fclose(fp);
      return;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == NULL) {
      fclose(fp);
      png_destroy_read_struct(&png_ptr, NULL, NULL);
      return;
    }    

    if (setjmp(png_jmpbuf(png_ptr))) {
      png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
      fclose(fp);
      return;
    }    
    png_init_io(png_ptr, fp);
    png_set_sig_bytes(png_ptr, PNG_BYTES_TO_CHECK);

    int bit_depth, pix_channels, color_type;
    png_read_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
    x_size = (uint)png_get_image_width(png_ptr, info_ptr);
    y_size = (uint)png_get_image_height(png_ptr, info_ptr);
    bit_depth = png_get_bit_depth(png_ptr, info_ptr);
    pix_channels = png_get_channels( png_ptr,  info_ptr);
    color_type=png_get_color_type(png_ptr, info_ptr);

    if (color_type == PNG_COLOR_TYPE_PALETTE)
      png_set_palette_to_rgb(png_ptr);
    //if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) 
     // png_set_gray_1_2_4_to_8(png_ptr);
    if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS)) 
      png_set_tRNS_to_alpha(png_ptr);
    
    fclose(fp);
    
    if(color_type!=PNG_COLOR_TYPE_GRAY && color_type!=PNG_COLOR_TYPE_RGB){
      cout << "Attention! strange color type "<< color_type << endl;
      
    }
  
    png_bytep *src_rows;
    src_rows = png_get_rows(png_ptr, info_ptr);
    
    if(pix_channels==1){
      init(y_size,x_size, UCHAR1);buftype=UCHAR1;
      for(uint row = 0; row < y_size; row++) {
	for(uint col = 0; col < x_size; col++) {
	  bel[row][col] = src_rows[row][col];
	}
      }
       
      coltype=GRAY;
    }else if(pix_channels==3){
      init(y_size,x_size,UCHAR3);buftype=UCHAR3;
       for(uint row = 0; row < y_size; row++) {
	 for(uint col = 0, k=0; col < x_size; col++) {
	  belr[row][col] = src_rows[row][k++] ;
	  belg[row][col] = src_rows[row][k++] ;
	  belb[row][col] = src_rows[row][k++] ;	    
	}  
       }
       coltype=RGB;
    
    }
 
    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
    
  }else{//PPM image
    ifstream image(name);
	char *buf1;
	char* bigbuf;

	buf1 = new char[4];
	bigbuf = new char [512];
	image >> buf1;
	image >> buf1;
	while( buf1[0] == '#') {  image.getline(bigbuf,512);  image >> buf1; } // skipping the comments

	 x_size = atoi(buf1);
	 image >>  buf1; 
	 y_size = atoi(buf1);

	image >> buf1; // number of graylevels
	image.getline(bigbuf,512);
	
	if(type ==1 || type ==2){
	  init(y_size,x_size,UCHAR1);	  
	  image.read((char*)bel[0],tsize);
	  buftype=UCHAR1;
	}else if(type == 3){
	  init(y_size,x_size,UCHAR3);
	  unsigned char *buff = new unsigned char [3*tsize];	  
	  image.read((char*)buff,3*tsize);
	  for (uint k = 0, i = 0 ; i <tsize ; i++){
	    belr[0][i] = buff[k++] ;
	    belg[0][i] = buff[k++] ;
	    belb[0][i] = buff[k++] ;	    
	  }
	    
	  buftype=UCHAR3;
	  delete [] buff;
	}
	delete [] buf1;
	delete [] bigbuf;
	    
	image.close();
  }
}

/****************************************************************/
// Ecriture d'une Image 
void ImageContent::write(const char* name){
    write(name,"# comments");
}


void ImageContent::writePNG(const char* name, unsigned char **buff){
  if(buftype==0)return;
  //cout << "writing " << name << endl;
  png_structp png_ptr;
  png_infop info_ptr;
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (png_ptr == NULL){cout << "png_create_write_struct error...  " << endl;exit(0);}
  info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL) {
    png_destroy_write_struct(&png_ptr, NULL);
    cout << "png_create_info_struct error... " << endl;exit(0);
  }
  if (setjmp(png_jmpbuf(png_ptr))) {
    png_destroy_write_struct(&png_ptr, &info_ptr);
    cout << "setjmp error...  " << endl;exit(0);
  }
  
  FILE *fp = fopen(name, "wb");
  if(fp == NULL)return;
  png_init_io(png_ptr, fp);
  int bit_depth=8;
  
  png_set_IHDR(png_ptr, info_ptr, x_size, y_size, bit_depth, PNG_COLOR_TYPE_GRAY,
		 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  png_set_rows(png_ptr, info_ptr, buff);
    //png_write_info(png_ptr, info_ptr);
    //png_write_image(png_ptr, bel);
    //png_write_end(png_ptr, info_ptr);
  png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
    //cout << "OK2 wrote " << buftype << endl;

  png_destroy_write_struct(&png_ptr, &info_ptr);
  fclose(fp);
}



void ImageContent::writePNG(const char* name){
  if(buftype==0)return;
  
  png_structp png_ptr;
  png_infop info_ptr;
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (png_ptr == NULL){cout << "png_create_write_struct error...  " << endl;exit(0);}
  info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL) {
    png_destroy_write_struct(&png_ptr, NULL);
    cout << "png_create_info_struct error...  " << endl;exit(0);
  }
  if (setjmp(png_jmpbuf(png_ptr))) {
    png_destroy_write_struct(&png_ptr, &info_ptr);
    cout << "setjmp error...  " << endl;exit(0);
  }
  
  FILE *fp = fopen(name, "wb");
  if(fp == NULL)return;
  png_init_io(png_ptr, fp);
  int bit_depth=8;
  
  if(buftype==UCHAR1 || buftype==UCHAR1FLOAT1 || buftype==FLOAT1){
    if(buftype==UCHAR1FLOAT1 || buftype==FLOAT1)convert(UCHAR1FLOAT1);
    png_set_IHDR(png_ptr, info_ptr, x_size, y_size, bit_depth, PNG_COLOR_TYPE_GRAY,
		 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_set_rows(png_ptr, info_ptr, bel);
    //png_write_info(png_ptr, info_ptr);
    //png_write_image(png_ptr, bel);
    //png_write_end(png_ptr, info_ptr);
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
    //cout << "OK2 wrote " << buftype << endl;
  
  }
  if(buftype==UCHAR3 || buftype==UCHAR3FLOAT1){
    png_set_IHDR(png_ptr, info_ptr, x_size, y_size, bit_depth, PNG_COLOR_TYPE_RGB,
		 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    
    
    unsigned char **bell = new unsigned char*[y_size];
    bell[0]=new unsigned char[3*tsize];
    for(uint i=1;i<y_size;i++)bell[i]=bell[0]+i*3*x_size;  	
 
    for (uint row = 0; row <y_size ; row++){
      for (uint col = 0, k = 0 ; col <x_size ; col++){
	bell[row][k++]=belr[row][col];
	bell[row][k++]=belg[row][col];
	bell[row][k++]=belb[row][col];
      }
    }
    png_set_rows(png_ptr, info_ptr, bell);
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
    delete [] bell[0];delete []  bell;bell=NULL;
  }
  png_destroy_write_struct(&png_ptr, &info_ptr);
  fclose(fp);
}


void ImageContent::write(const char* name,const char* comments)
{ 
  if(buftype==0)return;
  
  if(buftype==FLOAT1 || buftype==UCHAR1FLOAT1){
    convert(UCHAR1FLOAT1);
  }
  
  if(buftype==UCHAR1 || buftype==UCHAR1FLOAT1){    
    writePGM(name, bel[0], comments);
    if(buftype==UCHAR1FLOAT1){
      delete[] bel[0];delete[] bel;
      buftype=FLOAT1;
    }
  }
  
  if(buftype==UCHAR3 || buftype==UCHAR3FLOAT1){
    if(buftype==UCHAR3FLOAT1)convert(UCHAR1FLOAT1);
    ofstream theexit (name);
    if( !theexit)  { cerr << "Enable to open " << name << endl;   return;}
    
    theexit <<"P6"<<endl; 
    theexit << "# "<< comments <<endl;
    theexit << x_size << " " << y_size << endl;
    theexit << "255" << endl;
    
    unsigned char *buff = new unsigned char [3*tsize];
    for (uint k = 0, i = 0 ; i <tsize ; i++){
      buff[k++]=belr[0][i];
      buff[k++]=belg[0][i];
      buff[k++]=belb[0][i];
    }
    theexit.write((char *)buff,3*tsize);
    theexit.close();
    delete []buff;
    if(buftype==UCHAR1FLOAT1){
      delete [] belr[0];delete []  belr;belr=NULL;
      delete [] belg[0];delete []  belg;belg=NULL;
      delete [] belb[0];delete []  belb;belb=NULL;
      buftype=FLOAT3;
    }
  }
}

void ImageContent::writeF(const char* name){
    unsigned char **bell = new unsigned char*[y_size];
    bell[0]=new unsigned char[tsize];
    for(uint i=1;i<y_size;i++)bell[i]=bell[0]+i*x_size;  
    for(uint i=1;i<tsize;i++)bell[0][i]=(unsigned char)fel[0][i];
    writePNG(name,bell);
   delete [] bell[0];delete []  bell;bell=NULL;

}
void ImageContent::writeR(const char* name){
  writePNG(name,belr);
}
void ImageContent::writeG(const char* name){
  writePNG(name,belg);
}
void ImageContent::writeB(const char* name){
  writePNG(name,belb);
}


void ImageContent::writePGM(const char* name, unsigned char* buff, const char* comments){
  if(buftype==0)return;
  
    ofstream theexit(name);
    if( !theexit)  { cerr << "Enable to open " << name << endl;   return;}
    
    theexit <<"P5"<<endl; 
    theexit << "# " << comments <<endl;
    theexit << x_size << " " << y_size << endl;
    theexit << "255" << endl; 
    
    theexit.write((char*)buff,tsize);
    theexit.close();

}

void ImageContent::RGB2xyY(){
  if(buftype!=UCHAR3){
    cout << "ERROR ImageContent::RGB2xyYwrong buftype - convert to float"<< endl;exit(1);
  }
  for(uint i=0;i<tsize;i++){
     felr[0][i]/=255.0;//norm(inverse_gamma_correction(rgb->fel3[j][i]));
     felg[0][i]/=255.0;//norm(inverse_gamma_correction(rgb->fel2[j][i]));
     felb[0][i]/=255.0;//norm(inverse_gamma_correction(rgb->fel1[j][i]));  
  }
  
  float r,g,b;
  //float sum;
  for (uint i = 0 ; i < tsize ; i++){
      r=((float)belr[0][i])/255.0;
      g=((float)belg[0][i])/255.0;
      b=((float)belb[0][i])/255.0;
      belr[0][i] =(unsigned char)(255.0*(0.412453 * r + 0.357580 * g + 0.180423 * b)) ;
      belg[0][i] =(unsigned char)(255.0*(0.212671 * r + 0.715160 * g + 0.072169 * b)) ;
      belb[0][i] =(unsigned char)(255.0*(0.019334 * r + 0.119193 * g + 0.950227 * b)) ;
 /*     sum=r+g+b;
      if(sum!=0){
	belr[0][i] = sum ;
	belg[0][i] = sum ;
      }*/
  }
}


void ImageContent::RGB2opp(){
  if(buftype!=UCHAR3 && buftype!=UCHAR3FLOAT1){
    cout << "ERROR ImageContent::RGB2opp wrong buftype - convert to uchar3"<< endl;
    exit(1);
  }
  
  float r,g,b;
  for (uint i = 0 ; i < tsize ; i++){
    r=belr[0][i];
    g=belg[0][i];
    b=belb[0][i];
    belr[0][i] =(unsigned char)(0.40 * r + 0.35 * g + 0.25 * b) ;
    belg[0][i] =118+(unsigned char)(0.53 * r - 0.46 * g);//(0.62 * r - 0.54 * g);
    belb[0][i] =175+(unsigned char)((0.23 * r + 0.22 * g - b)/1.46) ;
  }
  coltype=OPP;
}


void ImageContent::RGB2lbrg(){
  if(buftype!=UCHAR3 && buftype!=UCHAR3FLOAT1){
    cout << "ERROR ImageContent::RGB2rgb wrong buftype - convert to uchar3"<< endl;
    exit(1);
  }
  float r,g,b,lg,lr,lb;
  for(uint i=0;i<x_size;i++){
    for(uint j=0;j<y_size;j++){
      r=(float)belr[j][i];
      g=(float)belg[j][i];
      b=(float)belb[j][i];
      lr=(g==0)?0:log(r);
      lg=(g==0)?0:log(g);
      lb=(b==0)?0:log(b);
      belb[j][i]=(unsigned char)(lg/5.55);
      belg[j][i]=(unsigned char)((lr-lg+5.55)/11.1);
      belr[j][i]=(unsigned char)(((lb-(lg+lr)/2.0)+5.55)/11.1);	    
    }
  }
}
void ImageContent::RGB2rgb(){
  if(buftype!=UCHAR3){
    cout << "ERROR ImageContent::RGB2rgb wrong buftype - convert to float"<< endl;
    exit(1);
  }
  float sum;
  for(uint i=0;i<x_size;i++){
    for(uint j=0;j<y_size;j++){
      sum=felr[j][i]+felg[j][i]+felb[j][i];
      if(sum!=0){
	felb[j][i]/=sum;//norm(inverse_gamma_correction(rgb->fel3[j][i]));
	felg[j][i]/=sum;//norm(inverse_gamma_correction(rgb->fel2[j][i]));
	felr[j][i]=sum;//norm(inverse_gamma_correction(rgb->fel1[j][i]));  
      }
    }
  }
}


/****************************************************************/
void ImageContent::interpolate(DARY *im_in, float m_x, float m_y, 
				 float scalex, float scaley, 
				 float angle){
    float lecos = cos(2*M_PI*angle/360);
    float lesin = sin(2*M_PI*angle/360);	        
    float vec0x = (1.0/scalex)*lecos; 
    float vec0y = (1.0/scaley)*lesin;
    float vec1x = -lesin*(1.0/scalex);
    float vec1y = lecos*(1.0/scaley);
    interpolate(im_in, m_x, m_y, vec0x , vec0y  ,vec1x ,vec1y);
}

/****************************************************************/
void ImageContent::interpolate(DARY *im_in, float m_x, float m_y, float vec0x, float vec0y,
			       float vec1x, float vec1y){
    if((buftype!=FLOAT1 && buftype!=UCHAR1FLOAT1 && buftype!=UCHAR3FLOAT1) || 
            (im_in->getType()!=FLOAT1 && im_in->getType()!=UCHAR1FLOAT1 && im_in->getType()!=UCHAR3FLOAT1)){
        cout << "ERROR ImageContent::interpolate wrong buftype " << buftype<< endl;exit(1);
    }
    float px, py, dx, dy;
    int arondx, arondy;
    int width_2  = x_size>>1;
    int height_2 = y_size>>1;
    if(!(x_size%2))width_2--;
    if(!(y_size%2))height_2--;
    int xim_in = im_in->x()-1;
    int yim_in = im_in->y()-1;
    for(int j=-height_2;j<=height_2; j++){
      for(int i=-width_2;i<=width_2; i++){
	px = m_x + i*vec0x + j*vec0y;
	py = m_y + i*vec1x + j*vec1y;
	    arondx = (int) floor(px);
	    arondy = (int) floor(py);
	    dx = px -  arondx;
	    dy = py -  arondy;     
	    if(arondx>=0 && arondy>=0 && arondx<xim_in && arondy<yim_in){
		fel[j+height_2][i+width_2] = 
		    ((1.0 - dy)*((1.0 - dx)*im_in->fel[arondy][arondx] +
				 dx*im_in->fel[arondy][arondx+1]) +
		     dy*((1.0 - dx)*im_in->fel[arondy+1][arondx] +
			 dx *im_in->fel[arondy+1][arondx+1]));
	    }else if(arondx<0 && arondy<0){
	      fel[j+height_2][i+width_2]= im_in->fel[0][0];
	    } else if(arondx>=xim_in && arondy>=yim_in){
	      fel[j+height_2][i+width_2]= im_in->fel[yim_in][xim_in];
	    } else if(arondx<0 && arondy>yim_in){
	      fel[j+height_2][i+width_2]= im_in->fel[yim_in][0];
	    } else if(arondx>=xim_in && arondy<0){
	      fel[j+height_2][i+width_2]= im_in->fel[0][xim_in];
	    } else if(arondx>xim_in && arondy>=0){
	      fel[j+height_2][i+width_2]= im_in->fel[arondy][xim_in];
	    } else if(arondx>=0 && arondy<0){
	      fel[j+height_2][i+width_2]= im_in->fel[0][arondx];
	    } else if(arondx>=0 && arondy>=yim_in){
	      fel[j+height_2][i+width_2]= im_in->fel[yim_in][arondx];
	    } else if(arondx<0 && arondy>=0){
	      fel[j+height_2][i+width_2]= im_in->fel[arondy][0];
	    } else if(arondx>=xim_in && arondy>=0){
	      fel[j+height_2][i+width_2]= im_in->fel[arondy][xim_in];
	      }
	}
    }  
}

/****************************************************************/
void ImageContent::scale_float(float**im_in, int xim_in, int yim_in, float**im_out, float scalex, float scaley){
  float px, py, dx, dy;
    int arondx, arondy;
    //cout << "scaling " << endl;
    yim_in=yim_in-1;
    xim_in=xim_in-1;
    for(uint j=0;j<y_size; j++){
	for(uint i=0;i<x_size; i++){
	  px = i*scalex; 
	  py = j*scaley;
	  arondx = (int) floor(px);
	  arondy = (int) floor(py);
	  dx = px -  arondx;
	  dy = py -  arondy;     
          //cout << i <<"  "<< j << " " << x_size<< " " <<  y_size << "  "<< arondx << "  "<< arondy<< "  "<< xim_in << "  "<< yim_in << " " << im_out[j][i]<< endl;
          //im_out[j][i]=1.0;
	  if(arondx<xim_in && arondy<yim_in){
            im_out[j][i] = (1.0 - dy)*((1.0 - dx)*im_in[arondy][arondx] +
				  dx*im_in[arondy][arondx+1]) +
	      dy*((1.0 - dx)*im_in[arondy+1][arondx] +
		  dx *im_in[arondy+1][arondx+1]);
              // cout << a << endl;
	  }
          else if(arondx==xim_in && arondy<yim_in){
            im_out[j][i] = (1.0 - dy)*im_in[arondy][arondx] +
	      dy*(im_in[arondy+1][arondx]);	  	      
	  } else if(arondx<xim_in && arondy==yim_in){
            im_out[j][i]= (1.0 - dx)*im_in[arondy][arondx] +
	      dx*im_in[arondy][arondx+1];	      
	  }else if(arondx==xim_in && arondy==yim_in){
            im_out[j][i]= im_in[arondy][arondx];
          }else if(arondx>xim_in && arondy<=yim_in){im_out[j][i]=im_in[arondy][xim_in];
          }else if(arondx<=xim_in && arondy>yim_in){im_out[j][i]=im_in[yim_in][arondx];
          }else im_out[j][i]=im_in[yim_in][xim_in]; 
	}
    } //cout << "scaling done" << endl;
}
/****************************************************************/
void ImageContent::scale_uchar(unsigned char**im_in, int xim_in, int yim_in, unsigned char**im_out, float scalex, float scaley){
      float px, py, dx, dy;
      int arondx, arondy;
      yim_in=yim_in-1;
      xim_in=xim_in-1;   
      for(uint j=0;j<y_size; j++){
        for(uint i=0;i<x_size; i++){
          px = i*scalex; 
          py = j*scaley;
          arondx = (int) floor(px);
          arondy = (int) floor(py);
          dx = px -  arondx;
          dy = py -  arondy;     
          if(arondx<xim_in && arondy<yim_in){
            im_out[j][i]=(unsigned char) ((1.0 - dy)*((1.0 - dx)*im_in[arondy][arondx] +
                dx*im_in[arondy][arondx+1]) +
                dy*((1.0 - dx)*im_in[arondy+1][arondx] +
                dx *im_in[arondy+1][arondx+1]));
          } else if(arondx==xim_in && arondy<yim_in){
            im_out[j][i] = (unsigned char) ((1.0 - dy)*im_in[arondy][arondx] +
                dy*(im_in[arondy+1][arondx]));	  	      
          } else if(arondx<xim_in && arondy==yim_in){
            im_out[j][i]= (unsigned char) ((1.0 - dx)*im_in[arondy][arondx] +
                dx*im_in[arondy][arondx+1]);	      
          }else if(arondx==xim_in && arondy==yim_in){
            im_out[j][i]= im_in[arondy][arondx];
          }else if(arondx>xim_in && arondy<=yim_in){im_out[j][i]=im_in[arondy][xim_in];
          }else if(arondx<=xim_in && arondy>yim_in){im_out[j][i]=im_in[yim_in][arondx];
          }else im_out[j][i]=im_in[yim_in][xim_in]; 
        }
      } 
}

// void ImageContent::decrease(uchar *, int rescWidth, int rescHeight){
// }
void ImageContent::decrease(DARY *im_in){
    if((im_in->getType()==UCHAR3 || im_in->getType()==UCHAR3FLOAT1) && (buftype==UCHAR3 || buftype==UCHAR3FLOAT1)){
     decrease(im_in->belr, im_in->x(), im_in->y(), belr);
     decrease(im_in->belg, im_in->x(), im_in->y(), belg);
     decrease(im_in->belb, im_in->x(), im_in->y(), belb);
   }
}


void ImageContent::decrease(unsigned char **ucharArray, int xim_in, int yim_in, unsigned char**im_out){
  int tmpSize = x_size*yim_in;
  float* tmpArray = new float[tmpSize];
  if (x_size < xim_in) {//decrease nb of columns
    bzero(tmpArray,tmpSize*sizeof(float));
    float scale = ((float)xim_in)/(float)x_size;
    for (int y = 0; y < yim_in; y++) {
      int framePointX = y*xim_in;
      int rescPointX = y*x_size;
      int frameRowEnd = framePointX+xim_in;
      int rescRowEnd = rescPointX+x_size;
      float tweight = scale;// scale > 1, tweight+pweight=scale
      float pweight = 1.0;
      do {
        if (tweight > 1.0) {
          tmpArray[rescPointX] += pweight*(float)ucharArray[0][framePointX];
	  framePointX++;
          tweight -= pweight;
          pweight = 1.0;
          if (tweight <= 0.0) {
            tweight = scale;
            rescPointX++;//start to sum up in row for new point with pweight
          }
        }else {
          tmpArray[rescPointX] += tweight*(float)ucharArray[0][framePointX];
	  rescPointX++;
          pweight = 1.0-tweight;//start to sum up in row for new point with pweight
          tweight = scale;
        }
      }
      while (framePointX < frameRowEnd && rescPointX < rescRowEnd);
    }
  } else {// it seems x_size >= frameWidth, so it must be y_size <= frameHeight

    //float* tmp = tmpArray;
    //tmpArray = ucharArray;
    //floatArray = tmp;
  }

  //delete[] floatArray; 
  int aSize = x_size*y_size;
  float *floatArray = new float[aSize];
  if (y_size < yim_in) {
    bzero(floatArray,aSize*sizeof(float));
    float scale = ((float)yim_in)/(float)y_size;
    for (int x = 0; x < x_size; x++) { //decrease nb of rows
      int framePointY = x;
      int rescPointY = x;
      int frameColEnd = yim_in*x_size+x;
      int rescColEnd = y_size*x_size+x;
      float tweight = scale;
      float pweight = 1.0;
      do {
        if (tweight > 1.0) {
          floatArray[rescPointY] += pweight*tmpArray[framePointY];
          framePointY += x_size;
          tweight -= pweight;
          pweight = 1.0;
         if (tweight <= 0.0) {
            rescPointY += x_size;//start to sum up in col for new point with pweight
            tweight = scale;
          }
        }
        else {
          floatArray[rescPointY] += tweight*tmpArray[framePointY];
          rescPointY += x_size;
          pweight = 1.0-tweight;
          tweight = scale;
        }
      }
      while (framePointY < frameColEnd && rescPointY < rescColEnd);
    }
  }
  else {
    //float* tmp = floatArray;
    //floatArray = tmpArray;
    //tmpArray = tmp;
  }
  // Normalize by the number of pixels summed up for each output.
  float norm = ((float)aSize)/(yim_in*xim_in);
  for (int i = 0; i < aSize; i++)im_out[0][i]=(uchar)(floatArray[i] * norm);

  //=x_size*y_size;
   delete[] tmpArray;
   delete[] floatArray;
}



/****************************************************************/
void ImageContent::scale(DARY *im_in, float scalex, float scaley){
  if((im_in->getType()==FLOAT1 || im_in->getType()==UCHAR1FLOAT1 || im_in->getType()==UCHAR3FLOAT1) && 
      (buftype==UCHAR1FLOAT1 || buftype==FLOAT1 || buftype==UCHAR3FLOAT1)){
    scale_float(im_in->fel, im_in->x(), im_in->y(), fel, scalex, scaley);
  }
  if((im_in->getType()==UCHAR1 || im_in->getType()==UCHAR1FLOAT1) && (buftype==UCHAR1 || buftype==UCHAR1FLOAT1)){
    scale_uchar(im_in->bel, im_in->x(), im_in->y(), bel, scalex, scaley);
  }
  if((im_in->getType()==UCHAR3 || im_in->getType()==UCHAR3FLOAT1) && (buftype==UCHAR3 || buftype==UCHAR3FLOAT1)){
    scale_uchar(im_in->belr, im_in->x(), im_in->y(), belr, scalex, scaley);
    scale_uchar(im_in->belg, im_in->x(), im_in->y(), belg, scalex, scaley);
    scale_uchar(im_in->belb, im_in->x(), im_in->y(), belb, scalex, scaley);
  }
    
}




void ImageContent::normalize(float min_in, float max_in){
  float sc=255.0/(max_in-min_in);
  if(min_in<0 || max_in >255){
    cout << " min " << min_in << " < 0 or max "<< max_in << " > 255, sc=" << sc<< endl;
  }
  if(buftype==FLOAT1 || buftype==UCHAR1FLOAT1 || buftype==UCHAR3FLOAT1){
    for (uint i = 0; i < tsize; i++){
      if(fel[0][i]<min_in)fel[0][i]=min_in;
      else if(fel[0][i]>max_in)fel[0][i]=max_in;
      fel[0][i] = ((fel[0][i]-min_in)*sc);      
    }   
  } 
  if(buftype==UCHAR1 || buftype==UCHAR1FLOAT1){
    for (uint col = 0; col < tsize; col++){
      if(bel[0][col]<min_in)bel[0][col]=(unsigned char)min_in;
      else if(bel[0][col]>max_in)bel[0][col]=(unsigned char)max_in;
      bel[0][col] = (unsigned char)((bel[0][col]-min_in)*sc);
    }    
  }
  if(buftype==UCHAR3 || buftype==UCHAR3FLOAT3 || buftype==UCHAR3FLOAT1){
    for (uint col = 0; col < tsize; col++){
      if(belr[0][col]<min_in)belr[0][col]=(unsigned char)min_in;
      else if(belr[0][col]>max_in)belr[0][col]=(unsigned char)max_in;
      belr[0][col] = (unsigned char)((belr[0][col]-min_in)*sc);
      if(belg[0][col]<min_in)belg[0][col]=(unsigned char)min_in;
      else if(belg[0][col]>max_in)belg[0][col]=(unsigned char)max_in;
      belg[0][col] = (unsigned char)((belg[0][col]-min_in)*sc);
      if(belb[0][col]<min_in)belb[0][col]=(unsigned char)min_in;
      else if(belb[0][col]>max_in)belb[0][col]=(unsigned char)max_in;
      belb[0][col] = (unsigned char)((belb[0][col]-min_in)*sc);
    }    
  }     
  if(buftype==FLOAT3 || buftype==UCHAR3FLOAT3){
    for (uint col = 0; col < tsize; col++){
      if(felr[0][col]<min_in)felr[0][col]=(unsigned char)min_in;
      else if(felr[0][col]>max_in)felr[0][col]=(unsigned char)max_in;
      felr[0][col] = (unsigned char)((felr[0][col]-min_in)*sc);
      if(felg[0][col]<min_in)felg[0][col]=(unsigned char)min_in;
      else if(felg[0][col]>max_in)felg[0][col]=(unsigned char)max_in;
      felg[0][col] = (unsigned char)((felg[0][col]-min_in)*sc);
      if(felb[0][col]<min_in)felb[0][col]=(unsigned char)min_in;
      else if(felb[0][col]>max_in)felb[0][col]=(unsigned char)max_in;
      felb[0][col] = (unsigned char)((felb[0][col]-min_in)*sc);
    }    
  }     
}

void ImageContent::normalize(){

  float max,min;
  if(buftype==FLOAT1 || buftype==UCHAR1FLOAT1 || buftype==UCHAR3FLOAT1){
    max=fel[0][0];min=fel[0][0];
    for (uint col = 0; col < tsize; col++){
      if(max<fel[0][col])max=fel[0][col];
      else if(min>fel[0][col])min=fel[0][col];
    }
    normalize( min,  max);
  }else if(buftype==UCHAR1){
    max=bel[0][0];min=bel[0][0];
    for (uint col = 0; col < tsize; col++){
      if(max<(float)bel[0][col])max=(float)bel[0][col];
      else if(min>(float)bel[0][col])min=(float)bel[0][col];
    }
    normalize( min,  max);
  }else if(buftype==UCHAR3){
    max=bel[0][0];min=bel[0][0];
    for (uint col = 0; col < tsize; col++){
      if(max<(float)belr[0][col])max=(float)belr[0][col];
      else if(min>(float)belr[0][col])min=(float)belr[0][col];
      if(max<(float)belg[0][col])max=(float)belg[0][col];
      else if(min>(float)belg[0][col])min=(float)belg[0][col];
      if(max<(float)belb[0][col])max=(float)belb[0][col];
      else if(min>(float)belb[0][col])min=(float)belb[0][col];
    }
    normalize( min,  max);
  }else{
    cout <<  "ERROR ImageContent::normalize no such buffer type " << buftype << endl;exit(0);          
  }
}       
        
void ImageContent::crop(DARY *img, int x, int y){
   crop(img, x, y, "non");
}
void ImageContent::crop(DARY *img, int x, int y, char *mode){
  
  
  
  int dx=x_size>>1;
  int dy=y_size>>1;
  int iy=y-dy;
  int ix=x-dx;
  uint jy=y+dy;
  uint jx=x+dx;
  int sx=(ix<0)?(-ix):0;
  int sy=(iy<0)?(-iy):0;
  int ex=(jx<img->x())?(x_size):(img->x()-ix);
  int ey=(jy<img->y())?(y_size):(img->y()-iy);
   sx=0;
   sy=0;
   ex=x_size;
   ey=y_size;
   iy=y;
   ix=x;

//cout << sx << " "<< sy << " "<< ex << " "<< ey << endl;
  if(buftype==FLOAT1){
    for (int i = sy; i < ey; i++){
      for (int j = sx; j < ex; j++){	
        fel[i][j]= img->fel[iy+i][ix+j];
      }   
    }
    if(!strcmp("aver",mode)){
      int t=0,l=0,r=0,b=0;//top,left,right,bottom
      int *hist = new int[256];
      bzero(hist,256*sizeof(int));
      int v=0;
      for (int i = sy; i < ey; i++){//get hist of the first column
        v=(int)fel[i][sx];
        if(v<256 && v>=0)hist[v]++;
      }
      int maxi=0,maxv=0;
      for (int i = 0; i < 256; i++){//get most freq grayvalue
        if(maxv<hist[i]){
          l=i;maxv=hist[i];
        }    
      }
      bzero(hist,256*sizeof(int));
      for (int i = sy; i < ey; i++){//get hist of the last column
        v=(int)fel[i][ex-1];
        if(v<256 && v>=0)hist[v]++;
      }
      maxi=0,maxv=0;
      for (int i = 0; i < 256; i++){//get most freq grayvalue
        if(maxv<hist[i]){
          r=i;maxv=hist[i];
        }    
      }
      bzero(hist,256*sizeof(int));
      for (int i = sx; i < ex; i++){//get hist of the first row
        v=(int)fel[sy][i];
        if(v<256 && v>=0)hist[v]++;
      }
      maxi=0,maxv=0;
      for (int i = 0; i < 256; i++){
        if(maxv<hist[i]){
          t=i;maxv=hist[i];
        }    
      }
      bzero(hist,256*sizeof(int));
      for (int i = sx; i < ex; i++){//get hist of the last row
        v=(int)fel[ey-1][i];
        if(v<256 && v>=0)hist[v]++;
      }
      maxi=0,maxv=0;
      for (int i = 0; i < 256; i++){
        if(maxv<hist[i]){
          b=i;maxv=hist[i];
        }    
      }
      delete []hist;
      for (int j =0; j < sy; j++){
        for (int i = 0; i < x_size; i++){
          fel[j][i]=t;
        }
      }
      for (int j =ey; j < y_size; j++){
        for (int i = 0; i < x_size; i++){
          fel[j][i]=b;
        }
      }
      for (int j =sy; j < ey; j++){
        for (int i = 0; i < sx; i++){
          fel[j][i]=l;
        }
        for (int i = ex; i < x_size; i++){
          fel[j][i]=r;
        }
      }
     
    }else  if(!strcmp("mirror",mode)){
      if(ix<0){
        for (int i = sy; i < ey; i++){
          for (int j = 0; j < sx; j++){	
            fel[i][j]= fel[i][sx+sx-j];
          }		   
        }			
      }
      if(iy<0){
        for (int i = 0; i < sy; i++){
          for (int j = sx; j < ex; j++){
            fel[i][j]= fel[sy+sy-i][j];
          }		   
        }			
      }
      if(jy>=img->y()){
        for (int i = ey; i < y_size; i++){
          for (int j = sx; j < ex; j++){	
            fel[i][j]= fel[ey+ey-i-1][j];
          }		   
        }				
      }		
      if(jx>=img->x()){
        for (int i = sy; i < ey; i++){
          for (int j = ex; j < x_size; j++){	
            fel[i][j]= fel[i][ex+ex-j-1];
          }		   
        }				
      }		
      if(ix<0 && iy<0){
        for (int i = 0; i < sy; i++){
          for (int j = 0; j < sx; j++){	
            fel[i][j]= fel[i][sx+sx-j];
          }		   
        }
      }
      if(ix<0 && jy>=img->y()){
        for (int i = ey; i < y_size; i++){
          for (int j = 0; j < sx; j++){	
            fel[i][j]= fel[i][sx+sx-j];
          }		   
        }				
      }
      if(jx>=img->x() && jy>=img->y()){
        for (int i = ey; i < y_size; i++){
          for (int j = ex; j < x_size; j++){	
            fel[i][j]= fel[i][ex+ex-j-1];
          }		   
        }				
      }		
      if(jx>=img->x() && iy<0){
        for (int i = 0; i < sy; i++){
          for (int j = ex; j < x_size; j++){	
            fel[i][j]= fel[i][ex+ex-j-1];
          }		   
        }				
      }		
    }else{//copy the boundaries 
      if(ix<0){
        for (int i = sy; i < ey; i++){
          for (int j = 0; j < sx; j++){	
            fel[i][j]= fel[i][sx];
          }		   
        }			
      }
      if(iy<0){
        for (int i = 0; i < sy; i++){
          for (int j = sx; j < ex; j++){
            fel[i][j]= fel[sy][j];
          }		   
        }			
      }
      if(jy>=img->y()){
        for (int i = ey; i < y_size; i++){
          for (int j = sx; j < ex; j++){	
            fel[i][j]= fel[ey-1][j];
          }		   
        }				
      }		
      if(jx>=img->x()){
        for (int i = sy; i < ey; i++){
          for (int j = ex; j < x_size; j++){	
            fel[i][j]= fel[i][ex-1];
          }		   
        }				
      }		
      if(ix<0 && iy<0){
        for (int i = 0; i < sy; i++){
          for (int j = 0; j < sx; j++){	
            fel[i][j]= fel[sy][sx];
          }		   
        }
      }
      if(ix<0 && jy>=img->y()){
        for (int i = ey; i < y_size; i++){
          for (int j = 0; j < sx; j++){	
            fel[i][j]= fel[ey-1][sx];
          }		   
        }				
      }
      if(jx>=img->x() && jy>=img->y()){
        for (int i = ey; i < y_size; i++){
          for (int j = ex; j < x_size; j++){	
            fel[i][j]= fel[i][ex+ex-j-1];
          }		   
        }				
      }		
      if(jx>=img->x() && iy<0){
        for (int i = 0; i < sy; i++){
          for (int j = ex; j < x_size; j++){	
            fel[i][j]= fel[sy][ex-1];
          }		   
        }				
      }		
    }
  }else if(buftype==UCHAR1){
    for (int i = sy; i < ey; i++){
      for (int j = sx; j < ex; j++){
        bel[i][j]= img->bel[iy+i][ix+j];
      }   
    }    
  }else if(buftype==UCHAR3){
    for (int i = sy; i < ey; i++){
      for (int j = sx; j < ex; j++){
        belr[i][j]= img->belr[iy+i][ix+j];
        belg[i][j]= img->belg[iy+i][ix+j];
        belb[i][j]= img->belb[iy+i][ix+j];
      }   
    }    
  }  
   
}






/****************************************************************/
void ImageContent::convert(int type){
  if(type==UCHAR3 || type==UCHAR3FLOAT1){//******************************************
    init(y_size,x_size,UCHAR3);
    if(buftype==UCHAR1){//UCHAR1->UCHAR3
      memcpy(belr[0],bel[0],tsize*sizeof(unsigned char)); 
      memcpy(belg[0],bel[0],tsize*sizeof(unsigned char)); 
      memcpy(belb[0],bel[0],tsize*sizeof(unsigned char)); 
      buftype=UCHAR3;     
      if(type==UCHAR3FLOAT1){//UCHAR1->UCHAR3FLOAT1
        init(y_size,x_size,FLOAT1);
        for(uint i=0;i<tsize;i++)fel[0][i]=(float)bel[0][i];           
        buftype=UCHAR3FLOAT1;     
      }
      delete [] bel[0];delete [] bel;bel=NULL;      
    }else if(buftype==FLOAT1 || buftype==UCHAR3FLOAT1 || buftype==UCHAR1FLOAT1){
      float max_in=fel[0][0];float min_in=fel[0][0];
      for (uint col = 0; col < tsize; col++){
        if(max_in<fel[0][col])max_in=fel[0][col];
        else if(min_in>fel[0][col])min_in=fel[0][col];
      }
      float sc=255.0/(max_in-min_in);
      //for(uint i=0;i<tsize;i++)bel[0][i]=(unsigned char)((fel[0][i]-min_in)*sc);      
      unsigned char c;
      for(uint i=0;i<tsize;i++){
        c=(unsigned char)((fel[0][i]-min_in)*sc);
        belr[0][i]=c; 
        belg[0][i]=c; 
        belb[0][i]=c; 
      }
      if(buftype==UCHAR1FLOAT1){//UCHAR1FLOAT1->UCHAR3
        delete [] bel[0];delete []  bel;bel=NULL;
      }
      if(type==UCHAR3 && (buftype==UCHAR3FLOAT1 || buftype==UCHAR1FLOAT1 || FLOAT1)){//UCHAR1FLOAT1->UCHAR3, UCHAR3FLOAT1->UCHAR3, FLOAT1->UCHAR3,
        delete [] fel[0];delete []  fel;fel=NULL;
        buftype=UCHAR3;    
      }else {
        buftype=UCHAR3FLOAT1; //UCHAR1FLOAT1->UCHAR3FLOAT1, UCHAR3FLOAT1->UCHAR3FLOAT1, FLOAT1->UCHAR3FLOAT1,
      }    
    }else if(buftype==UCHAR3 && type==UCHAR3FLOAT1){//UCHAR3->UCHAR3FLOAT1
      init(y_size,x_size,FLOAT1);
      if(coltype==RGB){
        for(uint i=0;i<tsize;i++)fel[0][i]=(float)((belr[0][i]+belg[0][i]+belb[0][i])/3.0);  
      }else if(coltype==OPP){
        for(uint i=0;i<tsize;i++)fel[0][i]=(float)(belr[0][i]);              
      }
      buftype=UCHAR3FLOAT1;    
    }else{
      cout << "ERROR  ImageContent::convert"<< type<< endl;exit(1);  
    }
  }else if(type==UCHAR1FLOAT1 || type==FLOAT1){//******************************************
    init(y_size,x_size,FLOAT1);
    if(buftype==UCHAR3 || buftype==UCHAR3FLOAT1){
      for(uint i=0;i<tsize;i++)fel[0][i]=(float)((belr[0][i]+belg[0][i]+belb[0][i])/3.0);      
      delete [] belr[0];delete []  belr;belr=NULL;
      delete [] belg[0];delete []  belg;belg=NULL;
      delete [] belb[0];delete []  belb;belb=NULL;
      if(type==UCHAR1FLOAT1){ 
        init(y_size,x_size,UCHAR1);
        for(uint i=0;i<tsize;i++)bel[0][i]=(unsigned char)fel[0][i]; 
        buftype=UCHAR1FLOAT1;    //UCHAR3->UCHAR1FLOAT1, UCHAR3FLOAT1->UCHAR1FLOAT1
      }else{
        buftype=FLOAT1;//UCHAR3->FLOAT1, UCHAR3FLOAT1->FLOAT1  
      }
    }else if(buftype==UCHAR1){
      for(uint i=0;i<tsize;i++)fel[0][i]=(float)(bel[0][i]); 
      if(type==FLOAT1){
        delete [] bel[0];delete []  bel;bel=NULL;
        buftype=FLOAT1;   // UCHAR1->FLOAT1 
      }else{
        buftype=UCHAR1FLOAT1;   // UCHAR1->UCHAR1FLOAT1      
      }      
    }else if((buftype==UCHAR1FLOAT1 || buftype==FLOAT1 ) && type==UCHAR1FLOAT1){
      init(y_size,x_size,UCHAR1);
      float max_in=fel[0][0];float min_in=fel[0][0];
      for (uint col = 0; col < tsize; col++){
        if(max_in<fel[0][col])max_in=fel[0][col];
        else if(min_in>fel[0][col])min_in=fel[0][col];
      }
      float sc=255.0/(max_in-min_in);
      for(uint i=0;i<tsize;i++)bel[0][i]=(unsigned char)((fel[0][i]-min_in)*sc);      
      buftype=UCHAR1FLOAT1;   // UCHAR1FLOAT1->UCHAR1FLOAT1 ,  FLOAT1->UCHAR1FLOAT1          
    }else if(buftype==UCHAR1FLOAT1 && type==FLOAT1){
      delete [] bel[0];delete []  bel;bel=NULL;      
      buftype=FLOAT1;   // UCHAR1FLOAT1->FLOAT1            
    }else{
      cout << "ERROR  ImageContent::convert FLOAT1 "<< type<< endl;exit(1);  
    }
  }else if(type==UCHAR1){//******************************************
    if(buftype==UCHAR3){
      init(y_size,x_size,UCHAR1);
      for(uint i=0;i<tsize;i++)bel[0][i]=(float)((belr[0][i]+belg[0][i]+belb[0][i])/3.0);      
      delete [] belr[0];delete []  belr;belr=NULL;
      delete [] belg[0];delete []  belg;belg=NULL;
      delete [] belb[0];delete []  belb;belb=NULL;
      buftype=UCHAR1;   // UCHAR1FLOAT1->UCHAR1FLOAT1 ,  FLOAT1->UCHAR1FLOAT1              
    }else if(buftype==UCHAR1)return;
    else {
    cout << "ERROR  ImageContent::convert not implemented yet"<< type<< endl;exit(1);  
    }
    
  }else if(type=OPP3FLOAT1){
    if(buftype==UCHAR3FLOAT1){
       for (uint i = 0 ; i < tsize ; i++){
    	 fel[0][i]=(float)belr[0][i];
       }
    }
  }
}



/****************************************************************/
ImageContent::~ImageContent(void){
  delete []filename;
  if(buftype==UCHAR1 || buftype==UCHAR1FLOAT1){    
    delete []  bel[0];delete []  bel;bel=NULL;
  }
  if(buftype==FLOAT1 || buftype==UCHAR1FLOAT1 || buftype==UCHAR3FLOAT1){
    delete [] fel[0];delete [] fel;fel=NULL;
  }
  if(buftype==UCHAR3 || buftype==UCHAR3FLOAT3 || buftype==UCHAR3FLOAT1){
    delete [] belr[0];delete []  belr;//rbuffer=NULL;
    delete [] belg[0];delete []  belg;//gbuffer=NULL;
    delete [] belb[0];delete []  belb;//bbuffer=NULL;
  }  
  if(buftype==FLOAT3 || buftype==UCHAR3FLOAT3){
    delete [] felr[0];delete [] felr;felr=NULL;
    delete [] felg[0];delete [] felg;felg=NULL;
    delete [] felb[0];delete [] felb;felb=NULL;
  }
  buftype=0;
}
