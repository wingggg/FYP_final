
#include "util.h"

void int24byte(uint val, uchar &b1, uchar &b2, uchar &b3, uchar &b4){
    b1=(val>>24)&0x0FF;b2=(val>>16)&0x0FF;b3=(val>>8)&0x0FF;b4=(val)&0x0FF;
}
void byte42int(uchar b1, uchar b2, uchar b3, uchar b4, uint &val ){
      val=((b1<<24)|(b2<<16)|(b3<<8)|(b4) );
}


int getScaleIndex(float lsc, vector<float> logscales, int middle_scale){
  int indl=-1,flag=1;
  float mindl=10000,dsl;
  //for(uint i=0;i<logscales.size();i++)cout << i << " allsc "<< sc << " " << logscales[i] << " " << middle_scale << endl;

  if(lsc>=0){
    for(uint i=middle_scale;i<logscales.size() && flag;i++){
      dsl=fabs(lsc-logscales[i]);
      //cout << "possc "<< lsc << " " << logscales[i] << " " << i << " dist " << dsl << endl;
      if(mindl>dsl){
	mindl=dsl;
	indl=i;      
      }     
    }  
  }else{
    for(int i=middle_scale-1;i>=0 && flag;i--){
      dsl=fabs(lsc-(logscales[i]));
      //cout << "negsc "<< sc << " " << logscales[i] << " " << dsl << endl;
      if(mindl>dsl){
	mindl=dsl;
	indl=i;      
      }
    }  
  }
  //cout << mindl << " " << indl << endl;getchar();
  return indl;
}



float interpPeak(float a, float b, float c)
{
	if (b < 0.0) {
		a = -a; b = -b; c = -c;
	}
    //assert(b >= a  &&  b >= c);
	if((b>a && b<c) || (b<a && b > c) || (b==a && b==c))return 0.0;
	return 0.5 * (a - c) / (a - 2.0 * b + c);
}


void getEigen(float a, float b, float c, float d, float &l1, float &l2){
       float trace = a+d;
       float delta1=(trace*trace-4*(a*d-b*c));
       if(delta1<0){	 
	 l1=1;l2=100;
	 return;
       } 
       float delta = sqrt(delta1);
		    
       l1 = (trace+delta)/2.0;
       l2 = (trace-delta)/2.0;
}


void invSqrRoot(float &a, float &b, float &c, float &r_l, float &r_m, float &ea){
  
  //cout << "in a " << a << " b " << b << " c " << c << endl;
  float cos_2t = a - c;
  float sin_2t =  2*b;
  float r = sqrt(cos_2t*cos_2t + sin_2t*sin_2t);
  float u,v,l,m;
  if (r == 0) {
    // We have a = c and b = 0, i.e. the matrix is diagonal.
    u = 1;
    v = 0;
    l = a;
    m = c;
    ea=0;
    return;
  }else{ 
    cos_2t /= r;
    sin_2t /= r;
    //if(cos_2t<-1)cout << cos_2t<< endl;
    if(cos_2t<-1)cos_2t=-1;
    else if(cos_2t>1)cos_2t=1;
    
    // use half-angle formulae:
    float cos_t = sqrt((1 + cos_2t)/2);
    float sin_t = sqrt(1 - cos_t*cos_t);
    if (sin_2t < 0)
      sin_t = -sin_t;
    u = cos_t;  
    v = sin_t;
    
    l = a*u*u + 2*b*v*u + c*v*v;
    m = a*v*v - 2*b*u*v + c*u*u;
    ea=atan2(u,-v);
    //ea=(ea>0)?ea:ea+M_PI; 
    //cout << " l "<< l << " m "<< m << " r " << r<< " cos_t " << cos_t<< " sin_t "<< sin_t<< " atan " << 180*atan2(u,-v)/M_PI << endl;  
  }
 
  if ((l >= 0) && (m >= 0)) {
    r_l = 1/sqrt(l);
    r_m = 1/sqrt(m);
    float x=sqrt(r_l*r_m);
    r_l/=x;
    r_m/=x;
    if(r_l/r_m<0.1){
      r_l=3;
      r_m=0.33;
    }
    a = r_l*u*u + r_m*v*v;
    b = r_l*u*v - r_m*u*v;
    c = r_l*v*v + r_m*u*u;
    //cout << "out a " << a << " b " << b << " c " << c << endl;
  }else{
    cout << "out a " << a << " b " << b << " c " << c << " r " << r << " cos_2t " << cos_2t<< " sin_2t "<< sin_2t<<  endl;
    cout << "errors "<<  " l "<< l << " m "<< m << " r " << r<< " cos_t " << u<< " sin_t "<< v<< " atan " << 180*atan2(u,-v)/M_PI <<endl;
  }  
}



 float distEuc(float* vec1,float* vec2,int size){
  float dist=0,d;    
    for(int i=0;i<size;i++){
      d=vec1[i]-vec2[i];
	dist+=d*d;
    }
    return dist;
}

float distEucSqrt(float* vec1,float* vec2,int size){
  float dist=0,d;    
    for(int i=0;i<size;i++){
      d=vec1[i]-vec2[i];
	dist+=d*d;
    }
    return sqrt(dist);
}




float distEuc(float* vec1,float* vec2,int size, float thres){
  float dist=0,d;    

  for(int i=0;i<10;i++){
      d=vec1[i]-vec2[i];
	dist+=d*d;
  }
  if(dist>thres)return dist;
  for(int i=10;i<25;i++){
    d=vec1[i]-vec2[i];
    dist+=d*d;
  }
  if(dist>thres)return dist;
  for(int i=25;i<size;i++){
    d=vec1[i]-vec2[i];
    dist+=d*d;
  }
  
  return dist;
}
 

float distMahalanobis(float *vec1, float *vec2, float **covMat, int dim){
   float dist=0;
       
   float *A = new float[dim];
   float *B = new float[dim];
   for(int a=0;a<dim;a++){
     A[a]=vec1[a]-vec2[a];
   }    
   //getchar();
   // Sin = 1/CovMat
   //Matrix *Sin=new Matrix(covMat->inverse());
   
   for(int row=0;row<dim;row++){
     B[row]=0;
     for(int col=0;col<dim;col++){
       B[row]+=covMat[row][col]*A[col];
     }
   }
   
   for(int col=1;col<dim;col++){
     dist+=A[col]*B[col];
   }
   delete []A; delete []B;//delete Sin;
   
   return (dist);    
} 


float square(float a){
  return a*a;
}

float distCSift(float* vec1,float* vec2, int size){
   float dist=0,tmp;
   float adist=0;
   for(int i=0;i<size;i+=2){
    adist=fabs(vec1[i+1]-vec2[i+1]);
    adist=(adist<M_PI)?adist:(M_2PI-fabs(vec1[i+1]-vec2[i+1]));
    adist/=M_PI;adist=1-adist;//adist=0.5;
    if(adist>1 || adist<0)cout << endl << "angle error in distCSift" << endl;
    tmp=((square(vec1[i]-adist*vec2[i])+square(adist*vec1[i]-vec2[i]))/(1+square(adist)));
    dist+=tmp;
    //cout << "o1 "<< vec1[i+1]<< " o2 " <<vec2[i+1] << "v1 "<< vec1[i]<< " v2 " <<vec2[i] << " adist " << adist<< endl;
   }
   //cout << "dist " << dist<< endl;getchar();
  return dist;
}

float distCSift1(float* vec1,float* vec2, int size){
   float dist=0, tmp=0;
   float adist=0;
   for(int i=0;i<size;i+=2){
    adist=fabs(vec1[i+1]-vec2[i+1]);
     adist=(adist<M_PI)?adist:(M_2PI-fabs(vec1[i+1]-vec2[i+1]));
    adist/=M_PI;adist=1-adist;
    if(adist>1 || adist<0)cout << endl << "angle error in distCSift" << endl;
    if(adist>0.8)tmp=square(vec1[i]-vec2[i]); 
    else tmp=(square(vec1[i])+square(vec2[i]));
    dist+=tmp;
  }
  return dist;
}


float distHistogram(float* histos1,float* histos2,int taillehisto){

  float distance;
 
  distance = 0.0;
 
  for(int i = 0 ; i < taillehisto ; ++i)
    {
      //  cerr << " histos1[i] " << histos1[i] << " histos2[i] " << histos2[i] << endl;
 
      distance += (int)fabs(histos1[i] - histos2[i]);
    }
 
  return distance;
}

float distChi(float* histos1,float* histos2,int taillehisto){

  float somme;
  float distance;
 
  distance = 0.0;
 
  for(int i = 0 ; i < taillehisto ; ++i)
    {
      //  cerr << " histos1[i] " << histos1[i] << " histos2[i] " << histos2[i] << endl;
 
      somme = (histos1[i] + histos2[i])*(histos1[i] + histos2[i]);
 
      if(somme)
        distance += ((histos1[i] - histos2[i])*(histos1[i] - histos2[i]))/somme;
    }
 
  return distance;
} 

float distL1(float* histos1,float* histos2,int taillehisto){

  float distance;
 
  distance = 0.0;
 
  for(int i = 0 ; i < taillehisto ; ++i)
    {
      //  cerr << " histos1[i] " << histos1[i] << " histos2[i] " << histos2[i] << endl;
 
      distance += fabs((histos1[i] - histos2[i]));
 
    }
 
  return distance;
}



float mean(vector<float> &vec){
  float mean=0;
  for(int i=0;i<(int)vec.size();i++)mean+=vec[i];
  return (mean/vec.size());
  
}

float max(vector<float> &vec){
  float max=vec[0];
  for(int i=1;i<(int)vec.size();i++){
    if(max<vec[i])max=vec[i];
  }
  return (max); 
}
float min(vector<float> &vec){
  float min=vec[0];
  for(int i=1;i<(int)vec.size();i++){
    if(min>vec[i])min=vec[i];
  }
  return (min); 
}

float median2(vector<float> &vec, float thres){
  
  float med=thres*vec.size();
  float fl=vec[(int)floor(med)];
  float cl=vec[(int)ceil(med)];
  float median= fl+(cl-fl)*(med-floor(med));
  return median;
}

float stdev(vector<float> &vec, float center){
  return sqrt(variance(vec,center));
}

float variance(vector<float> &vec, float center){
  float  var=0;
  float  tmp=0;
  for(int i=0;i<(int)vec.size();i++){
    tmp=(vec[i]-center);
    var+=(tmp*tmp);
  }  
  
  return (var/(vec.size()-1));
}

void sort(vector<float> &vec){
  float min;
  int imin;
  vector<float> vec_out;
  do{
    min=vec[0];imin=0;
    for(int i=1;i<(int)vec.size();i++){
      if(vec[i]<min){min=vec[i];imin=i;}
    } 
    vec_out.push_back(min);
    vec.erase(vec.begin()+imin);
  }while((int)vec.size()>0);
  vec=vec_out;
}

float median(vector<float> &vec){
  float min;
  int imin;
  vector<float> vec_out;
  int siz = (int)vec.size(); 
  if(siz<2)return 0;
  do{
    min=vec[0];imin=0;
    for(int i=1;i<(int)vec.size();i++){
      if(vec[i]<min){min=vec[i];imin=i;}
    }
    vec_out.push_back(min);
    vec.erase(vec.begin()+imin);
  }while((int)vec.size()>(siz>>1));
  float out = (vec_out[vec_out.size()-1]+vec_out[vec_out.size()-2])/2.0;
  vec=vec_out; 
  return out; 
}


float gamma(float n){    
    float tmp=0;
    n=rint(2*n)/2;
    if(rint(n)==n){
        tmp=1;
	for(int i=1;i<(int)n;i++)tmp*=i;
    }else if(rint(n)!=n){
	int ni=(int)floor(n);
        tmp=0.5;
	for(int i=1;i<ni;i++)tmp*=(i+0.5);
	tmp*=sqrt(M_PI);
    }
    return tmp;
}

float chi2(float x, float n){
  float a = 1/(sqrt(pow((double)2,(double)n))*gamma(0.5*n));    
    return (a*pow((double)x,(double)(0.5*(n-2)))*exp(-0.5*x));
}

float round(float x, int n){
    float mult=pow((double)10,(double)n);
    return rint(mult*x)/mult;
}

float probChi2(float thres, float n){
    float sum=0;
    float i=0;
    while(i<thres){
	sum+=(chi2(i,n));
	i+=0.1;
    }
    return sum/10.0;
}


float thresChi2(float prob, float n){
    float sum=0;
    float thres=0;
    prob*=10;
    while(sum<prob){
	sum+=(chi2(thres,n));
	thres+=0.1;
    }
    return thres;
}

//???
vector<float> hist(vector<float> &vec,float min, float max, int bins){
  float db = (max-min)/bins;
  int i=0;while(vec[i]<min)i++;
  float thres=min+db;
  float bin=0;
  vector<float> hist;
  while( vec[i]<max){
    while(vec[i]<thres && i<(int)vec.size() && vec[i]<max){bin++;i++;}
    hist.push_back(((bin/(int)vec.size())*100));
    bin=0;
    thres+=db;
  }
  return hist;
}


#define ELEM_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }
float median(float *arr, int n){
  int low, high ;
  int median;
  int middle, ll, hh;
  low = 0 ; high = n-1 ; median = (low + high) / 2;
  for (;;) {
    if (high <= low) /* One element only */
      return arr[median] ;
    if (high == low + 1) { /* Two elements only */
      if (arr[low] > arr[high])
        ELEM_SWAP(arr[low], arr[high]) ;
      return arr[median] ;
    }
    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high]) ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high]) ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low]) ELEM_SWAP(arr[middle], arr[low]) ;
    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;
    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
      do ll++; while (arr[low] > arr[ll]) ;
      do hh--; while (arr[hh] > arr[low]) ;
      if (hh < ll)
        break;
      ELEM_SWAP(arr[ll], arr[hh]) ;
    }
    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;
    /* Re-set active partition */
    if (hh <= median)
      low = ll;
    if (hh >= median)
      high = hh - 1;
  }
}
#undef ELEM_SWAP


float computeOverlap(float x1s, float y1s,float x1e, float y1e,float x2s, float y2s,float x2e, float y2e){
  int margin=0; 
  //is in
  float xs=(x1s>x2s)?x1s:x2s;
  xs=xs-margin;
  float ys=(y1s>y2s)?y1s:y2s;
  ys=ys+margin;
  float xe=(x1e<x2e)?x1e:x2e;
  xe=xe+margin;
  float ye=(y1e<y2e)?y1e:y2e;
  ye=ye+margin;
  
 
  float dx=xe-xs;
  float dy=ye-ys;
  float aintersection=0;
  float aunion=1;
  if(dx>0 && dy>0){
    aintersection=dx*dy;
    aunion = (x1e-x1s)*(y1e-y1s)+(x2e-x2s)*(y2e-y2s)-aintersection;
   // aintersection=aintersection/aunion;
    
/*  cout <<"point 1 "<<  x1s << " " << y1s << " " << x1e << " " << y1e << endl;
  cout <<"point 2 "<< x2s << " " << y2s << " " << x2e << " " << y2e << endl;
  cout << xs << " " << ys << " " << xe << " " << ye <<  " " << dx << " " << dy << " " << aintersection << " " << aunion << "  "<< aintersection/aunion<< endl;*/
  //getchar();

    aintersection=aintersection/aunion;

  } 
  //  cout <<mdl1->xs << " " << mdl1->xe << " " <<  mdl1->ys<< " " <<  mdl1->ye  << endl;
  // cout <<mdl2->xs << " " << mdl2->xe << " " <<  mdl2->ys<< " " <<  mdl2->ye  << endl;
  // cout <<xs << " e " << xe << "  s " <<  ys<< " e " <<  ye  << endl;
  // cout << dx << " x " << dy << " " <<  intersection <<endl;getchar();
  return sqrt(aintersection); 

}

float computeOverlapNorm(float x1, float y1,float s1, float x2, float y2,float s2, float sf){
  
  float scaling=sf/s1;
  
  s1=s1*scaling;
  s2=s2*scaling;
  
  float x1s=x1-s1;
  float y1s=y1-s1;
  float x1e=x1+s1;
  float y1e=y1+s1;
  float x2s=x2-s2;
  float y2s=y2-s2;
  float x2e=x2+s2;
  float y2e=y2+s2;
  

  return computeOverlap( x1s,  y1s, x1e,  y1e, x2s,  y2s, x2e,  y2e);

}


