#ifndef _util_h_
#define _util_h_
#include <sys/types.h>
#include <sys/times.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

#ifndef M_2PI
#define M_2PI 6.2831853072
#endif
#ifndef MAX_FLOAT
#define MAX_FLOAT 40000000000.0
#endif

//typedef unsigned long uint;
typedef unsigned char uchar;
typedef unsigned long ulong;

void int24byte(uint val, uchar &b1, uchar &b2, uchar &b3, uchar &b4);
void byte42int(uchar b1, uchar b2, uchar b3, uchar b4, uint &val );
int getScaleIndex(float sc, vector<float> logscale, int middle_scale);
void loadFileNames(const char *filein, vector<char *> &filenames);
void saveFileNames(vector<char *> filenames, const char *filein );

float interpPeak(float a, float b, float c);
void getEigen(float a, float b, float c, float d, float &l1, float &l2);
void invSqrRoot(float &a, float &b, float &c, float &r_l, float &r_m, float &ea);
float distEuc(float* vec1,float* vec2,int size);
float distEuc(float* vec1,float* vec2,int size, float thres);
float distEucSqrt(float* vec1,float* vec2,int size);
float distL1(float* histos1,float* histos2,int taillehisto);

inline float distMahalanobis(float *vec1, float *vec2, float **covMat, int dim);
inline float distHistogram(float* histos1,float* histos2,int taillehisto);

float distChi(float* histos1,float* histos2,int taillehisto);
inline float distCSift(float* vec1,float* vec2,int size);

inline float square(float a);
inline void sort(vector<float> &vec);
inline float variance(vector<float> &vec, float center);
inline float stdev(vector<float> &vec, float center);
 float mean(vector<float> &vec);
 float median(vector<float> &vec, float thres);
 float median(vector<float> &vec);
 float max(vector<float> &vec);
 float min(vector<float> &vec);
inline float gamma(float n);
inline float chi2(float x, float n);
inline float thresChi2(float prob, float n);
inline float probChi2(float thres, float n);
inline float round(float x, int n);
vector<float> hist(vector<float> &vec,float min, float max, int bins);
float median(float *arr, int n);
float computeOverlap(float x1s, float y1s,float x1e, float y1e,float x2s, float y2s,float x2e, float y2e);
float computeOverlapNorm(float x1, float y1,float s1, float x2, float y2,float s2, float sf);

class Timer{
  long tmu,tms;
  struct tms colin_time;
 public:
  
  Timer();
  void  stop ();
  float getTime ();
  void start ();
 
};

class Params{
 public:
  vector<char *> names;
  vector<double> fvalues;
  vector<char *> cvalues;
  void save(const char *filename);
  void load(const char *filename);
  void put( const char *name, double value);
  void put( const char *name, const char *value);
  double getValue(const char *name);
  char *getString(const char *name);
  ~Params(){
   names.clear();
   fvalues.clear();
   cvalues.clear();
  }
};

/*
class Hypothesis{
 public:
  int obj, x1,  y1,  x2,  y2,  x3,  y3,  x4,  y4,; 
  float weight;
  Hypothesis(){}
  Hypothesis(int x1_, int y1_, int x2_, int y2_, int x3_, int y3_, int x4_, int y4_, float weight_in, int obj_in){
    x1=x1_;
    y1=y1_;
    x2=x2_;
    y2=y2_;
    x3=x3_;
    y3=y3_;
    x4=x4_;
    y4=y4_;
    obj=obj_in;
    weight=weight_in;
  }

  ~Hypothesis(){}

  void Cout(){
    cout << x1 << " x " <<   y1 << " " <<   x2 << " x " <<   y2 << " " <<  x3 << " x " <<   y3 << " " <<   x4 << " x " <<   y4 << endl;
    //cout << "x "<< x<< " y "<< y << " width "<< width << " height " << height <<  " angle " << angle << " weight " << weight <<endl;
  }
}; 
 
*/



/* DEFINITIONS */
#define MAX_SIG_SIZE   100
#define MAX_ITERATIONS 500
//#define INFINITY       1e20
#define EPSILON        1e-6

/*****************************************************************************/
/* feature_t SHOULD BE MODIFIED BY THE USER TO REFLECT THE FEATURE TYPE      */
typedef unsigned char feature_t;
/*****************************************************************************/


typedef struct
{
  int n;                /* Number of features in the signature */
  feature_t *Features;  /* Pointer to the features vector */
  float *Weights;       /* Pointer to the weights of the features */
} signature_t;


typedef struct
{
  int from;             /* Feature number in signature 1 */
  int to;               /* Feature number in signature 2 */
  float amount;         /* Amount of flow from "from" to "to" */
} flow_t;



float emd(signature_t *Signature1, signature_t *Signature2,
	  float (*func)(feature_t *, feature_t *),
	  flow_t *Flow, int *FlowSize);



#endif
