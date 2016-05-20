#ifndef _feature_h_
#define _feature_h_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;
#include <vector>
#include "../ImageContent/imageContent.h"
#include "../gauss_iir/gauss_iir.h"
#include "../util/util.h" 
#include "../segment/segment.hh" 
#include "../mser/extrema/libExtrema.h"
#include "../matrix/matrix.h"
using namespace  extrema;

#define SCALE_FACTOR  1

#define FALSE 0
#define TRUE 1


#define  AFFINE 16 

#define  ADETECTOR           0x000FFFFF
#define  DETECTOR            0x000FFFFF
//#define   DETECTOR         0x00000CCCC
//#define  DETECTOR          0x00000AAAA

//two types depending on the extremum sign
#define  DHARHES             0x00000002
#define  DHESSIAN            0x00000008
#define  DHARRIS             0x00000020
#define  DSEDGE              0x00000080
#define  DMSER               0x00000200
#define  DUNIFORM            0x00000800
#define  DHAR                0x00002000
#define  DSEGM               0x00002000
#define  DHES                0x00008000
#define  FFACE               0x00020000
#define  OSIFT               0x00080000



#define NO_MOTION            -100

#define  DESCRIPTOR          0xFFF00000
#define  DSIFT               0x00100000
#define  DAUDIO              0x00200000
#define  DCC                 0x00400000
#define  DCOLOR              0x00800000
#define  DMOM                0x00800000
#define  DPCA                0x01000000
#define  DKOEN               0x01000000
#define  DSHAPE              0x01000000
#define  DMOTION             0x01000000
#define  DCSIFT              0x01000000
#define  DSPIN               0x01000000
#define  DGLOH               0x01000000
#define  DJLA                0x01000000
#define  DCLBP               0x01000000
//Feature displaying flags
#define  DN                3
#define  DP                0
#define  DC                1
#define  DE                2


#define OAFFINE 16

extern DARY *patch_mask;
extern int PATCH_SIZE;
//#define PATCH_SIZE 21
//#define PATCH_SIZE 33 

void initPatchMask(int size);
class FeatureDescriptor;
typedef FeatureDescriptor FD;


#include "../sift/sift.h"

#define OX 0
#define OY 1
#define OSCALE 2
#define OFSCALE 3
#define OANGLE 4
#define OFANGLE 5
#define OMOTION 6
#define ODIST 7
#define OWEIGHT 8
#define OAREA 9
#define OEL 10
#define OEA 11
#define OOBJ 12
#define OIMG 13
#define OPARAMS_NB 14

class BOcc{
  public:
    BOcc(){}    
    uchar x;//round(255*x/width) // polar radius
    uchar y;//round(255*y/width) // polar angle
    uchar s;//round(log(scale)*40)<255 scale
    uint o;// image/object
    uchar w;//weight    
    uchar a;//angle, round(255*a/PI)   
    uint df1;
    uint df2;
};

typedef vector<BOcc *> BOccList;
typedef vector<BOccList> WordbookBOcc;

class Occ{
  public:
    Occ(){}    
    int x;
    int y;
    int s;
    int obj;
    float weight;    
    float angle;    
    
};


class Occurence{
 public:
  float *par;  
  FD *feat;

  void allocPar(){
    par = new float[OPARAMS_NB];
    bzero(par,OPARAMS_NB*sizeof(float));
  }
  
  Occurence(int x_in, int y_in, int scale_in, float angle_in, float weight_in, int obj_in){
    allocPar();
    par[OX]=x_in;
    par[OY]=y_in;
    par[OSCALE]=scale_in;
    par[OFSCALE]=1;
    par[OANGLE]=angle_in;
    par[OFANGLE]=0;
    par[OMOTION]=NO_MOTION;
    par[ODIST]=0;
    par[OWEIGHT]=weight_in;
    par[OOBJ]=obj_in;
  }

  Occurence(){allocPar();}
  
  Occurence(Occurence *occ){
    allocPar();
    memcpy(par, occ->par, sizeof(float)*OPARAMS_NB);
  }
  
  Occurence(float x_in, float y_in, float scale_in,float angle_in,float fangle_in, float area_in,float motion_in, float dist_in, float weight_in, int obj_in, int img_in){
    allocPar();
    par[OX]=x_in;
    par[OY]=y_in;
    par[OSCALE]=scale_in;
    par[OFSCALE]=1;
    par[OANGLE]=angle_in;
    par[OFANGLE]=fangle_in;
    par[OMOTION]=motion_in;
    par[ODIST]=dist_in;
    par[OWEIGHT]=weight_in;
    par[OOBJ]=obj_in;
    par[OIMG]=img_in;
    par[OAREA]=area_in;
  
    feat=NULL;
  }
  ~Occurence(){
    delete []par;
  }

  void Cout(){
    cout << par[OX] << " "<< par[OY] << " s " << par[OFSCALE] << " " << par[OWEIGHT] << endl;
  }

  void out2(ofstream &out, int width_2, int height_2, float image_scale){
    //cout << x << " x "<< y << " "  << width_2 << " "<< fscale << " " << image_scale<< endl;
    out  <<  par[OOBJ] << " " << (par[OX]-(width_2*OFSCALE))/image_scale << " "<< (par[OY]-height_2*par[OFSCALE])/image_scale << " " << (par[OX]+(width_2*par[OFSCALE]))/image_scale << " "<< (par[OY]+height_2*par[OFSCALE])/image_scale << " " << par[OWEIGHT] << endl;
  }
};

#define SegOriBins 12
class Segment;
class Segment{
 public:
  float scale, ang_mean;
  int x_mean,y_mean;
  int valid;
  vector<int> x;
  vector<int> y;
  vector<float> angle;
  vector<float> grad;
  vector<Segment *> segs;
  float *hist;
  Segment(){
    ang_mean=0;
    valid=0;
    hist = new float[SegOriBins];
  }
  Segment(Segment * seg){
    copy(seg);
  }
  void copy(Segment * seg);

  ~Segment(){
    x.clear();
    y.clear();
    angle.clear();
    for(uint i=0;i<segs.size();i++)delete segs[i];
    segs.clear();
    delete []hist;
  }
};



#define PSIZE 0
#define DSIZE 1
#define XPOS 2
#define YPOS 3
#define FEATURNESS 4
#define C_SCALE 5
#define ANGLE 6
#define OBJ 7
#define TYPE 8
#define LAP 9
#define EXTR 10
#define MI11 11
#define MI12 12
#define MI21 13
#define MI22 14
#define IMG 15
#define D_SCALE 16
#define WEIGHT 17
#define VAR 18
#define SIM 19
#define RADIUS 20
#define AREA 21
#define NEG_WEIGHT 22
#define MEAN_DIST 23
#define NBF 24
#define L1 25
#define L2 26
#define TREE_LEV 27
#define INT_SIG 28
#define DER_SIG 29
#define EANGLE 30
#define INT_LEV 31
#define DER_LEV 32
#define ID 33
#define IM_ID 34
#define MOTION 35
#define A_MOTION 36
#define SYMM 37
#define PARAMS 38

//feature save formats
#define FBIN 5
#define COMPACT_FBIN 6
#define OXFORD_FTEXT 1
#define FTEXT 2
#define AUDIO_FTEXT 4



class  FeatureDescriptor{

 protected:
 public:
    float *vec;  
    float *par;  
    char *imagename;
    DARY *cluster_dist;
    vector<FeatureDescriptor *> features;
    vector<Occurence *> occurences;
    vector<uint> occ;
    vector<float> neg_match;
    vector<float> pos_weight;//per object

 public:
     FeatureDescriptor(){init();} 
     FeatureDescriptor(int psize){psize=(psize<PARAMS)?PARAMS:psize;init(psize);} //some internal parameters at positions beyond psize wouldn't be copied if <PARAMS
     void init();
     void init(int psize);
     
     FeatureDescriptor(float xin, float yin, float scale_in, float featureness_in);

     void copy(FeatureDescriptor* ds);

    /*****READ WRITE***/
    void readMSER( ifstream &input);
    void readCommon( ifstream &input, int size_in);
    void writeCommon( ofstream &output);
    void write( ofstream &output);
    void read( ifstream &input, int psize_in, int dsize_in);	
    void writeBin(FILE *fid, float *buf);
    void readBin(FILE *fid, int size, float *buf);

    void allocVec(int);
    inline float * getVec(void){return vec;}
    inline float * getPar(void){return par;}

    inline float getV(int i){ if(i<par[DSIZE])return (vec[i]);else return 0;}    
    inline void setV(int i, float val){if(i<par[DSIZE] && i>=0)vec[i]=val;}

    inline int const getSize() const {return (int)par[DSIZE];} 
    inline void setSize(int size_in){par[DSIZE]=size_in;} 

    inline int const getID() const {return (uint)par[ID];} 
    inline void setID(uint id_in){par[ID]=id_in;} 
  
    inline uint const getImID() const {return (uint)par[IM_ID];} 
    inline void setImID(uint id_in){par[IM_ID]=id_in;} 
  
    inline float const getSim() const {return par[SIM];} 
    inline void setSim(float sim_in){par[SIM]=sim_in;}   
 
    inline float const getRadius() const {return par[RADIUS];} 
    inline void setRadius(float radius_in){par[RADIUS]=radius_in;}    

    inline int const getTreeLevel() const {return (int)par[TREE_LEV];} 
    inline void setTreeLevel(int tree_lev_in){par[TREE_LEV]=tree_lev_in;}    

    inline float const getWeight() const {return par[WEIGHT];} 
    inline void setWeight(float weight_in){par[WEIGHT]=weight_in;}    

    inline float const getArea() const {return par[AREA];} 
    inline void setArea(float area_in){par[AREA]=area_in;}    

    inline float const getMeanDist() const {return par[MEAN_DIST];} 
    inline void setMeanDist(float d_in){par[MEAN_DIST]=d_in;}    

    inline uint const getNbf() const {return (uint)par[NBF];} 
    inline void setNbf(uint nbf_in){par[NBF]=nbf_in;}    

    inline unsigned int   const  getType() const { return ( unsigned int)par[TYPE];}
    inline void     setType(unsigned int type_in)  {par[TYPE]=type_in;}

    inline uint   const  getSymm() const { return (uint)par[SYMM];}
    inline void     setSymm(uint symm_in) {par[SYMM]=symm_in;}

    inline int   const  getObj() const { return (int)par[OBJ];}
    inline void     setObj(int obj_in)  {par[OBJ]=obj_in;}
    
    inline int   const  getImg() const { return (int)par[IMG];}
    inline void     setImg(int img_in)  {par[IMG]=img_in;}
    
  
    /****CORNERNESS***/
    inline float const getFeatureness(void) const { return par[FEATURNESS];}
    inline void setFeatureness(float featureness_in) {par[FEATURNESS]=featureness_in;}
    
    inline void setMi(float m11,float m12,float m21,float m22) 
    {par[MI11]=m11;par[MI12]=m12;par[MI21]=m21;par[MI22]=m22;}
    inline void setMi(float m11,float m12,float m21,float m22,float e1,float e2, float ea) 
      {par[MI11]=m11;par[MI12]=m12;par[MI21]=m21;par[MI22]=m22;par[L1]=e1;par[L2]=e2;par[EANGLE]=ea;}
    inline float   const  getMi11() const { return par[MI11];}
    inline float   const  getMi12() const { return par[MI12];}
    inline float   const  getMi21() const { return par[MI21];}
    inline float   const  getMi22() const { return par[MI22];}
    inline float   const  getL1() const { return par[L1];}
    inline float   const  getL2() const { return par[L2];}
    
    /****LOCALISATION***/
    inline float   const  getX() const { return par[XPOS];}
    inline float   const  getY() const { return par[YPOS];}
    inline void    setX_Y(float x_in,float y_in)  {par[XPOS]=x_in;par[YPOS]=y_in;}


    /****LAPLACIAN VALUE***/
    inline float const getLap() const { return par[LAP];}
    inline void setLap(float lap_in) {par[LAP]=lap_in;}
    inline void setExtr(int ex) {par[EXTR]=ex;}
    inline void setMax() {par[EXTR]=1;}
    inline void setMin() {par[EXTR]=-1;}
    inline int   const  getExtr() const { return (int)par[EXTR];}
    

    /****ANGLE*****/
    inline float   const  getAngle() const { return par[ANGLE];}
    inline void    setAngle(const float angle_in){ par[ANGLE]=angle_in;}

    /****MOTION*****/
    inline float   const  getMotion() const { return par[MOTION];}
    inline void    setMotion(const float motion_in){ par[MOTION]=motion_in;}

    inline float   const  getEAngle() const { return par[EANGLE];}
    inline void    setEAngle(const float eangle_in){ par[ANGLE]=eangle_in;}

    inline float   const  getVar() const { return par[VAR];}
    inline void    setVar(const float var_in){ par[VAR]=var_in;}

    /*****SCALE******/
    inline float const  getScale() const { return par[C_SCALE];}
    inline void setScale(const float scale_in){par[C_SCALE]=scale_in;}
    inline void setIntSig(const float sig_in){par[INT_SIG]=sig_in;}
    inline float const getIntSig() const { return par[INT_SIG];}
    inline void setDerSig(const float sig_in){par[DER_SIG]=sig_in;}
    inline float const  getDerSig() const { return par[DER_SIG];}
    inline int const  getIntLev() const { return (int)par[INT_LEV];}
    inline int const  getDerLev() const { return (int)par[DER_LEV];}
    inline void setIntLev(int lev_in) {par[INT_LEV]=lev_in;}
    inline void setDerLev(int lev_in) {par[DER_LEV]=lev_in;}
    
    
    inline const char*   getImageName() const { return imagename;}
    void setImageName(const char* name);
    
    //    void addOccurence(float x, float y, float scale,float angle,float fangle, float area, float dist, float sim, int obj, uint img){
    //      occurences.push_back(new Occurence(x,y,scale,angle,fangle,area,dist,sim,obj,img));
    //    }  
    void findNearestNeighbor(vector<FeatureDescriptor*> features, int dim, int &ind, float &sim);

    void changeBase(float* mat);
    void pca(int dim,  float *avg, float *base);
    void normalizeVect(); 
    void Cout(int dim=3);


    ~FeatureDescriptor();
   
};


typedef vector <FeatureDescriptor*> FeatureDescList;
typedef vector <FeatureDescList> FeatureDescSequence;



inline float square(float a){return a*a;}

char * getName( unsigned long  type);

void normalize(DARY * img,int x, int y, float radius);

void deleteDescriptors(vector<FeatureDescriptor *> &desc);

DARY* normalizeAffine(DARY *image, float x, float y, float c_scale, 
		      float angle, float mi11, float mi12,float mi21, float mi22, float &scal,
		      DARY *patch, float DESC_CALE);
void computeAffineDescriptor( DARY *imgbs, DARY *patch, float scal, int DESC,
			FeatureDescriptor *cor, vector<FeatureDescriptor *> &desc);
void computeAffineDescriptor( DARY *image, DARY *patch, int DESC, 
			FeatureDescriptor *cor, vector<FeatureDescriptor *> &desc, float DESC_CALE);
void computeAffineDescriptors(DARY *image,  vector<FeatureDescriptor *> &desc, int DESC, float DESC_CALE);

void computeDescriptor( DARY *gpatch, DARY *opatch, uint DESC,
			FeatureDescriptor *cor, vector<FeatureDescriptor *> &desc, int noangle);
void computeDescriptors(DARY *img, vector<FeatureDescriptor*>&features, Params *params);


void loadBinFeatures(const char* points_out, vector<FeatureDescriptor*> &cor);
void saveBinFeatures(vector<FeatureDescriptor*> cor, const char* points_out, int format=FBIN);

void saveFeatures( vector<FeatureDescriptor*> cor1, const char* points_out, int format=FBIN);
void loadFeatures( const char* points1, vector<FeatureDescriptor*> &cor1, int format=FBIN);
void loadFeatures(vector<char *> filenames, vector<FD *> &features, ulong max_number_features, ulong max_memory_features, int format=FBIN);
void loadBOcc(const char* points_out, vector<BOccList> &wbocc);
void saveBOcc(vector<BOccList> wbocc, const char* points_out);
void loadAndProjectFeatures(vector<char *> filenames, vector<FD *> &features, ulong max_number_features, ulong max_memory_features, int format, float *spca_mean, float *spca_base, uint desc_dim);
void loadAndProjectFeatures( const char* points1, vector<FeatureDescriptor*> &cor1, int format, float *spca_mean, float *spca_base, uint desc_dim);





void computeSiftOrig(DARY *dx, DARY *dy, DARY *grad, DARY *gori , FeatureDescriptor *ds, int oriSize, int locSize, float max_gori);
void computeSiftOrig(DARY *dx, DARY *dy ,DARY *grad, DARY *gori , FeatureDescriptor *ds);
void computeSift(DARY *dx, DARY *dy , DARY *grad, DARY *gori ,  FeatureDescriptor *ds);
void computeESift(DARY *dx, DARY *dy ,DARY *grad, DARY *gori , FeatureDescriptor *ds);
void computeShape(DARY *dx, DARY *dy ,DARY *grad, DARY *gori, FeatureDescriptor *ds);


//void computeJLA(DARY *patch, FeatureDescriptor *ds);
void symmetricSift(FD *fd, float x);
void computeSift(Segment *seg, FeatureDescriptor *ds);
void computeSift(DARY *patch, FeatureDescriptor *ds);
void computeESift(DARY *patch, FeatureDescriptor *ds);
//void computeMoments(DARY *patch, FeatureDescriptor *ds);
//void computeKoen(DARY *patch, FeatureDescriptor *ds);
//void computeCF(DARY *patch, FeatureDescriptor *ds);
void computeShape(DARY *patch, FeatureDescriptor *ds);
//void computeSpin(DARY *patch, FeatureDescriptor *ds);
//void computePcaSift(DARY *patch, FeatureDescriptor *ds);
//void computeCC(DARY *patch, FeatureDescriptor *ds);
void computeColorDesc(DARY *cimg, FD *f);
void computeColorSift(DARY *cimg, FD *f);

void extractColorDescriptor(const char *name, vector<FD *> &features,   Params *params);
void extractColorLBP(DARY *img, DARY *gray, vector<FD *> &features,  Params *params);
void extractSegmentSift(DARY *image, vector<FD *> &features,  Params *params);

void normalizeFeature(float *vec, int size);
void computeHistAngle(DARY *grad, DARY *ori, int x, int y, int size,float &angle);
void computeHistAngle(int x, int y, int radius, DARY *grad, DARY *ori,vector<float> &angles);


void displayFeatures(DARY *image, vector<FeatureDescriptor*> cor, char *filename, float color, int circle=0);
void displayFeatures(const char *filein, vector<FeatureDescriptor*> cor, char *filename, float color, int circle=0);

void separateFeatures(vector<FD *> clusters, vector<FD *> &feat);

void cannyEdges(DARY *img, DARY *edge,  float scale, float lower_threshold, float higher_threshold);
void cannyEdges(DARY *dx, DARY *dy, DARY *edge,  float lower_threshold, float higher_threshold);
void cannyEdges(DARY *img, DARY *edge,  float lower_threshold, float higher_threshold);
void cannyEdges(DARY *dx, DARY *dy, DARY *grad, DARY *edge,  float lower_threshold, float higher_threshold);
void cannyEdges(DARY *dx, DARY *dy, DARY *grad, DARY *edge, DARY *tmp_edge, float lower_threshold,
                float higher_threshold, int &edge_nb);

void cannySegments(DARY *dx, DARY *dy, DARY *grad,DARY *ori, DARY *cedge, vector<FD*> &segs, float lower_threshold,float higher_threshold);


void pca(vector<FD*> &features,float *mvec, float *base, int newdim);
void pcaBase(vector<FD*> &features,float *eigen,  float *mvec, float *base,  int newdim);
void meanVarBase(vector<FD*> &features, char * filename);
void savePCAVectors(const char *filename,  float *eigen, float *vmean, uint size1 ,float *vbase, uint size2, uint desc_dim);
int loadPCAVectors(const char *filename,  float *&eigen, float *&vmean, uint &mean_dim, float *&vbase, uint &base_dim, uint &desc_dim);
void meanBase(vector<FD*> features, float *mvec, int dim);
void varianceBase(vector<FD*> features, float *vvec, int dim);


void fastDense(DARY *img, vector<FeatureDescriptor*>&features, Params *params);

void fastHarris(DARY *img, vector<FeatureDescriptor*>&features, float threshold, float step, uint DESC, float reg_thes, int aff, int noangle, int dradius, int desc_nb);
void fastHarris(DARY *img, vector<FeatureDescriptor*>&features, Params *params);
void harris(DARY *img, vector<FeatureDescriptor*>&features, Params *params);

void fastHessian(DARY *img, vector<FeatureDescriptor*>&features, float threshold, float step, uint DESC, float reg_thes, int aff, int noangle, int dradius, int desc_nb, int hesdet);
void fastHessian(DARY *img, vector<FeatureDescriptor*>&features, Params *params);
void hessian(DARY *img, vector<FeatureDescriptor*>&features, Params *params);

void fastUniform(DARY *img, vector<FeatureDescriptor*>&features, Params *params);

void fastHarrisHessian(DARY *img, vector<FeatureDescriptor*>&corners, float harthres, float hesthres, float step, uint DESC, float reg_thes, int aff, int noangle, int dradius, int desc_nb);
void fastHarrisHessian(DARY *img, vector<FeatureDescriptor*>&corners, Params *params);


void fastSEdge(DARY* img, vector<FeatureDescriptor*>&features,float harthres, float hesthres,float edgeHthres,float edgeLthres,float step,uint DESC,int aff, int radius, float eigenratio, int noangle, int dradius);
void fastSEdge(DARY* img, vector<FeatureDescriptor*>&features, Params *params);

void fastEdge(DARY* img, vector<FeatureDescriptor*>&features, float edgeHthres, float edgeLthres, float step, uint DESC, int aff, int radius, int noangle, int dradius);
void fastEdge(DARY* img, vector<FeatureDescriptor*>&features, Params *params);

void fastSegm(DARY* img, vector<FeatureDescriptor*>&features, Params *params);
void fastSegm(DARY* img, vector<FeatureDescriptor*>&features, float edgeHthres, float edgeLthres, 
	      float step, uint DESC, int noangle);
void fastEGSF(DARY *image, vector<FD *> &features, Params *params);

void extractFeatures(DARY *image, vector<FD*> &features, Params *params);
void detectFeatures(const char *name, DARY *img, DARY *mask, vector<FD *> &features, int &width, int &height,  Params *params);
void detectFeatures(const char *filein, vector<FD *> &features,  int &width, int &height,  Params *params);
void regionFilter(vector<FeatureDescriptor*>&features, float dist_thres);

void KeySample(FD *fd, vector<edge> points, float **grad, float **ori, int oriSize, int locSize, float max_angle);

void medianFilter(DARY *img, int size);
void openingFilter(DARY *img, int size);
void gammaFilter(DARY *img);
void smoothHistogram(float *hist, int bins);
void findDominantAngles(float *hist, int OriBins, vector<float> &angles, float range_angle, float max_angle);
#endif
