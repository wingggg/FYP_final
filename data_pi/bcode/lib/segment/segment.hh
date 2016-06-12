#ifndef _segment_h_
#define _segment_h_
  
#include <iostream>
#include <iomanip>

#include <math.h>
#include <stdlib.h>
//#include "misc.h"
using namespace std;
#include <vector>
#
#define THRESHOLD(size, c) (c/size)
typedef unsigned char uchar;

typedef struct { uchar r, g, b; } rgb;


// template <class T>
//     inline T abs(const T &x) { return (x > 0 ? x : -x); };
// 
// template <class T>
//     inline int sign(const T &x) { return (x >= 0 ? 1 : -1); };
// 

// template <class T> inline T bound(const T &x, const T &min, const T &max) {
//   return (x < min ? min : (x > max ? max : x));
// }
// 
// template <class T> inline bool check_bound(const T &x, const T&min, const T &max) {
//       return ((x < min) || (x > max));
// }
// 
//     inline int vlib_round(float x) { return (int)(x + 0.5F); }
// 
//     inline int vlib_round(double x) { return (int)(x + 0.5); }
// 
//     inline double gaussian(double val, double sigma) {
//       return exp(-square(val/sigma)/2)/(sqrt(2*M_PI)*sigma);
//     }

       
typedef struct {
  float w;
  uint a, b;
  uint laba, labb;
} edge;

typedef struct{
    vector<edge> points;//segment points
    vector<edge> bpoints;//indexes of boundary points
    vector<int> blabels;//neighboring labels
    vector<float> bweights;//accumulated boundary weights
    vector<float> bsizes;//boundary length
    uint size;//nb of segment points
    int label;
    float w;
}segment;

typedef struct {
  int rank;
  int p;
  int size;
  float w;
} uni_elt;

class universe {
  public:
    universe(int elements);
    ~universe();
    int find(int x);  
    void join(int x, int y, float &w);
    int size(int x) const { return elts[x].size; }
    int num_sets() const { return num; }

  private:
    uni_elt *elts;
    int num;
};

 

void graph_segmentation(DARY *image, int k, float grad_thres, float bound_thres,  int min_size, int &num_ccs, DARY *output);
void segment_image(DARY *smooth_r, DARY *smooth_g, DARY *smooth_b, float k, float grad_thres, float bound_thres, int min_size, int &num_ccs, DARY *igrad, DARY *iori, DARY *labout, vector<segment *> &epoints);
void graph_segmentation(DARY *image, DARY *output, Params *params);







#endif
