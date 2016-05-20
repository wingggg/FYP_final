#ifndef _sift_h_
#define _sift_h_

/* From the standard C libaray: */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "../descriptor/feature.h"
#include "../ImageContent/imageContent.h"


/*--------------------------- Useful macros -----------------------------*/

/* Slight overestimate of PI, so orientations can be in range [0,PI]. */
#define PI 3.1415927

#define MAX(x,y)  (((x) > (y)) ? (x) : (y))
#define MIN(x,y)  (((x) < (y)) ? (x) : (y))
#define ABS(x)    (((x) < 0) ? (-(x)) : (x))

/* Following defines TRUE and FALSE if not previously defined. */
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/* Given the name of a structure, NEW allocates space for it and returns
     a pointer to the structure.
*/
#define NEW(s,pool) ((struct s *) MallocPool(sizeof(struct s),pool))

/* Following defines the maximum number of different storage pools. */
#define POOLNUM 10

/* Each of the following assigns a unique number less than POOLNUM to
   a particular pool of storage needed for this application. 
*/
#define PERM_POOL  0     /* Permanent storage that is never released. */
#define MODEL_POOL 1     /* Data for current set of object models. */
#define IMAGE_POOL 2     /* Data used only for the current image. */
#define TREE_POOL  3     /* Nodes for current k-d tree.. */


/*---------------------------- Structures -------------------------------*/

/* Data structure for a float image.
*/
typedef struct ImageSt {
    int rows, cols;          /* Dimensions of image. */
    float **pixels;          /* 2D array of image pixels. */
} *Image;


/* This structure describes a keypoint that has been found in an image.
*/
typedef struct KeypointSt {
    float row, col;      /* Row, column location relative to input image.  */
    float scale;         /* Scale (relative to input image). */
    float ori;           /* Orientation in radians (-PI to PI). */
    float extr;           /* Orientation in radians (-PI to PI). */
    struct KeypointSt *next;   /* Links keypoints into a global list. */
} *Keypoint;


/*------------------------------- Externals -----------------------------*/

extern const int Tracing, FirstOnly, SampleSize, VecLength;
extern int OutputAllKeys, DoubleImSize;
//extern const float InitSigma;


/*-------------------------- Function prototypes -------------------------*/
/* The are prototypes for the external functions that are shared
   between files.
*/

Keypoint GetKeypoints(Image image,char *name);
void DisplayKeypoint(Keypoint k, int smalltick);

/* Display routines from main.c */
void DisplayImage(Image im);
void NewEdgeDisplay();
void EndEdgeDisplay();
void DisplayEdge(float r1, float c1, float r2, float c2);

/* Following are from util.c */
void FatalError(char *msg);
void SysError(char *msg);
void *MallocPool(int size, int pool);
void FreeStoragePool(int pool);
float **AllocMatrix(int rows, int cols, int pool);
Image CreateImage(int rows, int cols, int pool);
void SolveLeastSquares(float *solution, int rows, int cols, float **jacobian,
		       float *errvec, float **sqarray);
void SolveLinearSystem(float *solution, float **sq, int size);
void FindKeypoints(Image image);

void writeImage(Image im, char *name);
void fastSIFT(DARY *image, vector<FD *> &features, Params *params);
#endif
