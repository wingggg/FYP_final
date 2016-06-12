 
#ifndef _gauss_
#define _gauss_

#ifndef M_2PI
//#define M_PI  3.1415926537
#define M_2PI 6.2831853072
#endif

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef MAC_OS
#include <malloc.h>
#else
#include <malloc/malloc.h>
#endif
#include "../ImageContent/imageContent.h"
using namespace std;

#define  GAUSS_SIZE 100
static float gauss_100[100]={ 1.0000, 0.9997, 0.9989, 0.9974, 0.9954, 0.9929, 0.9898, 0.9861, 0.9819, 0.9771, 0.9718, 0.9660, 0.9597, 0.9529, 0.9455, 0.9377, 0.9295, 0.9207, 0.9116, 0.9020, 0.8920, 0.8816, 0.8708, 0.8597, 0.8483, 0.8365, 0.8244, 0.8120, 0.7993, 0.7864, 0.7733, 0.7599, 0.7463, 0.7326, 0.7187, 0.7047, 0.6905, 0.6763, 0.6619, 0.6475, 0.6331, 0.6186, 0.6041, 0.5896, 0.5751, 0.5607, 0.5463, 0.5320, 0.5177, 0.5036, 0.4895, 0.4756, 0.4618, 0.4482, 0.4347, 0.4214, 0.4082, 0.3952, 0.3825, 0.3699, 0.3575, 0.3454, 0.3334, 0.3217, 0.3103, 0.2991, 0.2881, 0.2773, 0.2668, 0.2566, 0.2466, 0.2369, 0.2274, 0.2182, 0.2092, 0.2005, 0.1920, 0.1838, 0.1758, 0.1681, 0.1606, 0.1534, 0.1464, 0.1397, 0.1332, 0.1269, 0.1209, 0.1150, 0.1094, 0.1040, 0.0988, 0.0939, 0.0891, 0.0845, 0.0801, 0.0759, 0.0719, 0.0680, 0.0643, 0.0608};

/************************************************************************
   structure used to store the coefficients nii+, nii-, dii+, dii- 
   and the normalisation factor 'scale'
   see [1] p. 9 - 11; p. 13, 14 
*************************************************************************/

#define  GAUSS_CUTOFF 3

void  dXY9(DARY* image_in, DARY* smooth_image);
void  dXX9(DARY* image_in, DARY* smooth_image);
void  dYY9(DARY* image_in, DARY* smooth_image);
void  smooth9(DARY* image_in, DARY* smooth_image);
void  dXY7(DARY* image_in, DARY* smooth_image);
void  dXX7(DARY* image_in, DARY* smooth_image);
void  dYY7(DARY* image_in, DARY* smooth_image);
void  dX6(DARY* image_in, DARY* smooth_image);
void  dY6(DARY* image_in, DARY* smooth_image);
void  dX4(DARY* image_in, DARY* smooth_image);
void  dY4(DARY* image_in, DARY* smooth_image);
void dX2(DARY* image_in, DARY* dximage);
void dY2(DARY* image_in, DARY* dyimage);
void HorConv3(DARY *image,  DARY *result);
void VerConv3(DARY *image,  DARY *result);
void grad2(DARY* image_in, DARY* dyimage);
void dXX3(DARY* image_in, DARY* dximage);
void dYY3(DARY* image_in, DARY* dyimage);
void dXX5(DARY* image_in, DARY* dximage);
void dYY5(DARY* image_in, DARY* dyimage);
void dXX_YY3(DARY* image_in, DARY* dyimage);
void  smooth3(DARY* image_in, DARY* smooth_image);
void  smooth5(DARY* image_in, DARY* smooth_image);
void  smoothSqrt(DARY* image_in, DARY* smooth_image);

void gradAngle(DARY *im, DARY *grad, DARY *ori);
void gradAngle(DARY *dx, DARY *dy, DARY *grad, DARY *ori);
void gradAngle(DARY *img, DARY *dx, DARY *dy, DARY *grad, DARY *ori);

 float smooth(int x, int y, DARY* image_in, float scale);
 float dX(int x, int y, DARY* image_in, float scale);
 float dY(int x, int y, DARY* image_in, float scale);
 float dXX(int x, int y, DARY* image_in, float scale);
 float dYY(int x, int y, DARY* image_in, float scale);
 float dXY(int x, int y, DARY* image_in, float scale);

 void smooth (DARY* image_in, DARY* out_image, float scale);
 void dX (DARY* image_in, DARY* out_image, float scale);
 void dY (DARY* image_in, DARY* out_image, float scale);
 void dXX (DARY* image_in, DARY* out_image, float scale);
 void dXY (DARY* image_in, DARY* out_image, float scale);
 void dYY (DARY* image_in, DARY* out_image, float scale);
 void dXX_YY (DARY* image_in, DARY* out_image, float scale);
 void dX (DARY* image_in,  DARY* image_out, float scalex, float scaley);	
 void dY (DARY* image_in,  DARY* image_out, float scalex, float scaley);	
 float smooth(int x, int y, DARY* image_in, float scalex, float scaley);
 void  smooth(DARY* image_in, DARY* smooth_image, float scalex, float scaley);
 void HorConvSqrt2(DARY *image,  DARY *result);
 void VerConvSqrt2( DARY *image,  DARY *result);


 float smoothf(int x, int y, DARY* image_in, float scale);
 float dXf(int x, int y, DARY* image_in, float scale);
 float dYf(int x, int y, DARY* image_in, float scale);
 float dXXf(int x, int y, DARY* image_in, float scale);
 float dXYf(int x, int y, DARY* image_in, float scale);
 float dYYf(int x, int y, DARY* image_in, float scale);
 float dXX_YYf(int x, int y, DARY* image_in, float scale);
 float dXXXf(int x, int y, DARY* image_in, float scale);
 float dXXYf(int x, int y, DARY* image_in, float scale);
 float dXYYf(int x, int y, DARY* image_in, float scale);
 float dYYYf(int x, int y, DARY* image_in, float scale);
 float dXXXXf(int x, int y, DARY* image_in, float scale);
 float dXXXYf(int x, int y, DARY* image_in, float scale);
 float dXXYYf(int x, int y, DARY* image_in, float scale);
 float dXYYYf(int x, int y, DARY* image_in, float scale);
 float dYYYYf(int x, int y, DARY* image_in, float scale);
void drawGauss(DARY* image_in,int x, int y, float scale);



#endif




