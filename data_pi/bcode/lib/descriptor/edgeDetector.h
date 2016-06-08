#ifndef _edgeDetector_h_
#define _edgeDetector_h_

#include "../ImageContent/imageContent.h"
#include "../gauss_iir/gauss_iir.h"
using namespace std;
#include<vector>

void cannyEdges(DARY *dx,DARY *dy, DARY *edge,  float lower_threshold, float higher_threshold);
void cannyEdges(DARY *img, DARY *edge,  float scale, float lower_threshold, float higher_threshold);
void cannyEdges(DARY *dy, DARY *dx, DARY *grad, DARY *edge,  float lower_threshold,
		float higher_threshold);
void cannyEdges(DARY *img, DARY *edge,  float lower_threshold, float higher_threshold);


#endif
