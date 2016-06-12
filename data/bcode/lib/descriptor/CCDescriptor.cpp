#include "feature.h"
#include "cc_base.h"
#include "cc_pca_base.h"
#include "../gauss_iir/gauss_iir.h"

static int CCSize=9;

void computeCC(DARY *patch,FeatureDescriptor *ds){
	DARY * simgn = new DARY(PATCH_SIZE,PATCH_SIZE);   
	int size=CCSize*CCSize;
	ds->allocVec(size);
	float *vec = ds->getVec();
	int step = PATCH_SIZE/CCSize;
	smooth(patch,simgn,step);
	for(int j=0,yi=step;j<CCSize;j++,yi+=step){
	  for(int i=0,xi=step;i<CCSize;i++,xi+=step){
	    vec[j*CCSize+i]=simgn->fel[yi][xi];
	  }
	}

	float sum=0;
	for(int j=0;j<size;j++){
	  sum+=vec[j];
	}
	
	for(int j=0;j<size;j++){
	  vec[j]-=sum;
	}
	
	delete simgn;
	//int cc_pca_size=81;
	//ds->pca(cc_pca_size,cc_pca_avg,cc_pca_base);	

} 

