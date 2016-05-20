#include "feature.h"
#include "gaussFilters.h"
#include "jla_base.h"
#include "jla_pca_base.h"

float Xf(DARY *img){
  float sum=0;
  for (int i=0;i<PATCH_SIZE;i++){
    for (int j=0;j<PATCH_SIZE;j++){
      sum+=img->fel[i][j]*LX[i][j];
    }
  }
  return sum;
}
float Yf(DARY *img){
  float sum=0;
  for (int i=0;i<PATCH_SIZE;i++){
    for (int j=0;j<PATCH_SIZE;j++){
      sum+=img->fel[i][j]*LX[j][i];
    }
  }
  return sum;
}
float XXf(DARY *img){
  float sum=0;
  for (int i=0;i<PATCH_SIZE;i++){
    for (int j=0;j<PATCH_SIZE;j++){
      sum+=img->fel[i][j]*LXX[i][j];
    }
  }
  return sum;
}
float YYf(DARY *img){
  float sum=0;
  for (int i=0;i<PATCH_SIZE;i++){
    for (int j=0;j<PATCH_SIZE;j++){
      sum+=img->fel[i][j]*LXX[j][i];
    }
  }
  return sum;
}
float XYf(DARY *img){
  float sum=0;
  for (int i=0;i<PATCH_SIZE;i++){
    for (int j=0;j<PATCH_SIZE;j++){
      sum+=img->fel[i][j]*LXY[i][j];
    }
  }
  return sum;
}
float XXXf(DARY *img){
  float sum=0;
  for (int i=0;i<PATCH_SIZE;i++){
    for (int j=0;j<PATCH_SIZE;j++){
      sum+=img->fel[i][j]*LXXX[i][j];
    }
  }
  return sum;
}
float YYYf(DARY *img){
  float sum=0;
  for (int i=0;i<PATCH_SIZE;i++){
    for (int j=0;j<PATCH_SIZE;j++){
      sum+=img->fel[i][j]*LXXX[j][i];
    }
  }
  return sum;
}
float XXYf(DARY *img){
  float sum=0;
  for (int i=0;i<PATCH_SIZE;i++){
    for (int j=0;j<PATCH_SIZE;j++){
      sum+=img->fel[i][j]*LXXY[i][j];
    }
  }
  return sum;
}
float XYYf(DARY *img){
  float sum=0;
  for (int i=0;i<PATCH_SIZE;i++){
    for (int j=0;j<PATCH_SIZE;j++){
      sum+=img->fel[i][j]*LXXY[j][i];
    }
  }
  return sum;
}
float XXXXf(DARY *img){
  float sum=0;
  for (int i=0;i<PATCH_SIZE;i++){
    for (int j=0;j<PATCH_SIZE;j++){
      sum+=img->fel[i][j]*LXXXX[i][j];
    }
  }
  return sum;
}
float YYYYf(DARY *img){
  float sum=0;
  for (int i=0;i<PATCH_SIZE;i++){
    for (int j=0;j<PATCH_SIZE;j++){
      sum+=img->fel[i][j]*LXXXX[j][i];
    }
  }
  return sum;
}
float XXXYf(DARY *img){
  float sum=0;
  for (int i=0;i<PATCH_SIZE;i++){
    for (int j=0;j<PATCH_SIZE;j++){
      sum+=img->fel[i][j]*LXXXY[i][j];
    }
  }
  return sum;
}
float XYYYf(DARY *img){
  float sum=0;
  for (int i=0;i<PATCH_SIZE;i++){
    for (int j=0;j<PATCH_SIZE;j++){
      sum+=img->fel[i][j]*LXXXY[j][i];
    }
  }
  return sum;
}
float XXYYf(DARY *img){
  float sum=0;
  for (int i=0;i<PATCH_SIZE;i++){
    for (int j=0;j<PATCH_SIZE;j++){
      sum+=img->fel[i][j]*LXXYY[i][j];
    }
  }
  return sum;
}


void computeJLA(DARY *img, FeatureDescriptor *ds){
	ds->allocVec(14);
	float *vec = ds->getVec();
	vec[0]=Xf(img);
	vec[1]=Yf(img);

	vec[2]=XXf(img);
	vec[3]=XYf(img);
	vec[4]=YYf(img);

	vec[5]=XXXf(img);
	vec[6]=XXYf(img);
	vec[7]=XYYf(img);
	vec[8]=YYYf(img);
	
	vec[9]=XXXXf(img);
	vec[10]=XXXYf(img);
	vec[11]=XXYYf(img);
	vec[12]=XYYYf(img);
	vec[13]=YYYYf(img);


	//int jla_pca_size=10;
	//ds->pca(jla_pca_size, jla_pca_avg, jla_pca_base); 
}

