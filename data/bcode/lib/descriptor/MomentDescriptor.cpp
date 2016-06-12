#include "feature.h"
#include "mom_base.h"
#include "mom_pca_base.h"
#include "../gauss_iir/gauss_iir.h"

float computMoment000(DARY *img, int rad){
  float moment=0;
  for(int i=-rad;i<=rad;i++){
    for(int j=-rad;j<=rad;j++){
      if(patch_mask>0){
	moment++;
      }
    }
  }
  return moment;
}

float computMoment001(DARY *img, int rad){
  float moment=0;
  for(int i=-rad;i<=rad;i++){
    for(int j=-rad;j<=rad;j++){
      if(patch_mask>0){
	moment+=img->fel[rad+i][rad+j];
      }
    }
  }
  return moment;
}
float computMoment002(DARY *img, int rad){
  float moment=0;
  for(int i=-rad;i<=rad;i++){
    for(int j=-rad;j<=rad;j++){
      if(patch_mask>0){
	moment+=img->fel[rad+i][rad+j]*img->fel[rad+i][rad+j];
      }
    }
  }
  return moment;
}

float computMoment101(DARY *img, int rad){
  float moment=0;
  for(int i=-rad;i<=rad;i++){
    for(int j=-rad;j<=rad;j++){
      if(patch_mask>0){
	moment+=i*img->fel[rad+i][rad+j];
      }
    }
  }
  return moment;
}
float computMoment011(DARY *img, int rad){
  float moment=0;
  for(int i=-rad;i<=rad;i++){
    for(int j=-rad;j<=rad;j++){
      if(patch_mask>0){
	moment+=j*img->fel[rad+i][rad+j];
      }
    }
  }
  return moment;
}

float computMoment111(DARY *img, int rad){
  float moment=0;
  for(int i=-rad;i<=rad;i++){
    for(int j=-rad;j<=rad;j++){
      if(patch_mask>0){
	moment+=i*j*img->fel[rad+i][rad+j];
      }
    }
  }
  return moment;
}
float computMoment102(DARY *img, int rad){
  float moment=0;
  for(int i=-rad;i<=rad;i++){
    for(int j=-rad;j<=rad;j++){
      if(patch_mask>0){
	moment+=i*img->fel[rad+i][rad+j]*img->fel[rad+i][rad+j];
      }
    }
  }
  return moment;
}
float computMoment012(DARY *img, int rad){
  float moment=0;
  for(int i=-rad;i<=rad;i++){
    for(int j=-rad;j<=rad;j++){
      if(patch_mask>0){
	moment+=j*img->fel[rad+i][rad+j]*img->fel[rad+i][rad+j];
      }
    }
  }
  return moment;
}
float computMoment201(DARY *img, int rad){
  float moment=0;
  for(int i=-rad;i<=rad;i++){
    for(int j=-rad;j<=rad;j++){
      if(patch_mask>0){
	moment+=i*i*img->fel[rad+i][rad+j];
      }
    }
  }
  return moment;
}
float computMoment021(DARY *img, int rad){
  float moment=0;
  for(int i=-rad;i<=rad;i++){
    for(int j=-rad;j<=rad;j++){
      if(patch_mask>0){
	moment+=j*j*img->fel[rad+i][rad+j];
      }
    }
  }
  return moment;
}
float computMoment121(DARY *img, int rad){
  float moment=0;
  for(int i=-rad;i<=rad;i++){
    for(int j=-rad;j<=rad;j++){
      if(patch_mask>0){
	moment+=i*j*j*img->fel[rad+i][rad+j];
      }
    }
  }
  return moment;
}
float computMoment211(DARY *img, int rad){
  float moment=0;
  for(int i=-rad;i<=rad;i++){
    for(int j=-rad;j<=rad;j++){
      if(patch_mask>0){
	moment+=i*i*j*img->fel[rad+i][rad+j];
      }
    }
  }
  return moment;
}
float computMoment301(DARY *img, int rad){
  float moment=0;
  for(int i=-rad;i<=rad;i++){
    for(int j=-rad;j<=rad;j++){
      if(patch_mask>0){
	moment+=i*i*i*img->fel[rad+i][rad+j];
      }
    }
  }
  return moment;
}

float computMoment031(DARY *img, int rad){
  float moment=0;
  for(int i=-rad;i<=rad;i++){
    for(int j=-rad;j<=rad;j++){
      if(patch_mask>0){
	moment+=j*j*j*img->fel[rad+i][rad+j];
      }
    }
  }
  return moment;
}

float computMoment112(DARY *img, int rad){
  float moment=0;
  for(int i=-rad;i<=rad;i++){
    for(int j=-rad;j<=rad;j++){
      if(patch_mask>0){
	moment+=i*j*img->fel[rad+i][rad+j]*img->fel[rad+i][rad+j];
      }
    }
  }
  return moment;
}
float computMoment202(DARY *img, int rad){
  float moment=0;
  for(int i=-rad;i<=rad;i++){
    for(int j=-rad;j<=rad;j++){
      if(patch_mask>0){
	moment+=i*i*img->fel[rad+i][rad+j]*img->fel[rad+i][rad+j];
      }
    }
  }
  return moment;
}
float computMoment022(DARY *img, int rad){
  float moment=0;
  for(int i=-rad;i<=rad;i++){
    for(int j=-rad;j<=rad;j++){
      if(patch_mask>0){
	moment+=j*j*img->fel[rad+i][rad+j]*img->fel[rad+i][rad+j];
      }
    }
  }
  return moment;
}



void computeMoments(DARY *img, FeatureDescriptor *ds){
	DARY * dx = new DARY(PATCH_SIZE,PATCH_SIZE);   
	DARY * dy = new DARY(PATCH_SIZE,PATCH_SIZE);   
	dX2(img,dx);
	dY2(img,dy);
	int rad=dx->x()>>1;
	ds->allocVec(20);//EIndexSize*EIndexSize*OriSize;
	float *vec = ds->getVec();
	float norm=computMoment000(dx,rad);
	//vec[0]=computMoment001(dx,rad)/norm;
	//float mean=computMoment002(grad,rad);

	vec[0]=computMoment101(dx,rad)/norm;
	vec[1]=computMoment011(dx,rad)/norm;
 	vec[2]=computMoment102(dx,rad)/norm;
	vec[3]=computMoment012(dx,rad)/norm;
      
 	vec[4]=computMoment111(dx,rad)/norm;
	vec[5]=computMoment201(dx,rad)/norm;
	vec[6]=computMoment021(dx,rad)/norm;
	vec[7]=computMoment112(dx,rad)/norm;
	vec[8]=computMoment202(dx,rad)/norm;
	vec[9]=computMoment022(dx,rad)/norm;

	//norm=computMoment000(dy,rad);
	vec[10]=computMoment101(dy,rad)/norm;
	vec[11]=computMoment011(dy,rad)/norm;
 	vec[12]=computMoment102(dy,rad)/norm;
	vec[13]=computMoment012(dy,rad)/norm;
      
 	vec[14]=computMoment111(dy,rad)/norm;
	vec[15]=computMoment201(dy,rad)/norm;
	vec[16]=computMoment021(dy,rad)/norm;
	vec[17]=computMoment112(dy,rad)/norm;
	vec[18]=computMoment202(dy,rad)/norm;
	vec[19]=computMoment022(dy,rad)/norm;
	//vec[0]=0;
	delete dx;delete dy;
	ds->changeBase(mom_base);
	//int mom_pca_size=20;
	//ds->pca(mom_pca_size, mom_pca_avg, mom_pca_base);	

}

