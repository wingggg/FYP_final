#include "feature.h"
#include "koen_base.h"
#include "koen_pca_base.h"



void computeKoen(DARY *img, FeatureDescriptor *ds){
	

        computeJLA(img,ds);
	float *vec = ds->getVec();  
	int nb_inv=12;
	double *val = new double[nb_inv];
	/* invariant f */
	/* invariant LwLw=fx*fx + fy*fy */
	val[ 0] = (vec[1]*vec[1]+vec[2]*vec[2]);
	float temp = sqrt(val[0]);
	temp=temp*temp*temp;
	/* invariant Lvv+Lww = fxx + fyy */
	val[ 1] = (vec[3]+vec[5]);
	
	/*invariant Lww=fx*fxx*fx+2*fx*fxy*fy+fy*fyy*fy */
	val[ 2] = (vec[1] * vec[1] * vec[3] + 
		   2 * vec[1] * vec[2] * vec[4] + vec[2] * vec[2] * vec[5]);
	
	/*invariant fxx*fxx+2*fxy*fxy+fyy*fyy */
	val[ 3] = (vec[3] * vec[3] + 2 * vec[4] * vec[4] + 
		   vec[5] * vec[5]);
	
	/* invariant of third order Lvvv (see p.321 paper of Romeny) */
	/* (fxxxfyfyfy-3fxxyfxfyfy+3fxyyfxfxfy-fyyyfxfxfx)/temp */
	val[ 4] = ((vec[6]*vec[2]*vec[2]*vec[2]-3*vec[7]*
                    vec[1]*vec[2]*vec[2]+ 3*vec[8]*vec[1]*vec[1]*vec[2]-vec[9]*
                    vec[1]*vec[1]*vec[1]))/temp ; 
	
	/* invariant of third order Lvvw (see p.321 paper of Romeny) */
	/* (fxxxfxfyfy+fxxy(-2fxfxfy+fyfyfy)+fxyy(-2fxfyfy+fxfxfx)+fyyyfxfxfy)/temp */
	val[ 5] = ((vec[6]*vec[1]*vec[2]*vec[2] + vec[7]*
                    (-2*vec[1]*vec[1]*vec[2]+vec[2]*vec[2]*vec[2]) + 
                    vec[8]*(-2*vec[1]*vec[2]*vec[2]+vec[1]*vec[1]*vec[1]) + 
                    vec[9]*vec[1]*vec[1]*vec[2]))/temp ;
    
	/* invariant of third order Lvww (see p.321 paper of Romeny) */
	/* (fxxy(-fxfxfx+2fxfyfy)+fxyy(-2fxfxfy+fyfyfy)-fyyyfxfyfy+fxxxfxfxfy)/temp */
	val[ 6] = ((vec[7]*(-vec[1]*vec[1]*vec[1]+ 2*vec[1]*vec[2]*vec[2])+
                    vec[8]*(-2*vec[1]*vec[1]*vec[2]+vec[2]*vec[2]*vec[2])-
                    vec[9]*vec[1]*vec[2]*vec[2]+vec[6]*vec[1]*vec[1]*vec[2]))/temp ;
	
	/* invariant of third order Lwww (see p.321 paper of Romeny) */
	/* (fxxxfxfxfx+3fxxyfxfxfy+3fxyyfxfyfy+fyyyfyfyfy)/temp */
	val[ 7] = ((vec[6]*vec[1]*vec[1]*vec[1]+
                    3*vec[7]*vec[1]*vec[1]*vec[2]+3*vec[8]*vec[1]*vec[2]*vec[2]+
                    vec[9]*vec[2]*vec[2]*vec[2]))/temp ;
	
	
	/* invariant of third order Lvw (see p.321 paper of Romeny) */
	/* (fxxfxfy+fxyfyfy-fxyfxfx-fyyfxfy)/temp */
	val[ 8] = (vec[3]*vec[1]*vec[2]+vec[4]*vec[2]*vec[2]-vec[4]*vec[1]*vec[1]-
		   vec[5]*vec[1]*vec[2])/temp ;

	/* invariant of third order Lvv (see p.321 paper of Romeny) */
	/* (fxxfyfy-2fxyfxfy+fyyfxfx)/temp */
	val[ 9] = (vec[3]*vec[2]*vec[2]-2*vec[4]*vec[1]*vec[2]+vec[5]*vec[1]*vec[1])/temp ;

	/* invariant of third order Lvvvv (see p.321 paper of Romeny) */
	/* (fxxxxfyfyfyfy-4fxyyyfxfxfxfy-4fxxxyfxfyfyfy+6fxxyyfxfxfyfy+fxfxfxfxfyyyy)/temp */
	val[ 10] = (vec[10]*vec[2]*vec[2]*vec[2]*vec[2]-4*vec[13]*vec[1]*vec[1]*vec[1]*vec[2]-
		    4*vec[11]*vec[1]*vec[2]*vec[2]*vec[2]+6*vec[12]*vec[1]*vec[1]*vec[2]*vec[2]+
		    vec[14]*vec[1]*vec[1]*vec[1]*vec[1])/temp ;

	/* invariant of third order Lwwww (see p.321 paper of Romeny) */
	/* (fxxxxfxfxfxfx+4fxyyyfxfyfyfy+4fxxxyfxfxfxfy+6fxxyyfxfxfyfy+fyfyfyfyfyyyy)/temp */
	val[ 11] = (vec[10]*vec[1]*vec[1]*vec[1]*vec[1]+4*vec[13]*vec[1]*vec[2]*vec[2]*vec[2]+
		    4*vec[11]*vec[1]*vec[1]*vec[1]*vec[2]+6*vec[12]*vec[1]*vec[1]*vec[2]*vec[2]+
		    vec[14]*vec[2]*vec[2]*vec[2]*vec[2])/temp ;


    ds->allocVec(nb_inv); 
    vec = ds->getVec();
    for(int i=0;i<nb_inv;i++)vec[i]=val[i];
    delete[] val;
    
    ds->changeBase(koen_base); 
    int koen_pca_size=12;
    //  ds->pca(koen_pca_size,  koen_pca_avg, koen_pca_base);	


}
