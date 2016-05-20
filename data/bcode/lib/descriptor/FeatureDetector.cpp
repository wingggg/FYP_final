#include "feature.h"

#define ALPHA 0.06
#define SCALE_LAP_THRES 3
#define MIN_SCALE 2
#define MAX_SCALE 5
#define BORDER  3

/**
**/

void findAffineRegion(vector<DARY *> image,vector<DARY *> laplacian, vector<float> scale ,
                      vector<FeatureDescriptor*> &cor, int lmax);


void thresholdFeatures(vector<FeatureDescriptor*> &cor) {

    vector<FeatureDescriptor*> cor_tmp;
    for (int i=0;i<(int)cor.size();i++) {
        if (cor[i]->getFeatureness()>2000 && fabs(cor[i]->getLap())>15 && cor[i]->getExtr())
            cor_tmp.push_back(cor[i]);
    }
    cor.clear();
    cor=cor_tmp;
}

void findAffineRegion(vector<DARY *> image,vector<DARY *> laplacian, vector<float> scale ,
                      vector<FeatureDescriptor*> &cor, int lmax);

int check(int x, int y, int width, int height, int border) {
    if (x>border && y>border && x<width-border && y<height-border)return 1;
    else return 0;

}


/******************MULTI HESSIAN DETECTOR**************/

/******************MULTI HESSIAN DETECTOR**************/

/******************MULTI HESSIAN DETECTOR**************/

/******************MULTI HESSIAN DETECTOR**************/

/******************MULTI HESSIAN DETECTOR**************/

/******************MULTI HESSIAN DETECTOR**************/

/******************MULTI HESSIAN DETECTOR**************/

void max_detection(DARY *map,DARY * har11,DARY * har12,DARY * har22, vector<FeatureDescriptor*>&features, float threshold, float eratio, uint type, uint dradius) {

    FeatureDescriptor *cor;
    float act_pixel;
    float a,b,c,l1,l2,ea,el=1,thres=threshold;
    for (uint row =  BORDER; row+BORDER  < map->y(); row++) {
        for (uint col =  BORDER; col+BORDER  < map->x(); col++) {
            act_pixel=map->fel[row][col];
            if ( act_pixel > thres &&
                    act_pixel > map->fel[row-1][col-1] &&
                    act_pixel > map->fel[row-1][col] &&
                    act_pixel > map->fel[row-1][col+1] &&
                    act_pixel > map->fel[row][col-1] &&
                    act_pixel > map->fel[row][col+1] &&
                    act_pixel > map->fel[row+1][col-1] &&
                    act_pixel > map->fel[row+1][col] &&
                    act_pixel > map->fel[row+1][col+1]) {
                cor=new FeatureDescriptor((float)col,(float)row, 1.0, act_pixel);
                cor->setRadius(dradius);
                cor->setScale(dradius);
                cor->setType(type);
                if (type==DHARRIS || type==DSEDGE || type==DHARHES) {
                    a=har11->fel[row][col];
                    b=har12->fel[row][col];
                    c=har22->fel[row][col];
                    invSqrRoot(a,b,c,l2,l1,ea);
                    el=l2/l1;
                    cor->setMi(a,b,b,c,l1,el,ea);
                } else {
                    cor->setMi(1,0,0,1);
                }
                features.push_back(cor);
            }
        }

    }

}




void laplacian(DARY *sm, DARY * lap) {

    DARY *fxx = new DARY(sm->y(), sm->x(),FLOAT1);
    DARY *fyy = new DARY(sm->y(), sm->x(),FLOAT1);

    dXX9(sm,fxx);
    dYY9(sm,fyy);
    for (uint x=0;x<fxx->x();x++) {
        for (uint y=0;y<fxx->y();y++) {
            lap->fel[y][x]=fxx->fel[y][x]+fyy->fel[y][x];
        }
    }
    delete fxx;
    delete fyy;
    //lap->writePNG("lap.png");cout << "ci"<< endl;getchar();
}



/****** HESSIAN DETECTOR***********/
void hessian(DARY *image_in, DARY *hes, DARY *lap) {


    int col_nb, row_nb;
    float A, B, AB, determinant;
    float C, C2;
    row_nb = image_in->y();
    col_nb = image_in->x();

    //printf("- Starting to calculate gradients \n");
    DARY  *fxx= new DARY(row_nb,col_nb,FLOAT1);
    DARY  *fyy= new DARY(row_nb,col_nb,FLOAT1);
    DARY  *fxy= new DARY(row_nb,col_nb,FLOAT1);
    dXX9(image_in, fxx);//fxx->writePNG("fxx.png");
    dYY9(image_in, fyy);//fyy->writePNG("fyy.png");
    dXY7(image_in, fxy);//fxy->writePNG("fxy.png");

    int row, col;

    for ( row = 0; row < row_nb; row++)
        for ( col = 0; col < col_nb; col++)
        {
            /*        A = B = C = determinant = trace = 0.0;*/
            A = fxx->fel[row][col];
            B = fyy->fel[row][col];
            C = fxy->fel[row][col];
            C2=110*(C*C);// scaling factor to make equal amplitudes of fxx and fxy
            AB=(A*B);
            determinant = AB - (C2);
            lap->fel[row][col]=A+B;
            hes->fel[row][col] =(determinant);
        }
    //hes->writePNG("har.png");cout << "ci"<< endl;getchar();
    delete fxx;
    delete fyy;
    delete fxy;
}

void harris(DARY *img,DARY *har,DARY *har11,DARY *har12,DARY *har22) {

    int col_nb, row_nb;
    float A, B, C, determinant, trace, t1,t2;

    row_nb = img->y();
    col_nb = img->x();

    //   printf("- Starting to calculate gradients %d %d\n",col_nb, row_nb);
    DARY  *fx= new DARY(row_nb,col_nb,FLOAT1);
    dX6(img, fx);
    DARY  *fy= new DARY(row_nb,col_nb,FLOAT1);
    dY6(img, fy);
    DARY *fxy  = new DARY(row_nb, col_nb,FLOAT1);
    int row, col;
    for ( row = 0; row < row_nb; row++)
        for ( col = 0; col < col_nb; col++) {
            t1 = fx->fel[row][col];
            t2 = fy->fel[row][col];
            fx->fel[row][col] = t1*t1;
            fy->fel[row][col] = t2*t2;
            fxy->fel[row][col] = t1*t2;
        }

    smooth9(fx,har11);
    delete fx;
    smooth9(fy,har22);
    delete fy;
    smooth9(fxy,har12);
    delete fxy;

    for ( row = 0; row < row_nb; row++)
        for ( col = 0; col < col_nb; col++)
        {
            /*        A = B = C = determinant = trace = 0.0;*/
            A = har11->fel[row][col];
            B = har22->fel[row][col];
            C = har12->fel[row][col];
            determinant = A * B - (C*C);
            trace = A + B;
            har->fel[row][col] = (determinant - ALPHA * (trace*trace));
        }
}


void harris(DARY *dx,DARY *dy,DARY *har,DARY *har11,DARY *har12,DARY *har22) {

    int col_nb, row_nb;
    float A, B, C, determinant, trace, t1,t2;

    row_nb = dx->y();
    col_nb = dx->x();

    //   printf("- Starting to calculate gradients %d %d\n",col_nb, row_nb);
    DARY  *fx= new DARY(row_nb,col_nb,FLOAT1);
    DARY  *fy= new DARY(row_nb,col_nb,FLOAT1);
    DARY *fxy  = new DARY(row_nb, col_nb,FLOAT1);
    int row, col;
    for ( row = 0; row < row_nb; row++)
        for ( col = 0; col < col_nb; col++) {
            t1 = dx->fel[row][col];
            t2 = dy->fel[row][col];
            fx->fel[row][col] = t1*t1;
            fy->fel[row][col] = t2*t2;
            fxy->fel[row][col] = t1*t2;
        }

    smooth9(fx,har11);
    delete fx;
    smooth9(fy,har22);
    delete fy;
    smooth9(fxy,har12);
    delete fxy;

    for ( row = 0; row < row_nb; row++)
        for ( col = 0; col < col_nb; col++)
        {
            /*        A = B = C = determinant = trace = 0.0;*/
            A = har11->fel[row][col];
            B = har22->fel[row][col];
            C = har12->fel[row][col];
            determinant = A * B - (C*C);
            trace = A + B;
            har->fel[row][col] = (determinant) - ALPHA * (trace*trace);
        }
}




float inline interpScale(float a, float b, float c, int i, vector<float> scale) {
    float ds=interpPeak(a,b,c),sc;
    //return (1+ds*(scale[i]/scale[i-1]-1.0))*scale[i];
    return (pow((scale[i]/scale[i-1]),ds)*scale[i]);
    //return scale[i];
}

int getLapMax(vector<DARY *> lap, vector<float> scale, int level, int minlev, int maxlev, float x, float y,float &sc, int &extr) {

    vector<float> llap;
    float fx,fy;
    int flag=1;
    for (int i=0; i<(int)lap.size() && i<(level+maxlev) && flag;i++) {
        fx=x/scale[i];
        fy=y/scale[i];
        if (fx>BORDER && fx<(lap[i]->x()-BORDER) && fy>BORDER && fy<(lap[i]->y()-BORDER)) {
            llap.push_back(lap[i]->getValue(fx,fy));
            //cout << llap[llap.size()-1]<< " ";
        } else flag=0;
    }
    if (llap.size()<=level) {
        return -1;
    }

    //cout << endl;

    if (level-minlev<2)minlev=level-2;
    //cout << level-minlev<< " max " << level+maxlev<< endl;
    for (int i=level-minlev;i<(level+maxlev) && i<(int)(llap.size()-2);i++) {
        //cout << " li-1 "<< llap[i-1] << " "<< llap[i]<< " "<< llap[i+1]<< endl;
        //for local maximum or minimum of laplacian
        if (llap[i]>SCALE_LAP_THRES && llap[i-1]>= llap[i-2] && llap[i]> llap[i-1] && llap[i]>llap[i+1] && llap[i+1]>=llap[i+2]) {
            sc=interpScale(llap[i-1],llap[i],llap[i+1],i,scale);
            extr=1;
            llap.clear();
            return i;
        } else if (llap[i]<-SCALE_LAP_THRES && llap[i-1]<= llap[i-2] && llap[i]< llap[i-1] && llap[i]<llap[i+1] && llap[i+1]<=llap[i+2]) {
            sc=interpScale(llap[i-1],llap[i],llap[i+1],i,scale);
            extr=-1;
            llap.clear();
            return i;
        }
    }

    if (level<2)level=2;
    if (llap[level]>SCALE_LAP_THRES) {
        sc=scale[level];//sc=interpScale(llap[level-1],llap[level],llap[level+1],level,scale);
        extr=2;
        llap.clear();
        return level;
    } else if (llap[level]<-SCALE_LAP_THRES) {
        sc=scale[level];//interpScale(llap[level-1],llap[level],llap[level+1],level,scale);
        extr=-2;
        llap.clear();
        return level;
    } else {
        llap.clear();
        return -1;
    }
    return -1;
}

int getLapMax(vector<DARY *> lap, vector<float> scale, float step, int level, int minlev, int maxlev, float x, float y, vector<int> &levels, vector<float> &sc, vector<int> &extr, vector<float> &lapv) {

    vector<float> llap;
    float fx,fy;
    int flag=1;
    for (int i=0; i<(int)lap.size() && i<(level+maxlev) && flag;i++) {
        fx=x/scale[i];
        fy=y/scale[i];
        if (fx>BORDER && fx<(lap[i]->x()-BORDER) && fy>BORDER && fy<(lap[i]->y()-BORDER)) {
            llap.push_back(lap[i]->getValue(fx,fy));
            //cout << llap[llap.size()-1]<< " ";
        } else flag=0;
    }
    if (llap.size()<=level) {
        return -1;
    }

    float ds;
    if (level-minlev<1)minlev=level-1;
    //cout << " li-1 "<< (level-minlev) << " " << (level+maxlev)<< " ";
    //for(int i=level-minlev;i<(level+maxlev) && i<(int)(llap.size()-1);i++)cout << " "<< llap[i];cout << endl;
    for (int i=level-minlev;i<(level+maxlev) && i<(int)(llap.size()-1);i++) {
        //cout << " li-1 "<< i << " " << (level+maxlev)<< " " << llap[i-1] << " "<< llap[i]<< " "<< llap[i+1]<< endl;
        //for local maximum or minimum of laplacian
        if (llap[i]>SCALE_LAP_THRES && llap[i]> llap[i-1] && llap[i]>llap[i+1] ) {
//      if(llap[i]>SCALE_LAP_THRES && llap[i-1]>= llap[i-2] && llap[i]> llap[i-1] && llap[i]>llap[i+1] && llap[i+1]>=llap[i+2]){
            sc.push_back(pow(step,interpPeak(llap[i-1],llap[i],llap[i+1]))*scale[i]);
            //sc.push_back(interpScale(llap[i-1],llap[i],llap[i+1],i,scale));
            //cout << llap[i] << endl;
            extr.push_back(1);
            levels.push_back(i);
            lapv.push_back(llap[i]);
        } else if (llap[i]<-SCALE_LAP_THRES && llap[i]< llap[i-1] && llap[i]<llap[i+1]) {
//   }else if(llap[i]<-SCALE_LAP_THRES && llap[i-1]<= llap[i-2] && llap[i]< llap[i-1] && llap[i]<llap[i+1] && llap[i+1]<=llap[i+2]){
            sc.push_back(pow(step,interpPeak(llap[i-1],llap[i],llap[i+1]))*scale[i]);
            //sc.push_back(interpScale(llap[i-1],llap[i],llap[i+1],i,scale));
            extr.push_back(-1);
            lapv.push_back(llap[i]);
            //cout << lapv[lapv.size()-1] << endl;
            levels.push_back(i);
        }
    }

    if (levels.size()) {
        llap.clear();
        return 1;
    }

    if (level<2)level=2;
    if (llap.size()<=level && llap.size()>0)level=llap.size()-1;
    else if (llap.size()==0)return -1;
    if (llap[level]>SCALE_LAP_THRES) {
        sc.push_back(scale[level]);//sc=interpScale(llap[level-1],llap[level],llap[level+1],level,scale);
        extr.push_back(2);
        lapv.push_back(llap[level]);
        levels.push_back(level);
    } else if (llap[level]<-SCALE_LAP_THRES) {
        sc.push_back(scale[level]);//interpScale(llap[level-1],llap[level],llap[level+1],level,scale);
        extr.push_back(-2);
        lapv.push_back(llap[level]);
        levels.push_back(level);
    } else {
        llap.clear();
        return -1;
    }
    if (levels.size()) {
        llap.clear();
        return 1;
    }  else {
        llap.clear();
        return -1;
    }
}



void setMask(int x, int y, int r, DARY *mask) {
    int xs=((x-r)>=0)?(x-r):0;
    int ys=((y-r)>=0)?(y-r):0;
    int xe=((x+r)<(int)mask->x())?(x+r):(mask->x()-1);
    int ye=((y+r)<(int)mask->y())?(y+r):(mask->y()-1);

    //cout << xs << " " << ys << " " << xe << " " << ye << endl;

    for (int j=ys;j<=ye;j++) {
        for (int i=xs;i<=xe;i++) {
            mask->fel[j][i]=0;
        }
    }
}



void sedge_lap(vector<DARY *> edge,vector<DARY *> mask,vector<DARY *> har11,vector<DARY *> har12,vector<DARY *> har22, vector<DARY *> lap, vector<float> scale, vector<FeatureDescriptor*>&features, float threshold, int radius, float eratio, uint type, int dradius) {

    FeatureDescriptor *cor;
    float act_pixel,fx,fy,fcol,frow,lev,fscale;
    int level,sx=0,sy=0;
    vector<float> sc;
    vector<int> levels;
    vector<int> extr;
    vector<float> lapv;//(int)rint(GAUSS_CUTOFF*scale+2);
    //float sc;int extr;
    float a,b,c,l1,l2,ea,el=1,thres=threshold;
    int cbad=0,cgood=0;
    float step=scale[1]/scale[0];
    float level_red=2;
    float rescale=scale[level_red]/scale[0];

    for (uint i=0;i<edge.size()-1;i++) {
        cgood=0;
        for (uint row =  BORDER; row+BORDER  < edge[i]->y(); row++) {
            for (uint col =  BORDER; col+BORDER  < edge[i]->x(); col++) {
                act_pixel=edge[i]->fel[row][col];
                //     if(col==148 && row==61)
                //cout << "actpixel "<< act_pixel << "  "<< mask[i]->fel[row][col]<< endl;
                if ( act_pixel > thres && mask[i]->fel[row][col]>0 &&
                        ((type!=DMSER && type!=DSEDGE &&
                          act_pixel > edge[i]->fel[row-1][col-1] &&
                          act_pixel > edge[i]->fel[row-1][col] &&
                          act_pixel > edge[i]->fel[row-1][col+1] &&
                          act_pixel > edge[i]->fel[row][col-1] &&
                          act_pixel > edge[i]->fel[row][col+1] &&
                          act_pixel > edge[i]->fel[row+1][col-1] &&
                          act_pixel > edge[i]->fel[row+1][col] &&
                          act_pixel > edge[i]->fel[row+1][col+1]
                          /*&& edge[i]->fel[row][col+1] > thres &&
                          edge[i]->fel[row+1][col] > thres &&
                          edge[i]->fel[row][col-1] > thres &&
                          edge[i]->fel[row-1][col] > thres */
                         ) || type==DMSER ||type==DSEDGE)
                   ) {
                    fcol=(float)col+interpPeak(edge[i]->fel[row][col-1],act_pixel,edge[i]->fel[row][col+1]);
                    frow=(float)row+interpPeak(edge[i]->fel[row-1][col],act_pixel,edge[i]->fel[row+1][col]);
                    fx=scale[i]*fcol;//get coordinates at scale level 0
                    fy=scale[i]*frow;
                    levels.clear();
                    sc.clear();
                    extr.clear();
                    lapv.clear();
                    level=getLapMax(lap, scale, step, i, MIN_SCALE, MAX_SCALE, fx, fy, levels,sc, extr, lapv);//level=i+1;sc=scale[level];
                    //level=getLapMax(lap, scale, i, MIN_SCALE, MAX_SCALE, fx, fy,sc, extr);//level=i+1;sc=scale[level];
                    //cout << i << "  "  << cgood<< " " << fcol << " " << frow << " " << act_pixel << " " << levels.size()<< endl;cgood++;
                    if (mask[i]->fel[row][col]>0) {// && features.size()<10){
                        for (uint l=0;l<levels.size();l++) {
                            //if(level>0){
                            //	cout << sc.size()<< "   sc   " << sc[l] << " "<< extr[l] << " "<< levels[l] << " lapv "<< lapv[l] << endl;
                            //cout << sc.size()<< "fx "<< fx << " fy "<< fy <<  "   sc   " << sc[l] << " "<< extr[l] << " "<< levels[l] << " lapv "<< lapv[l] << endl;
//                            cor=new FeatureDescriptor(fx,fy, sc[l]*dradius/1.44, act_pixel);
			    lev=levels[l]-level_red;
			    lev=(lev<0)?0:lev;
			    fscale=dradius*sc[l]/rescale;
                            cor = new FeatureDescriptor(fx,fy, fscale, act_pixel);
                            cor->setIntLev(lev);
                            //cor->allocVec(192);
                            cor->setExtr(extr[l]);
                            cor->setRadius(dradius);
                            cor->setDerLev(i);
                            cor->setLap(lapv[l]);
//                            cor->setIntLev(levels[l]);
                            if (extr[l]>=0) {
                                cor->setType(type);
                            } else {
                                cor->setType(type | (type>>1));
                            }
                            setMask(col,row,radius,mask[i]);
                            if (type==DHARRIS || type==DSEDGE || type==DHARHES) {
                                a=har11[i]->fel[row][col];
                                b=har12[i]->fel[row][col];
                                c=har22[i]->fel[row][col];
                                invSqrRoot(a,b,c,l2,l1,ea);
                                el=l2/l1;
                                cor->setMi(a,b,b,c,l1,el,ea);
                            } else {
                                cor->setMi(1,0,0,1);
                            }
                            if (el>eratio || type==DMSER || type==DSEDGE)
                                features.push_back(cor);
                            else delete cor;
                        }
                    }

                }
            }

        }
//   cout << "good " << cgood << endl;  getchar();
    }
    levels.clear();
    sc.clear();
    extr.clear();
    lapv.clear();
}

void buildScaleSpace(vector<DARY *> sm, vector<DARY *> &grad, vector<DARY *> &gori) {
    for (uint i=0;i<sm.size();i++) {

        grad.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1,0.0));
        gori.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1,0.0));
        gradAngle(sm[i], grad[i], gori[i]);
        //cout << i <<" scale of " << sm.size()<< endl;
    }
}


void buildScaleSpace(DARY *img, vector<DARY *> &sm, vector<float> &sc, float step) {
    int sx=img->x(),sy=img->y();
    sc.push_back(1);
    sm.push_back(new DARY(sy,sx,FLOAT1));
    while ((sx/step)>12 && (sy/step)>12) {
        sx=(int)(0.5+sx/step);
        sy=(int)(0.5+sy/step);
        sm.push_back(new DARY(sy,sx,FLOAT1));
        sc.push_back(step*sc[sc.size()-1]);
    }
    
    sm[0]->set(img);
    //smooth5(img,sm[0]);
    sm[1]->scale(sm[0],step,step);
    float sc2=step*step;
    for (uint i=2;i<sm.size();i+=2) {
        //cout << i <<" scale of " << sm.size()<< endl;
        sm[i]->scale(sm[i-2],sc2,sc2);
        if ((i+1)<sm.size())
            sm[i+1]->scale(sm[i-1],sc2,sc2);
    }

}

void buildScaleSpaceFast(vector<DARY *> sm, vector<DARY *> &dx, vector<DARY *> &dy) {
    for (uint i=0;i<sm.size();i++) {		
        dx.push_back(new DARY(sm[i]->y(),sm[i]->x(),UCHAR3,0.0));
        dy.push_back(new DARY(sm[i]->y(),sm[i]->x(),UCHAR3,0.0));
//        dX2b(sm[i]->belr[0], dx[i]->belr[0], sm[i]->x(), sm[i]->y());
//		dX2b(sm[i]->belg[0], dx[i]->belg[0], sm[i]->x(), sm[i]->y());
//		dX2b(sm[i]->belb[0], dx[i]->belb[0], sm[i]->x(), sm[i]->y());
//        dX2b(sm[i]->belr[0], dy[i]->belr[0], sm[i]->x(), sm[i]->y());
//		dX2b(sm[i]->belg[0], dy[i]->belg[0], sm[i]->x(), sm[i]->y());
//		dX2b(sm[i]->belb[0], dy[i]->belb[0], sm[i]->x(), sm[i]->y());
		
        //cout << i <<" scale of " << sm.size()<< endl;
    }
}

void buildScaleSpaceFast(DARY *img, vector<DARY *> &sm, vector<float> &sc, float step) {
    int sx=img->x(),sy=img->y();
    sc.push_back(1);
    sm.push_back(new DARY(sy,sx,UCHAR3));
    while ((sx/step)>12 && (sy/step)>12) {
        sx=(int)(0.5+sx/step);
        sy=(int)(0.5+sy/step);
        sm.push_back(new DARY(sy,sx,UCHAR3));
        sc.push_back(step*sc[sc.size()-1]);
    }
	sm[0]->set(img);
    //smooth5(img,sm[0]);
    sm[1]->scale(sm[0],step,step);
    float sc2=step*step;
    for (uint i=2;i<sm.size();i+=2) {
        //cout << i <<" scale of " << sm.size()<< endl;
        sm[i]->scale(sm[i-2],sc2,sc2);
        if ((i+1)<sm.size())
            sm[i+1]->scale(sm[i-1],sc2,sc2);
    }
	
}


void saveScaleSpace(const char* filename, float step){
    vector<DARY *> sm;
    vector<DARY *> grad;
    vector<DARY *> gori;
    vector<float> sc;
	
	DARY *image = new DARY(filename);
	if (image->getType()==UCHAR1)image->convert(UCHAR3);
	
	if (image->getType()==UCHAR3) {
		image->RGB2opp();
		//image->convert(UCHAR3FLOAT1);
	}
	
	
	buildScaleSpaceFast(image, sm, sc, step);
    buildScaleSpaceFast(sm, grad, gori);
	
	
}

void loadScaleSpace(const char* filename, vector<DARY*> &grad, vector<DARY*> &gori,vector<float> &sc){ 
	
	
}


float computeOverlap(FD *fd1, FD *fd2, float overlap_thres) {//between 0-1, 0=perfect overlap, 1=no overlap
    //is in
    int xs=((fd1->getX()-fd1->getScale())>(fd2->getX()-fd2->getScale()))?(fd1->getX()-fd1->getScale()):(fd2->getX()-fd2->getScale());
    int ys=((fd1->getY()-fd1->getScale())>(fd2->getY()-fd2->getScale()))?(fd1->getY()-fd1->getScale()):(fd2->getY()-fd2->getScale());
    int xe=((fd1->getX()+fd1->getScale())<(fd2->getX()+fd2->getScale()))?(fd1->getX()+fd1->getScale()):(fd2->getX()+fd2->getScale());
    int ye=((fd1->getY()+fd1->getScale())<(fd2->getY()+fd2->getScale()))?(fd1->getY()+fd1->getScale()):(fd2->getY()+fd2->getScale());

    float s1=fd1->getScale();
    float s2=fd2->getScale();
    
    
    return 1.0-computeOverlap(fd1->getX()-s1, fd1->getY()-s1, fd1->getX()+s1,  fd1->getY()+s1, fd2->getX()-s2, fd2->getY()-s2, fd2->getX()+s2,  fd2->getY()+s2);
    
    


}

int selectFeature(vector<FeatureDescriptor*>&features,vector<uint> fsim) {
    //float ratio;
//  vector<FD*>feats;
    float maxlap=0;
    int maxl=-1;
    for (uint i=0;i<fsim.size();i++) {
        //  feats.push_back(features[fsim[i]]);
        //cout << features[fsim[i]]->getLap()<< endl;
        if (maxlap<fabs(features[fsim[i]]->getLap()) && fabs(features[fsim[i]]->getExtr())==1) {
            // if(maxlap<fabs(features[fsim[i]]->getLap())){
            maxlap=fabs(features[fsim[i]]->getLap());
            maxl=i;
        }
        //ratio=(float)pow(1.2,features[fsim[i]]->getIntLev())/(float)pow(1.2,features[fsim[i]]->getDerLev());
    }
    if(maxl<0){
      
      for (uint i=0;i<fsim.size();i++) {
	  if (maxlap<fabs(features[fsim[i]]->getLap())) {
            // if(maxlap<fabs(features[fsim[i]]->getLap())){
	      maxlap=fabs(features[fsim[i]]->getLap());
	      maxl=i;
	   }
      }
      
    }
    
    /*  cout << "write "<< maxl << " " << fsim.size() <<endl;
      DARY *image=new DARY(680,850,0.0);
      displayFeatures(image, feats, "feat.png", 255,DC);
      feats.clear();
      feats.push_back(features[fsim[maxl]]);
      displayFeatures(image, feats, "feat1.png", 100,DC);
      cout << "feats nb "<< feats.size()<< endl;
      getchar();*/

    return maxl;
}

void regionFilter(vector<FeatureDescriptor*>&features, float dist_thres) {
    vector<uint> fsim;
    DARY *dist=new DARY(features.size(),features.size(),FLOAT1,0.0);
// ofstream output("rat.txt");
    float lap,fet,extr;    
    for (uint i=0;i<features.size();i++) {
        //float rat=(float)pow(1.2,features[i]->getIntLev())/(float)pow(1.2,features[i]->getDerLev());
        //  output<< rat << " ";
        FD *f1=features[i];
        for (uint j=i+1;j<features.size();j++) {
            dist->fel[i][j]=(computeOverlap(f1,features[j],0.6));//0=perfect overlap
	}
    }
    //cout << features.size()<< endl;
    int *flag=new int[features.size()];
    bzero(flag,features.size()*sizeof(int));
    
    for (uint i=0;i<features.size();i++) {
        //float rat=(float)pow(1.2,features[i]->getIntLev())/(float)pow(1.2,features[i]->getDerLev());
        //  output<< rat << " ";
       // if(features[i]!=NULL){
	  FD *f1=features[i];
	  lap=fabs(f1->getLap());
	  fet=fabs(f1->getFeatureness());	
	  extr=fabs(f1->getExtr());	
	  //cout << i << endl;
	  for (uint j=i+1;j<features.size();j++) {	    
	    // if(features[i]!=NULL && features[j]!=NULL && dist->fel[i][j]<dist_thres){
	      if( dist->fel[i][j]<dist_thres){
		if(extr!=0 && features[j]->getExtr()!=0){				
		  if(lap > fabs(features[j]->getLap())){
		    //delete features[j];
		    flag[j]=1;		    
		  }else{
		    //delete features[i];
		    flag[i]=1;		    		    
		  }
		}else if(extr!=0){
		    //delete features[j];
		    flag[j]=1;		  
		}else if(features[j]->getExtr()!=0){
		    //delete features[i];
		    flag[i]=1;		    		    		  
		}else{		  
		  if(lap > fabs(features[j]->getLap())){
		    //delete features[j];
		    flag[j]=1;		    		  
		  }else{
		    //delete features[i];
		    flag[i]=1;		    		    		  		  
		  }		
		}
	      }
	  }
	//}
    }
    vector<FD*> fsel;
    for (uint i=0;i<features.size();i++) {
            if(flag[i]==0){
	      fsel.push_back(features[i]);
	    }  else {
	     delete  features[i];
	    }
    }
   features.clear();
   features=fsel;
    delete dist;
  //  cout << features.size()<< endl;
    
   return;
  //   delete dist;return;
// output.close();
    //vector<FD*> fsel;
    float *vect=new float[500];
    float med;
    int nb;
    cout << "rf " << features.size() << endl;

    for (uint i=0;i<features.size();i++) {
        if (features[i]!=NULL && features[i]->getX()>600 && features[i]->getX()<622 && features[i]->getY()>469 && features[i]->getY()<484) {
            cout << i << " ";features[i]->Cout();
            nb=0;
            fsim.push_back(i);
            vect[nb++]=features[i]->getScale();
            for (uint j=i+1;j<features.size();j++) {//find overlaping features to i
                if (features[j]!=NULL && dist->fel[i][j] < dist_thres) {
                    fsim.push_back(j);
                    if (nb<500)vect[nb++]=features[j]->getScale();
                    //if(round(features[i]->getX())==154 && round(features[i]->getY())==600)
		    //features[j]->Cout(1);
		    if(nb>1){cout << i << " " << j << "  "<< nb<< " " << dist->fel[i][j];features[j]->Cout(1);}
                }
            }
            cout << i << " " << nb << "  " << dist->fel[6607][6867]<< endl;
            if(fsim.size()>1){
	      med=selectFeature(features,fsim);//select best
	      if (med>=0){
		cout << "select " << fsim[med]<< endl; 
		fsel.push_back(features[fsim[med]]);
		for(uint j=0;j<fsim.size();j++){
		   if(med!=j)delete features[fsim[j]];
		   features[fsim[j]]=NULL;
		}
	      }
	    }else if(fsim.size()==1){
	      //cout << i << " " << fsim[0];features[fsim[0]]->Cout(1);
	      fsel.push_back(features[fsim[0]]);
	      features[fsim[0]]=NULL;
	      
	    }
	    fsim.clear();

        }
    }


    features.clear();
    features=fsel;
    delete dist;
    delete []vect;
}



/*****************used by fastSedge ********************************/
void computeDescriptors(vector<DARY*> grad, vector<DARY*> gori, vector<FeatureDescriptor*>&features, vector<float> sc, vector<float> logsc, uint DESC,
                        int middle_scale,  int noangle, int oriSize, int locSize, float max_gori) {

 //   cout << "computeDescriptors " << endl;
  
  int level=0;
    int xi,yi;
    float x,y,angle,scal,xf,yf,sf;
    if (features.size()==0)return;
    float radius=features[0]->getRadius();
    float sradius=((int)radius)>>1;
    initPatchMask((int)(2*radius+1));
    //DARY *patch = new DARY((2*radius+1),(2*radius+1),FLOAT1);
    FD *feat,*nfeat;
    uint fsize=features.size();
    unsigned long type;
    if (0) {//compute for small radius
        for (uint i=0;i<fsize;i++) {
            feat = new FD();
            feat->copy(features[i]);
            feat->setIntLev(feat->getDerLev());
            feat->setScale(sc[feat->getDerLev()]);
            feat->setRadius(sradius);
            //type=feat->getType();
            //feat->setType(type | (type>>2));
            features.push_back(feat);
        }
    }
    vector<float> angles;
    vector<FeatureDescriptor*> desc;
    uint nb_angle=noangle;
    if (noangle<2)nb_angle=1;
    float scale;
    for (uint i=0;i<features.size();i++) {
        //if(i==fsize)initPatchMask((int)(2*sradius+1));//for small radius
        feat = features[i];
        //feat->Cout(0);cout << " DSIZE " << feat->par[DSIZE]<< endl;
        sf = feat->getScale();
        //scal = sf/radius;
        //if(scal<1)scal=1;
        //level = getScaleIndex(scal, logsc, middle_scale);
        level=feat->getIntLev();
        radius=feat->getRadius();
 
        if (level<0){	  
	  level=0;
	}else if (level>=grad.size()) {
          level=grad.size()-1;
        }
	scale= sc[level];
	
	
        //cout<< "rad " << radius << " scale  " << sf  << " sc[level]  " << sc[level] << " level "<< level << " srad " << round(sf/(sc[level])) << endl;//getchar();
        radius=round(sf/scale);
        angle=feat->getAngle();
        xf=feat->getX();
        yf=feat->getY();
        x=(xf/scale);
        y=(yf/scale);

        feat->setX_Y(x,y);
        feat->setScale(radius);
        //patch->crop(grad[level],x,y);patch->writePNG("patch.png");cout << "OK " << endl;getchar();
        if (noangle!=0 && angle>M_2PI)
            computeHistAngle(x, y, radius, grad[level], gori[level], angles);
        else if (noangle==0) {
            angles.clear();
            angles.push_back(0.0);
        }
        if (angle<M_2PI)
            angles.insert(angles.begin(),angle);
        for (uint a=0;a<angles.size() && a<nb_angle;a++) {
            nfeat = new FD();
            nfeat->copy(feat);
            nfeat->setAngle(angles[a]);
            //cout << xf << " x " << yf << " sc " << sf << " scal " << scal << " level " << level << " xs " << x << " y "<< y << " r " << radius << endl;
            computeSiftOrig(grad[level], grad[level], grad[level], gori[level], nfeat,oriSize/*OriSize*/,locSize/*LocSize*/,max_gori);
            nfeat->setX_Y(xf,yf);
            nfeat->setScale(sf);
            //nfeat->setRadius(radius*sc[level]);
            //nfeat->setScale(radius*sf);
            desc.push_back(nfeat);
        }
	//cout << "done "<< endl;
        angles.clear();
    }
    deleteDescriptors(features);
    features=desc;

}


/*************************************************/
void computeDescriptors(DARY *img, vector<FeatureDescriptor*>&features, Params *params) {

    float step=params->getValue("scale_space_step.float");
    int noangle=params->getValue("descriptor_noangle.int");
    int color=params->getValue("descriptor_color.int");
    int dradius=(int)params->getValue("descriptor_radius.int");
    int OriSize=(int)params->getValue("OriSize.int");
    int LocSize=(int)params->getValue("LocSize.int");
    int OriSizeColor=(int)params->getValue("OriSizeColor.int");
    int LocSizeColor=(int)params->getValue("LocSizeColor.int");

    vector<DARY *> sm;
    vector<DARY *> grad;
    vector<DARY *> gori;
    vector<float> sc;
    vector<float> logsc;
    //img->writePNG("test.png");
    if (img->getType()!=UCHAR1FLOAT1 && img->getType()!=UCHAR3FLOAT1 ) {
        cout << "ERROR:  computeDescriptors, wrong image format " << img->getType() << endl;
        exit(1);
    }
    
     if((int)params->getValue("rect_ex.int")!=0){
     int sx = (int)params->getValue("rect_sx.int");
     int sy = (int)params->getValue("rect_sy.int");
     int ex = (int)params->getValue("rect_ex.int");
     int ey = (int)params->getValue("rect_ey.int");
     for(uint i=0;i<features.size();i++){
         if(features[i]->getX()<sx || features[i]->getX()>ex || features[i]->getY()<sy || features[i]->getY()>ey){
	    features.erase(features.begin()+i);
            i--;
	 }
     }    
   }

    
    
    
    
    buildScaleSpace(img, sm, sc, step);
    buildScaleSpace(sm, grad, gori);
    int middle_scale=0;
    for (uint i=0;i<sc.size();i++) {
        logsc.push_back(log(sc[i]));
        if (logsc[i]<0.00001)middle_scale=i;
    }
    int minsc=100;
    for (uint i=0;i<features.size();i++) {
        if (1 || features[i]->getIntLev()<0) {
            if (features[i]->getRadius()<1)
                features[i]->setRadius(dradius);
            float rad=features[i]->getRadius();

            int level = getScaleIndex(log(features[i]->getScale()/rad), logsc, middle_scale);
            //if(level<minsc)minsc=level;
            //if((features[i]->getType()&DMSER)!=0)
            // level=(level<=10)?0:(level-10);
            features[i]->setIntLev(level);
            //features[i]->Cout();
        }
    }
    //cout <<"minsc "<< minsc << endl;
 
    computeDescriptors( grad, gori, features, sc, logsc, DSIFT, middle_scale, noangle, OriSize, LocSize,M_2PI);
    // cout << "computing uniform 1" << endl;
    for (uint i=0;i<sm.size();i++) {
        //gori[i]->writePNG("grad.png");cout << "done " << endl;getchar();
        delete sm[i];
        delete grad[i];
        delete gori[i];
    }
    sm.clear();
    grad.clear();
    gori.clear();
    sc.clear();
    //*******************************************************************
    //*******************************************************************
    //*******************************************************************

    if (img->getType()==UCHAR3FLOAT1 && color) {//color
        //cout << "doing Color now"<< endl;
        //must copy image first
        //DARY *test = new DARY(img->y(),img->x(),UCHAR1);
        for (uint i=0;i<img->size();i++) {
            img->fel[0][i]=(float)img->belg[0][i];//copy RG channel to fel
            //test->bel[0][i]=img->belg[0][i];//copy RG channel to fel

        }
        //test->writePNG("rg.png");
        buildScaleSpace(img, sm, sc, step);
        buildScaleSpace(sm, grad, gori);
        vector<FeatureDescriptor*> rgfeats;
        for (uint i=0;i<features.size();i++) {
            FD *fd = new FD();
            fd->copy(features[i]);
            rgfeats.push_back(fd);
            //fd->Cout();getchar();
        }
        //cout << "noangle   "<< noangle << endl;
        computeDescriptors( grad, gori, rgfeats, sc, logsc, DSIFT, middle_scale, noangle, OriSizeColor, LocSizeColor, M_2PI);
        //for(uint i=0;i<grad.size();i++)grad[i]->set(1.0);
        //computeDescriptors( grad, sm, rgfeats, sc, logsc, DSIFT, middle_scale, noangle, 8/*OriSize*/, 2/*LocSize*/,255);

       // cout << " rg feats "<< rgfeats.size()<< endl;
        for (uint i=0;i<sm.size();i++) {
            //gori[i]->writePNG("grad.png");cout << "done " << endl;getchar();
            delete sm[i];
            delete grad[i];
            delete gori[i];
        }
        sm.clear();
        grad.clear();
        gori.clear();
        sc.clear();
        //*******************************************************************

        //must copy image first
        for (uint i=0;i<img->size();i++) {
            img->fel[0][i]=(float)img->belb[0][i];//copy BY channel to fel
            //test->bel[0][i]=img->belb[0][i];//copy RG channel to fel
        }
        //test->writePNG("by.png");
        buildScaleSpace(img, sm, sc, step);
        buildScaleSpace(sm, grad, gori);
        vector<FeatureDescriptor*> byfeats;
        //MUST copy features first
        for (uint i=0;i<features.size();i++) {
            FD *fd = new FD();
            fd->copy(features[i]);
            byfeats.push_back(fd);
        }
        computeDescriptors( grad, gori, byfeats, sc, logsc, DSIFT, middle_scale, noangle, OriSizeColor, LocSizeColor ,M_2PI);
        //for(uint i=0;i<grad.size();i++)grad[i]->set(1.0);
        //computeDescriptors( grad, sm, byfeats, sc, logsc, DSIFT, middle_scale, noangle, 8/*OriSize*/, 2/*LocSize*/,255);

       // cout << " by feats "<< byfeats.size()<< endl;

        for (uint i=0;i<sm.size();i++) {
            //gori[i]->writePNG("grad.png");cout << "done " << endl;getchar();
            delete sm[i];
            delete grad[i];
            delete gori[i];
        }
        sm.clear();
        grad.clear();
        gori.clear();
        sc.clear();
        //*******************************************************************
        //must combine all feats into features
        vector<FeatureDescriptor*> feats;
        for (uint i=0;i<features.size();i++) {
            FD *fd = new FD();
            fd->copy(features[i]);
            fd->allocVec(features[i]->getSize()+rgfeats[i]->getSize()+byfeats[i]->getSize());
            int s=features[i]->getSize();
            for (uint j=0;j<s;j++) {
                fd->vec[j]=features[i]->vec[j];
            }
            int s1=rgfeats[i]->getSize();
            for (uint j=0;j<s1;j++) {
                fd->vec[j+s]=(((int)rgfeats[i]->vec[j]));
            }
            int s2=byfeats[i]->getSize();
            for (uint j=0;j<s2;j++) {
                fd->vec[j+s+s1]=(((int)byfeats[i]->vec[j]));
            }
            feats.push_back(fd);
        }



        deleteDescriptors(features);
        deleteDescriptors(rgfeats);
        deleteDescriptors(byfeats);
        features.clear();
        features=feats;

    }
    
     logsc.clear();


}
/*************************************************/
/*****************used by fastSedge ********************************/
void ncomputeDescriptors( vector<DARY *> grad,vector<DARY *> gori,
                          vector<float> sc, vector<FeatureDescriptor*>&features, uint DESC, int noangle) {

    vector<float> logsc;
    int middle_scale=0;
    for (uint i=0;i<sc.size();i++) {
        logsc.push_back(log(sc[i]));
        if (logsc[i]<0.00001)middle_scale=i;
    }
    computeDescriptors(grad, gori, features, sc, logsc,  DESC,  middle_scale, noangle, 8/*OriSize*/, 4/*LocSize*/,M_2PI);
    logsc.clear();

}

/**************fastEdge, fastHarrisHessian, fastHessian, fastHarris***********************************/
void computeDescriptors( vector<DARY *> sm, vector<DARY *> dx, vector<DARY *> dy,
                         vector<float> sc, vector<FeatureDescriptor*>&features, uint DESC, int noangle) {


    int level=0;
    int xi,yi;
    float x,y,angle,scal;
    vector<FeatureDescriptor*> desc;
    if (features.size()==0)return;
    int psize=(int)features[0]->getRadius();
    initPatchMask(psize);
    DARY * dxpatch = new DARY(psize,psize,FLOAT1,0.0);
    DARY * dypatch = new DARY(psize,psize,FLOAT1,0.0);
    DARY * patch = new DARY(psize,psize,FLOAT1,0.0);
    float c_scale=floor(psize/2.0);

    float scale_factor=(psize/2.88);//2*scale/(step*step)=2*1.44
    FD *feat;
    //cout << "noangle " << noangle << endl;

    uint fsize=features.size();

    //for(uint i=0;i<fsize;i++){
    //feat = new FD();
    //feat->copy(features[i]);
    //feat->setIntLev(feat->getIntLev()+3);
    //feat->setScale(feat->getScale()*sc[3]);
    //features.push_back(feat);
    //}

    for (uint i=0;i<features.size();i++) {
        feat=features[i];
        level=feat->getIntLev()-2;    //scale/(step*step)=1.44
        angle=feat->getAngle();

        x=(feat->getX()/sc[level]);
        y=(feat->getY()/sc[level]);

        feat->setScale(feat->getScale()*scale_factor);//(patch_size/2)*scale/(step*step)

        if ((x+c_scale)<dx[level]->x() && (x-c_scale)>=0 && (y+c_scale)<dx[level]->y() && (y-c_scale)>0 ) {
            if ((DSIFT&DESC)==DSIFT) {
                xi=(int)(0.5+x);
                yi=(int)(0.5+y);
                dxpatch->set(0.0);
                dypatch->set(0.0);
                dxpatch->crop(dx[level],xi,yi);
                dypatch->crop(dy[level],xi,yi);
                //patch_mask->fel[0][0]=-0.00001;patch_mask->writePNG("patch_mask.png");
                //patch->set(0.0);patch->crop(sm[level],xi,yi);patch->writePNG("mask.png");cout << feat->getIntLev() << " "  << feat->getX() << " " << feat->getY() << " " << feat->getScale()<< endl;getchar();
                //if(features[i]->getScale()>7){dxpatch->writePNG("patch2.png");cout <<i <<  " OK 2 "<<features[i]->getX() << " " << features[i]->getY() <<" " << features[i]->getScale() <<  endl;getchar();}
                //if(!(i%100))cout << "\rdescriptor " << i << " of " << features.size()<< "    "<< flush;
                computeDescriptor( dxpatch, dypatch, DESC, feat, desc, noangle);
                //feat->Cout(128);

            } else if ((DESC&DESCRIPTOR)!=DCOLOR) {

                DARY *imgbs= normalizeAffine(sm[level], x, y, c_scale,
                                             angle, feat->getMi11(), feat->getMi12(),
                                             feat->getMi21(),feat->getMi22(),
                                             scal,dxpatch, 1.0);
                //patch->writePNG("patch1.png");cout <<i <<  " OK1 "<<features[i]->getX() << " " << features[i]->getY() <<  endl;getchar();
                //if(!(i%100))cout << "\raff_descriptor " << i << " of " << features.size()<< "    "<< flush;
                computeAffineDescriptor(imgbs, dxpatch, scal,DESC,feat,desc);
                delete imgbs;
            }
        } else {
            features.erase(features.begin()+i);
            i--;

        }

    }
    //cout << desc.size()<< endl;
    delete dxpatch;
    delete dypatch;
    //  cout << "OK2 " << (DESC&DESCRIPTOR)<< " " << DCOLOR << " " << DSIFT << endl;
    if ((DESC&DESCRIPTOR)!=DCOLOR) {
        deleteDescriptors(features);
        features=desc;
        //cout << "2OK" << features.size()<<  endl;
    }
}


void fastEdge(DARY* img, vector<FeatureDescriptor*>&features, Params *params) {
    float edgeHthres=params->getValue("edge_Hthreshold.float");
    float edgeLthres=params->getValue("edge_Lthreshold.float");
    float step=params->getValue("scale_space_step.float");
    uint DESC= (uint)params->getValue("feature_type.int");
    int aff=(int)params->getValue("affine_iterations.int");
    int noangle=(int)params->getValue("descriptor_noangle.int");
    int radius=(int)params->getValue("detector_mask_radius.int");
    int dradius=(int)params->getValue("descriptor_radius.int");

    fastEdge( img, features,  edgeHthres,  edgeLthres, step,  DESC,  aff, radius, noangle, dradius);
}

void fastEdge(DARY* img, vector<FeatureDescriptor*>&features, float edgeHthres, float edgeLthres,
              float step, uint DESC, int aff, int radius, int noangle, int dradius) {

    vector<DARY *> sm;
    vector<DARY *> lap;
    vector<DARY *> edge;
    vector<DARY *> dx;
    vector<DARY *> dy;
    vector<DARY *> grad;
    vector<DARY *> ori;
    vector<DARY *> mask;
    vector<float> sc;

    buildScaleSpace(img, sm, sc,step);
    cout << " sm " << sm.size()<< endl;
    for (uint i=0;i<sm.size();i++) {
        lap.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));
        laplacian(sm[i],lap[i]);
        dx.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
        dy.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
        if (i>=0 && i<sm.size()-1) {
            cout<< i << endl;
            grad.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
            ori.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
            gradAngle(sm[i],dx[i],dy[i],grad[i],ori[i]);
            edge.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1,0.0));//
            mask.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1,255.0));//
            cannyEdges(dx[i],dy[i], grad[i], edge[i], edgeLthres, edgeHthres);

            char nm[512];
            sprintf(nm,"sm%d.png",i);
            sm[i]->writePNG(nm);
            sprintf(nm,"grad%d.png",i);
            grad[i]->writePNG(nm);
            sprintf(nm,"edge%d.png",i);
            edge[i]->writePNG(nm);
            sprintf(nm,"lap%d.png",i);
            lap[i]->writePNG(nm);
            getchar();
        } else {
            grad.push_back(new DARY(1,1,FLOAT1,0.0));//
            ori.push_back(new DARY(1,1,FLOAT1));//
            edge.push_back(new DARY(1,1,FLOAT1,0.0));//
            mask.push_back(new DARY(1,1,FLOAT1));//
        }
    }
    sedge_lap(edge,mask,edge,edge,edge,lap,sc,features,edgeLthres,radius,1,DMSER,dradius);
    cout<< "edge2 " << features.size()<< endl;

    if ((DESC&DESCRIPTOR)!=0)computeDescriptors( sm, dx, dy, sc, features, DESC, noangle);


    for (uint i=0;i<sm.size();i++) {
        delete sm[i];
        delete lap[i];
        delete edge[i];
        delete dx[i];
        delete dy[i];
        delete grad[i];
        delete ori[i];
        delete mask[i];
    }
    sc.clear();
    sm.clear();
    lap.clear();
    edge.clear();
    dx.clear();
    dy.clear();
    grad.clear();
    ori.clear();
    mask.clear();

}



void fastSegm(DARY* img, vector<FeatureDescriptor*>&features, float edgeHthres, float edgeLthres, float step, uint DESC, int noangle) {

    vector<DARY *> sm;
    vector<float> sc;
    vector<FD*> desc;
    DARY *dx,*dy,*grad,*ori,*edge;
    //cout << "SEGM "<< endl;
    buildScaleSpace(img, sm, sc,step);
    for (uint i=0;i<sm.size();i+=2) {
        dx=new DARY(sm[i]->y(),sm[i]->x(),FLOAT1);//
        dy=new DARY(sm[i]->y(),sm[i]->x(),FLOAT1);//
        grad=new DARY(sm[i]->y(),sm[i]->x(),FLOAT1);//
        ori=new DARY(sm[i]->y(),sm[i]->x(),FLOAT1);//
        edge=new DARY(sm[i]->y(),sm[i]->x(),FLOAT1);//

        gradAngle(sm[i],dx,dy,grad,ori);
        cannySegments(dx,dy, grad,ori, edge, desc,  edgeLthres, edgeHthres);
        for (uint j=0;j<desc.size();j++) {
            //cout << desc[j]->getX() << "  "  << desc[j]->getY()<< endl;
            desc[j]->setX_Y(sc[i]*desc[j]->getX(),sc[i]*desc[j]->getY());
            desc[j]->setScale(sc[i]*desc[j]->getScale());
        }

        features.insert(features.end(),desc.begin(),desc.end());
        desc.clear();
        delete dx;
        delete dy;
        delete grad;
        delete ori;
        delete edge;
        //char nm[512];sprintf(nm,"sm%d.png",i);sm[i]->writePNG(nm);sprintf(nm,"grad%d.png",i);grad[i]->writePNG(nm);
        //sprintf(nm,"edge%d.png",i);edge[i]->writePNG(nm);sprintf(nm,"lap%d.png",i);lap[i]->writePNG(nm);
    }


    for (uint i=0;i<sm.size();i++) {
        delete sm[i];
    }
    sc.clear();
    sm.clear();
}

void harEdge(DARY* edge,DARY *har) {

    for (uint row = 0 ; row  < edge->y(); row++) {
        for (uint col =  0; col  < edge->x(); col++) {
            if (edge->fel[row][col]<=0 || har->fel[row][col]<=100)
                har->fel[row][col]=0;
        }
    }
}
void fastSegm(DARY* img, vector<FeatureDescriptor*>&features, Params *params) {

    float harthres=params->getValue("harris_threshold.float");
    float hesthres=params->getValue("hessian_threshold.float");
    float edgeHthres=params->getValue("edge_Hthreshold.float");
    float edgeLthres=params->getValue("edge_Lthreshold.float");
    float step=params->getValue("scale_space_step.float");
    uint DESC= (uint)params->getValue("feature_type.int");
    int aff=(int)params->getValue("affine_iterations.int");
    int noangle=(int)params->getValue("descriptor_noangle.int");
    int radius=(int)params->getValue("detector_mask_radius.int");
    int dradius=(int)params->getValue("descriptor_radius.int");
    float eigenratio=params->getValue("detector_eigenratio_threshold.float");

    fastSegm( img, features, edgeHthres, edgeLthres, step, DESC, noangle);

}
void corner_detector(DARY *image,DARY *labels, vector<FD *> &features, int mask, int type, float dradius, float thres) {

    if (image->getType()!=FLOAT1) {
        cout << "wrong image type in FeatureDetector:corner_detector" << endl;
        exit(0);
    }
    int width=labels->x();
    int height=labels->y();
    float maxval= (2*mask+1)*(2*mask+1);
//  DARY *map= new DARY(height,width,maxval/2.0);
    DARY *map= new DARY(height,width,FLOAT1,0.0);
    float diff;
    int drad=3;
    for (int j=mask;j<height-mask-1;j++) {
        for (int i=mask;i<width-mask-1;i++) {
            {//      if(1||labels->fel[j][i-drad]!=labels->fel[j][i+drad] || labels->fel[j-drad][i]!=labels->fel[j+drad][i] ||
                // labels->fel[j-drad][i-drad]!=labels->fel[j+drad][i+drad] || labels->fel[j+drad][i-drad]!=labels->fel[j-drad][i+drad] ){
                float lab=labels->fel[j][i];
                //float pix=image->fel[j][i];
                float nb=0,nb2=0;
                for (int m=-mask;m<=mask;m++) {
                    float *fel=labels->fel[j+m]+i;
                    //float *iel=image->fel[j+m]+i;
                    for (int n=-mask;n<=mask;n++) {
                        //if(lab==fel[n])nb++;
                        //else nb2++;
                        //diff=fabs(pix-iel[n]);
                        if (lab==fel[n])nb++;//=diff;
                        else nb2++;//=diff;
                    }
                }
                map->fel[j][i]=maxval-nb;//(nb<nb2)?nb:nb2;
            }
        }
    }

    DARY *smap= new DARY(map);
    smooth5(map,smap);
    FD *cor;
    vector<pair<float,int> > tmpVote;
    vector<FD *> feats;
    for (int j=mask;j<height-mask-1;j++) {
        for (int i=mask;i<width-mask-1;i++) {
            float val=smap->fel[j][i];
            if (val>thres && val>smap->fel[j-1][i-1] &&
                    val>smap->fel[j-1][i] &&
                    val>smap->fel[j-1][i+1] &&
                    val>smap->fel[j][i-1] &&
                    val>smap->fel[j][i+1] &&
                    val>smap->fel[j+1][i-1] &&
                    val>smap->fel[j+1][i] &&
                    val>smap->fel[j+1][i+1]) {
                cor=new FeatureDescriptor((float)i,(float)j, 1.0, val);
                cor->setRadius(dradius);
                cor->setScale(dradius);
                cor->setType(type);
                tmpVote.push_back(make_pair(val,feats.size()));
                feats.push_back(cor);

//         for(int c=-2;c<=2;c++){
//           map->fel[j+c][i+c]=255; map->fel[j-c][i+c]=255;
//         }

            }

        }
    }
    sort(tmpVote.begin(),tmpVote.end(),greater< pair<float,int> >());
    for (uint i=0;i<tmpVote.size();i++) {
        features.push_back(feats[tmpVote[i].second]);
    }
    tmpVote.clear();
    feats.clear();
    //smap->writePNG("smap.png");
// map->writePNG("map.png");
    delete map;
    delete smap;
}

void fastEGSF(DARY *image, vector<FD *> &features, Params *params) {
    //DARY *output= new DARY(image->y(),image->x(),"3uchar");
    DARY *output= new DARY(image->y(),image->x(),FLOAT1,0.0);
    graph_segmentation(image, output, params);
    //output->writePNG("segment.png");exit(0);
    int size = (int)params->getValue("susan_size.int");
    float thres = 10;//(int)params->getValue("susan_threshold.float");
    int dradius=(int)params->getValue("descriptor_radius.int");
    corner_detector(image, output, features, size, DSEGM, dradius, thres);
    delete output;
}


void exportAffVector(vector<RLERegion> &rle_vector,vector<FD*> &features, double factor, int extr, int dradius, int type)
{
    size_t i;

    sort(rle_vector.begin(), rle_vector.end());

    for (i=0; i < rle_vector.size(); i++)
    {
        const RLERegion *r = &rle_vector[i];
        double barX, barY, sumX2, sumY2, sumXY;
        RLE2Ellipse(r->rle, barX, barY, sumX2, sumXY, sumY2);
        Matrix C(2,2,0.0);
        C(1,1)=sumX2;
        C(1,2)=sumXY;
        C(2,1)=sumXY;
        C(2,2)=sumY2;

        Matrix Vi, V, D;
        C.svd(Vi,D,V);
        //cout << barX << " " << barY<< " "<< D(1,1)<< " " <<  D(2,2)<< "  " << a<< endl;getchar();
        D.tabMat[2][2]=factor*sqrt(D.tabMat[2][2]);
        D.tabMat[1][1]=factor*sqrt(D.tabMat[1][1]);
        float a=sqrt(D(2,2)*D(1,1));
        D.tabMat[2][2]/=a;
        D.tabMat[1][1]/=a;
        C=V*D*V.transpose();
        //C = C.inverse();
        FD *feat = new FD();
        float *par=feat->getPar();
        feat->setX_Y(barX,barY);
        feat->setScale(a);
        feat->setExtr(extr);
        feat->setRadius(dradius);

        feat->setType(type);
        //par[C_SCALE]=a;

        par[MI11]=C(1,1);
        par[MI12]=C(1,2);
        par[MI21]=C(2,1);
        par[MI22]=C(2,2);
        features.push_back(feat);
        //C=C*(factor*factor);
        //A = C.inverse();


        //fprintf(fid, "%g %g %g %g %g %d\n", barX, barY,
        //         A[0][0], A[0][1], A[1][1], r->margin);
    }
}

void fastMSERsingle(DARY* img, vector<FeatureDescriptor*>&features, Params *params) {
    ExtremaImage im;
    ExtremaParams p;
    size_t w, h, c;
    double scale_factor=3;
    int dradius=(int)params->getValue("descriptor_radius.int");
    p.preprocess=0;// = OptionInt("pre", 0, "image preprocessing type");
    p.max_area=0.01;// = OptionDouble("per", 0.01, "maximum relative area");
    p.min_size=20;// = OptionInt("ms", 30, "minimum size of output region");
    p.min_margin=10;// = OptionInt("mm", 10, "minimum margin");
    p.relative=0;// = OptionToggle("rel", 0, "use relative margins")!=0;
    p.verbose=0;// = OptionToggle("v", 0, "verbose output");
    p.debug=0;// = OptionInt("d", 0, "debug outputs");
    int type = DMSER;
    if (img->getType()==UCHAR1 || img->getType()==FLOAT1 || img->getType()==UCHAR1FLOAT1)
        im.channels = 1;
    else if (img->getType()==UCHAR3 || img->getType()==UCHAR3FLOAT1) {
        im.channels = 3;
        cout << "fastMSERsingle: wrong image type "<<endl;
        exit(1);
    }
    im.width = img->x();
    im.height = img->y();
    uint tsize=im.width*im.height*im.channels;
    im.data = new uchar[tsize];
    if (img->getType()==UCHAR1 || img->getType()==UCHAR1FLOAT1)memcpy(im.data,img->bel[0],tsize*sizeof(unsigned char));
    else if (img->getType()==FLOAT1)for (uint i=0;i<tsize;i++)im.data[i]=(uchar)img->fel[0][i];
    else {

        cout << "fastMSER format problem "<<  img->getType() << endl;
        exit(0);
    }
    //int ret = read_image(image_fname, im.data, w, h, c);
    if (!p.preprocess)
    {
        if (im.channels<3)
            p.preprocess = PREPROCESS_CHANNEL_none;
        else
            p.preprocess = PREPROCESS_CHANNEL_intensity;
    }
    //BoundaryExtrema result;
    //result = getBoundaryExtrema(p, im);
    RLEExtrema result;
    result = getRLEExtrema(p, im);
    //exportAffVector(out, result.MSERplus, scale_factor, 0);
    exportAffVector(result.MSERplus, features, dradius*scale_factor, 1,dradius, type);
    //exportAffVector(out, result.MSERmin, scale_factor,);
    exportAffVector(result.MSERmin, features, dradius*scale_factor,-1,dradius, (type | (type>>1)));

    delete []im.data;

}

void fastMSER(DARY* img, vector<FeatureDescriptor*>&features, Params *params) {

    float scaled=(int)params->getValue("mser_scaled.float");
    int mser_color=(int)params->getValue("mser_color.int");
    float reg_thes=params->getValue("region_filter.float");
    float desc_nb=params->getValue("descriptor_density_per_image.int");
    DARY *image = new DARY(img->y(),img->x(),UCHAR1);
    if (img->getType()==UCHAR3FLOAT1) {
        image->set(img->belr[0],img->y(),img->x());
    } else if (img->getType()==UCHAR1FLOAT1) {
        image->set(img->bel[0],img->y(),img->x());
    }
    fastMSERsingle(image, features, params);

    if (mser_color && (img->getType()==UCHAR3FLOAT1)) {
        image->set(img->belg[0],img->y(),img->x());
        fastMSERsingle(image, features, params);
        image->set(img->belb[0],img->y(),img->x());
        fastMSERsingle(image, features, params);
    }
    delete image;

    if (scaled>0.5 && scaled<2.5) {
        int xs=img->x()/scaled;
        int ys=img->y()/scaled;
        DARY *im = new DARY(ys,xs,img->getType());
        im->scale(img,scaled,scaled);

        image = new DARY(im->y(),im->x(),UCHAR1);
        if (im->getType()==UCHAR3FLOAT1) {
            image->set(im->belr[0],im->y(),im->x());
        } else if (im->getType()==UCHAR1FLOAT1) {
            image->set(im->bel[0],im->y(),im->x());
        }

        vector<FeatureDescriptor*>feats;
        fastMSERsingle(image, feats, params);
//   image->writePNG("mser1.png");
        if (mser_color && (img->getType()==UCHAR3FLOAT1)) {
            image->set(im->belg[0],im->y(),im->x());
            fastMSERsingle(image, feats, params);
//   image->writePNG("mser2.png");
            image->set(im->belb[0],im->y(),im->x());
            fastMSERsingle(image, feats, params);
//   image->writePNG("mser3.png");
        }
        delete image;
        //im->writePNG("mser.png");
        delete im;

        //cout << "OK1 "<< endl;

        for (uint i=0;i<feats.size();i++) {
            feats[i]->setX_Y(feats[i]->getX()*scaled,feats[i]->getY()*scaled);
            feats[i]->setScale(feats[i]->getScale()*scaled);
            //feats[i]->Cout();
        }
        features.insert(features.end(),feats.begin(),feats.end());
        feats.clear();
    }

    if (reg_thes!=0)regionFilter(features, reg_thes);
//thresFilter(features, desc_nb, mserthres);


    uint DESC=(uint)params->getValue("feature_type.int");
    if ((DESC&DESCRIPTOR)!=0) {
        computeDescriptors(img, features, params);
    }

}



void fastSEdge(DARY* img, vector<FeatureDescriptor*>&features, Params *params) {

    float harthres=params->getValue("harris_threshold.float");
    float hesthres=params->getValue("hessian_threshold.float");
    float edgeHthres=params->getValue("edge_Hthreshold.float");
    float edgeLthres=params->getValue("edge_Lthreshold.float");
    float step=params->getValue("scale_space_step.float");
    uint DESC= (uint)params->getValue("feature_type.int");
    int aff=(int)params->getValue("affine_iterations.int");
    int noangle=(int)params->getValue("descriptor_noangle.int");
    int radius=(int)params->getValue("detector_mask_radius.int");
    int dradius=(int)params->getValue("descriptor_radius.int");
    float eigenratio=params->getValue("detector_eigenratio_threshold.float");
    fastSEdge( img, features, harthres,  hesthres, edgeHthres, edgeLthres, step, DESC, aff, radius, eigenratio, noangle, dradius);
}

void fastSEdge(DARY* img, vector<FeatureDescriptor*>&features,float harthres, float hesthres,float edgeHthres,float edgeLthres,float step,uint DESC,int aff, int radius, float eigenratio, int noangle, int dradius) {

    vector<DARY *> sm;
    vector<DARY *> lap;
    vector<DARY *> edge;
    vector<DARY *> hes;
    vector<DARY *> har;
    vector<DARY *> har11;
    vector<DARY *> har12;
    vector<DARY *> har22;
    vector<DARY *> dx;
    vector<DARY *> dy;
    vector<DARY *> grad;
    vector<DARY *> ori;
    vector<DARY *> mask;
    vector<float> sc;
    vector<float> logsc;
    vector<FD*> segmdesc;
    vector<FD*> mserdesc;
    vector<FD*> tmpdesc;
    buildScaleSpace(img, sm, sc, step);
    //img->writePNG("median.png");
// cout << "scale space "<< sm.size() << endl;
    for (uint i=0;i<sc.size();i++) {
        logsc.push_back(log(sc[i]));
    }
    //laplacian(sm[0],lap[0]);
    for (uint i=0;i<sm.size();i++) {
        lap.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));
        hes.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
        hessian(sm[i],hes[i],lap[i]);//
        dx.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
        dy.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
        grad.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
        ori.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
        gradAngle(sm[i],dx[i],dy[i],grad[i],ori[i]);
        mask.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1,255.0));//
//cout << i << " " << sm[i]->y() << endl;
        edge.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1,0.0));//
        //cannyEdges(dx[i],dy[i], grad[i], edge[i], edgeLthres, edgeHthres);
        har.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
        har11.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
        har12.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
        har22.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
        harris(dx[i],dy[i],har[i],har11[i],har12[i],har22[i]);

        if (!(i%4) && sm[i]->y()>50 && sm[i]->x()>50) {
            cannySegments(dx[i], dy[i], grad[i], ori[i], edge[i], tmpdesc, edgeLthres, edgeHthres);

            //edge[i]->writePNG("edge.png");cout << "OK edge " << endl; getchar();
            /*	for(uint j=0;j<tmpdesc.size();j++){
            	  //cout << desc[j]->getX() << "  "  << desc[j]->getY()<< endl;
            	  tmpdesc[j]->setX_Y(sc[i]*tmpdesc[j]->getX(),sc[i]*tmpdesc[j]->getY());
            	  tmpdesc[j]->setScale(sc[i]*tmpdesc[j]->getScale());
            	}
            	segmdesc.insert(segmdesc.end(),tmpdesc.begin(),tmpdesc.end());
            	tmpdesc.clear();*/
        }


        //harEdge(edge[i],har[i]);
        //char nm[512];sprintf(nm,"sm%d.png",i);sm[i]->writePNG(nm);//sprintf(nm,"grad%d.png",i);grad[i]->writePNG(nm);
        //sprintf(nm,"edge%d.png",i);edge[i]->writePNG(nm);//sprintf(nm,"lap%d.png",i);lap[i]->writePNG(nm);
        //sprintf(nm,"har%d.png",i);har[i]->writePNG(nm);
        //sprintf(nm,"hes%d.png",i);hes[i]->writePNG(nm);
        //sprintf(nm,"har2%d.png",i);har[i]->writePNG(nm);

    }
    //cout << " space maps "<< endl;

    cout << "radius " << radius << "  eigenratio " << eigenratio << "  harthres "  << harthres << "  hesthres "<< hesthres << endl;
    sedge_lap(hes,mask,har11,har12,har22,lap,sc,features,hesthres,radius,eigenratio,DHARHES,dradius);
    cout<< "hes " << features.size()<< endl;
    sedge_lap(har,mask,har11,har12,har22,lap,sc,features,harthres,radius/*radius*/,eigenratio,DHARRIS,dradius);
    cout<< "har " << features.size()<< endl;
    //cout << " space segde "<< endl;

    sedge_lap(edge,mask,har11,har12,har22,lap,sc,features,100,radius<<1/*radius*/,eigenratio,DSEDGE,dradius);
    cout<< "edge " << features.size()<< endl;


    /*  for(uint i=0;i<sm.size();i++){
        if((!(i%2) || i==0) && sm[i]->x()>50);
          fastMSER(sm[i], lap, mserdesc ,sc,logsc, i,  DMSER);

      }*/
    //features.clear();
    //features.insert(features.end(),mserdesc.begin(),mserdesc.end());
    //noangle=1;
    //computeDescriptors( sm, dx, dy, sc, features, DESC, noangle);



    if ((DESC&DESCRIPTOR)!=0)ncomputeDescriptors(grad, ori,  sc, features, DESC, noangle);


    //  cout<< "mser " << features.size()<< endl;
    //for(uint i=0;i<features.size();i++){
    //  cout << i<< " pos " << features[i]->getType()<< endl;
    // }//if(features[i]->getV(0)<0 || features[i]->getV(0)>1000){features[i]->Cout(5);getchar();}

    //features.insert(features.end(),segmdesc.begin(),segmdesc.end());
    cout << "sift points "<< features.size()<< endl;

    //displayFeatures(sm[0], features, "disp.png", 255, 0);getchar();

    for (uint i=0;i<sm.size();i++) {
        delete sm[i];
        delete lap[i];
        delete dx[i];
        delete dy[i];
        delete grad[i];
        delete ori[i];
        delete har[i];
        delete har11[i];
        delete har12[i];
        delete har22[i];
        delete hes[i];
        delete mask[i];
        delete edge[i];
    }
    sc.clear();
    logsc.clear();
    sm.clear();
    lap.clear();
    edge.clear();
    dx.clear();
    dy.clear();
    grad.clear();
    ori.clear();

    hes.clear();
    har.clear();
    har11.clear();
    har12.clear();
    har22.clear();
    mask.clear();
    segmdesc.clear();
    mserdesc.clear();
}

void thresFilter(vector<FeatureDescriptor*>&features, int desc_nb, float threshold) {
    if (features.size()<desc_nb)return;
    vector<pair<float,int> > tmpVote;

    for (uint i=0;i<features.size();i++) {
        tmpVote.push_back(make_pair(features[i]->getFeatureness(),i));
    }
    sort(tmpVote.begin(),tmpVote.end(),greater< pair<float,int> >());
    vector<FeatureDescriptor*> feats;
    for (uint i=0;i<tmpVote.size();i++) {
        feats.push_back(features[tmpVote[i].second]);
        //feats[i]->Cout();
    }
    tmpVote.clear();
    features.clear();
    features.insert(features.begin(),feats.begin(),feats.begin()+desc_nb);
    for (uint i=desc_nb;i<feats.size();i++) {
        delete feats[i];
        //feats[i]->Cout();
    }
    feats.clear();


    return;

    /********************************************************************************/
    if (features.size()<desc_nb) {
        return;
    } else if (features.size()>(3.5*desc_nb)) {
        threshold=3.0*threshold;
    } else if (features.size()>(2.0*desc_nb)) {
        threshold=threshold;
    } else if (features.size()>(1.5*desc_nb)) {
        threshold=0.75*threshold;
    }
    //cout <<"thres " <<  threshold << endl;
    for (uint i=0;i<features.size();i++) {
        if (features[i]->getFeatureness()>threshold) {
            feats.push_back(features[i]);
        } else {
            delete features[i];
        }
    }
    features.clear();
    features=feats;
}

void fastHarrisHessian(DARY *img, vector<FeatureDescriptor*>&corners, Params *params) {
    float harthres=params->getValue("harris_threshold.float");
    float hesthres=params->getValue("hessian_threshold.float");
    float step=params->getValue("scale_space_step.float");
    uint DESC= (uint)params->getValue("feature_type.int");
    int aff=(int)params->getValue("affine_iterations.int");
    int noangle=(int)params->getValue("descriptor_noangle.int");
    int dradius=(int)params->getValue("descriptor_radius.int");
    float reg_thes=params->getValue("region_filter.float");
    float desc_nb=params->getValue("descriptor_density_per_image.int");
    fastHarrisHessian(img, corners,  harthres,  hesthres,  step,  DESC, reg_thes, aff,  noangle, dradius,desc_nb);
    if ((DESC&DESCRIPTOR)!=0)
        computeDescriptors(img, corners, params);
}

void fastHarrisHessian(DARY *img, vector<FeatureDescriptor*>&corners, float harthres, float hesthres, float step, uint DESC, float reg_thes, int aff, int noangle, int dradius, int desc_nb) {

    vector<DARY *> sm;
    vector<DARY *> lap;
    vector<DARY *> hes;
    vector<float> sc;
    vector<DARY *> dx;
    vector<DARY *> dy;
    vector<DARY *> har;
    vector<DARY *> har11;
    vector<DARY *> har12;
    vector<DARY *> har22;
    vector<DARY *> mask;

    buildScaleSpace(img, sm, sc,step);
    for (uint i=0;i<sm.size();i++) {
        lap.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
        hes.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
        hessian(sm[i],hes[i],lap[i]);//
        dx.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
        dy.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
        dX2(sm[i], dx[i]);
        dY2(sm[i], dy[i]);
        //dX6(sm[i], dx[i]);dY6(sm[i], dy[i]);

//   if(i>0 && i< sm.size()-1){
        har.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));
        har11.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));
        har12.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));
        har22.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));
        harris(dx[i],dy[i],har[i],har11[i],har12[i],har22[i]);
        mask.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1,255.0));//
        /*    }else {
              har.push_back(new DARY(1,1));
              har11.push_back(new DARY(1,1));
              har12.push_back(new DARY(1,1));
              har22.push_back(new DARY(1,1));
              mask.push_back(new DARY(1,1));
            }*/
        //char name[512];sprintf(name,"hes%d.png",i+2);hes[i]->writePNG(name);
        //cout <<i << " "<< sc[i] << endl;
    }
    //sedge_lap(har,mask,har11,har12,har22,lap,sc,corners,harthres,1/*radius*/,0.01/*eigenratio*/,DHARRIS,dradius);
    //if (reg_thes!=0)regionFilter(corners, reg_thes);//cout<< "har " << corners.size()<< endl;
    vector<FD*> feats;
    sedge_lap(hes,mask,har11,har12,har22,lap,sc,feats,hesthres,1/*radius*/,0.01/*eigenratio*/,DHARHES,dradius);
    if (reg_thes!=0)regionFilter(feats, reg_thes); // cout<< "hes " << feats.size()<< endl;

    desc_nb=(img->x()*img->y())/desc_nb;
    uint tot=corners.size()+feats.size();
    if (tot>0) {
        uint prop=(desc_nb*corners.size())/tot;
        uint prop2=(desc_nb*feats.size())/tot;
       // thresFilter(corners, prop, harthres);
        //cout<< "har " << corners.size()<< " " << prop<< endl;
        thresFilter(feats, prop2, hesthres);
       // cout<< "hes " << feats.size()<< " " << desc_nb<<  " " << prop2<<endl;
    }

    corners.insert(corners.end(),feats.begin(),feats.end());

    if (aff>1)findAffineRegion(sm, lap, sc, corners, aff);


    //cout << "Harris-Hessian-Laplace(affine) interest points "<< corners.size()<< " " << harthres << endl;



    //corners.clear();
    //fastMSER(img, lap, corners ,sc, 0, DMSER);

    //  deleteDescriptors(corners);
    /* if((DESC&DESCRIPTOR)!=0){
       vector<DARY *> grad;
       vector<DARY *> ori;
       for(uint i=0;i<sm.size();i++){
         grad.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
         ori.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
         gradAngle(dx[i],dy[i],grad[i],ori[i]);
       }
       ncomputeDescriptors(grad, ori,  sc, corners, DESC, noangle);
       //computeDescriptors( sm, dx, dy,  sc, corners, DESC, noangle);
       for(uint i=0;i<sm.size();i++){delete grad[i];delete ori[i];}
       grad.clear();ori.clear();
       }
    */

    for (uint i=0;i<sm.size();i++) {
        delete sm[i];
        delete lap[i];
        delete hes[i];
        delete har[i];
        delete har11[i];
        delete har12[i];
        delete har22[i];
        delete mask[i];
        delete dx[i];
        delete dy[i];
    }
    sc.clear();
    sm.clear();
    lap.clear();
    har.clear();
    hes.clear();
    har11.clear();
    har12.clear();
    har22.clear();
    hes.clear();
    dx.clear();
    dy.clear();
    mask.clear();
}

void fastUniform(DARY *img, vector<FeatureDescriptor*>&features, Params *params) {
 // cout << "computing uniform " << endl;
    float number=params->getValue("uniform_number.int");
    int scale_levels=params->getValue("uniform_scale_levels.int");
    int dradius=(int)params->getValue("descriptor_radius.int");
    float step=2;//params->getValue("scale_space_step.float");
     uint DESC= (uint)params->getValue("feature_type.int");
    
    float ratio = (float)img->x()/(float)img->y();
    int nb_per_width = 1+sqrt((number*ratio)/scale_levels);
    int nb_per_height = 1+sqrt(number/(scale_levels*ratio));
    vector<float> sc;
    int stepwidth;
    int stepheight;
    float sr;
    FD *fd;
    sc.push_back(1);
    for(uint i=1;i<scale_levels;i++)sc.push_back(step*sc[sc.size()-1]);       
    for(uint i=0;i<sc.size();i++){
      stepwidth=(img->x()/sc[i])/nb_per_width;
      stepheight=(img->y()/sc[i])/nb_per_height;
      sr=(stepwidth<stepheight)?stepwidth:stepheight;   
      //cout<< stepwidth << "   " << stepheight<<  "  "<<  sc[i] << "  " << sr  << endl;
         for(uint y=1;y<nb_per_height;y++){
           for(uint x=1;x<nb_per_width;x++){
	      fd=new FD();
	      fd->setX_Y(sc[i]*x*stepwidth,sc[i]*y*stepheight);
	      fd->setRadius(dradius);	      
	      fd->setScale(sr);
	      fd->setAngle(2*M_2PI);
	      fd->setType(DHESSIAN);
	      features.push_back(fd);
	      //fd->Cout();
	   }
	 }    
    }
  cout << "Uniform interest points " << features.size()  << endl;
	  
    if ((DESC&DESCRIPTOR)!=0)
        computeDescriptors(img, features, params);
   
  
}

void fastHessian(DARY *img, vector<FeatureDescriptor*>&features, Params *params) {
    float threshold=params->getValue("hessian_threshold.float");
    float step=params->getValue("scale_space_step.float");
    uint DESC= (uint)params->getValue("feature_type.int");
    int aff=(int)params->getValue("affine_iterations.int");
    int noangle=(int)params->getValue("descriptor_noangle.int");
    int dradius=(int)params->getValue("descriptor_radius.int");
    float reg_thes=params->getValue("region_filter.float");
    float desc_nb=params->getValue("descriptor_density_per_image.int");
    int hesdet=params->getValue("hessian_det.int");

    fastHessian(img, features,  threshold,  step,  DESC,  reg_thes, aff,  noangle, dradius,desc_nb, hesdet);
    if ((DESC&DESCRIPTOR)!=0)
        computeDescriptors(img, features, params);
}

void fastHessian(DARY *img, vector<FeatureDescriptor*>&features, float threshold, float step, uint DESC, float reg_thes, int aff, int noangle, int dradius, int desc_nb, int hesdet) {

    vector<DARY *> sm;
    vector<DARY *> lap;
    vector<DARY *> hes;
    vector<float>  sc;
    vector<DARY *> dx;
    vector<DARY *> dy;
    vector<DARY *> mask;

    buildScaleSpace(img, sm, sc,step);
    for (uint i=0;i<sm.size();i++) {
        lap.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
        hes.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
        hessian(sm[i],hes[i],lap[i]);//
        dx.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
        dy.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
        dX6(sm[i], dx[i]);
        dY6(sm[i], dy[i]);
        mask.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1,255.0));//
        //char name[512];sprintf(name,"hes%d.png",i);hes[i]->writePNG(name);sprintf(name,"lap%d.png",i);lap[i]->writePNG(name);sprintf(name,"sm%d.png",i);sm[i]->writePNG(name);cout << i<<": " <<  10*sc[i]/1.44<< ", ";
        //cout <<i << " "<< sc[i] << endl;
    }

    if (hesdet)sedge_lap(hes,mask,hes,hes,hes,hes,sc,features,threshold,1/*mask_radius*/,0.01/*eigenratio*/,DHESSIAN,dradius);
    else sedge_lap(hes,mask,hes,hes,hes,lap,sc,features,threshold,1/*mask_radius*/,0.01/*eigenratio*/,DHESSIAN,dradius);

    if (reg_thes!=0)regionFilter(features, reg_thes);
    desc_nb=(img->x()*img->y())/desc_nb;
    thresFilter(features, desc_nb, threshold);
    if (aff>1)findAffineRegion(sm, lap, sc, features, aff);

    cout << "Hessian-Laplace(affine) interest points "<< features.size()<< endl;


//    if((DESC&DESCRIPTOR)!=0){
//     vector<DARY *> grad;
//     vector<DARY *> ori;
//     for(uint i=0;i<sm.size();i++){
//       grad.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
//       ori.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
//       gradAngle(dx[i],dy[i],grad[i],ori[i]);
//     }
//     ncomputeDescriptors(grad, ori,  sc, features, DESC, noangle);
//     //computeDescriptors( sm, dx, dy,  sc, corners, DESC, noangle);
//     for(uint i=0;i<sm.size();i++){delete grad[i];delete ori[i];}
//     grad.clear();ori.clear();
//   }



    for (uint i=0;i<sm.size();i++) {
        delete sm[i];
        delete lap[i];
        delete hes[i];
        delete dx[i];
        delete dy[i];
        delete mask[i];
    }
    sc.clear();
    sm.clear();
    lap.clear();
    hes.clear();
    dx.clear();
    dy.clear();
    mask.clear();
}


void fastHarris(DARY *img, vector<FeatureDescriptor*>&features, Params *params) {

    float threshold=(int)params->getValue("harris_threshold.float");
    float step=params->getValue("scale_space_step.float");
    uint DESC= (uint)params->getValue("feature_type.int");
    int aff=(int)params->getValue("affine_iterations.int");
    int noangle=(int)params->getValue("descriptor_noangle.int");
    int dradius=(int)params->getValue("descriptor_radius.int");
    float reg_thes=params->getValue("region_filter.float");
    float desc_number=params->getValue("descriptor_density_per_image.int");
    fastHarris(img, features,  threshold,  step,  DESC,reg_thes, aff, noangle, dradius,desc_number);
    if ((DESC&DESCRIPTOR)!=0)
        computeDescriptors(img, features, params);
}

void fastHarris(DARY *img, vector<FeatureDescriptor*>&features, float threshold, float step, uint DESC, float reg_thes, int aff, int noangle, int dradius, int desc_nb) {

    vector<DARY *> sm;
    vector<DARY *> lap;
    vector<DARY *> har;
    vector<DARY *> har11;
    vector<DARY *> har12;
    vector<DARY *> har22;
    vector<DARY *> dx;
    vector<DARY *> dy;
    vector<float> sc;
    vector<DARY *> mask;
    buildScaleSpace(img, sm, sc,step);
    for (uint i=0;i<sm.size();i++) {
        lap.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));
        laplacian(sm[i],lap[i]);
        dx.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));
        dy.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));
        //dX2(sm[i], dx[i]);dY2(sm[i], dy[i]);
        dX2(sm[i], dx[i]);
        dY2(sm[i], dy[i]);

//      if(i>0 && i< sm.size()-1){
        har.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));
        har11.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));
        har12.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));
        har22.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));
        mask.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1, 255.0));
        harris(dx[i],dy[i],har[i],har11[i],har12[i],har22[i]);
        /*      }else {
        	har.push_back(new DARY(1,1));
        	har11.push_back(new DARY(1,1));
        	har12.push_back(new DARY(1,1));
        	har22.push_back(new DARY(1,1));
        	mask.push_back(new DARY(1,1));
              } */
        //char nm[512];sprintf(nm,"sm%d.png",i);lap[i]->writePNG(nm);
        //sprintf(nm,"har%d.png",i);har[i]->writePNG(nm);//sprintf(nm,"lap%d.png",i);lap[i]->writePNG(nm);
        //sprintf(nm,"har12%d.png",i);har12[i]->writePNG(nm);//sprintf(nm,"lap%d.png",i);lap[i]->writePNG(nm);
        //sprintf(nm,"har11%d.png",i);har11[i]->writePNG(nm);//sprintf(nm,"lap%d.png",i);lap[i]->writePNG(nm);
        //sprintf(nm,"har22%d.png",i);har22[i]->writePNG(nm);//sprintf(nm,"lap%d.png",i);lap[i]->writePNG(nm);
        //sprintf(nm,"dx%d.png",i);dx[i]->writePNG(nm);//sprintf(nm,"lap%d.png",i);lap[i]->writePNG(nm);
        //sprintf(nm,"dy%d.png",i);dy[i]->writePNG(nm);//sprintf(nm,"lap%d.png",i);lap[i]->writePNG(nm);
    }

    sedge_lap(har,mask,har11,har12,har22,lap,sc,features,threshold,1/*radius*/,0.01/*eigenratio*/,DHARRIS,dradius);
    if (reg_thes!=0)regionFilter(features, reg_thes);
    cout<< "har " << features.size()<< endl;
    desc_nb=(img->x()*img->y())/desc_nb;
    thresFilter(features, desc_nb, threshold);


    if (aff>1)findAffineRegion(sm,lap,sc,features,aff);

    cout << " Harris-Laplace(affine) interest points "<< features.size()<< endl;

//     features.push_back(new FD());
//     features[features.size()-1]->allocVec(192);
//     features[features.size()-1]->setType(DHARRIS);
//     features.push_back(new FD());
//     features[features.size()-1]->allocVec(192);
//     features[features.size()-1]->setType(DHARRIS);

//   if((DESC&DESCRIPTOR)!=0){
//     vector<DARY *> grad;
//     vector<DARY *> ori;
//     for(uint i=0;i<sm.size();i++){
//       grad.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
//       ori.push_back(new DARY(sm[i]->y(),sm[i]->x(),FLOAT1));//
//       gradAngle(dx[i],dy[i],grad[i],ori[i]);
//     }
//     ncomputeDescriptors(grad, ori,  sc, features, DESC, noangle);
//     //computeDescriptors( sm, dx, dy,  sc, corners, DESC, noangle);
//     //computeDescriptors(img, features, params);
//     for(uint i=0;i<sm.size();i++){delete grad[i];delete ori[i];}
//     grad.clear();ori.clear();
//   }


    for (uint i=0;i<sm.size();i++) {
        delete sm[i];
        delete lap[i];
        delete har[i];
        delete har11[i];
        delete har12[i];
        delete har22[i];
        delete mask[i];
        delete dx[i];
        delete dy[i];
    }
    sc.clear();
    sm.clear();
    lap.clear();
    har.clear();
    har11.clear();
    har12.clear();
    har22.clear();
    dx.clear();
    dy.clear();
    mask.clear();
}

void fastDense(DARY *img, vector<FeatureDescriptor*>&features, Params *params) {


    int dx=(int)params->getValue("dense_sampling_dx.int");
    int dy=(int)params->getValue("dense_sampling_dy.int");
    int dradius=(int)params->getValue("descriptor_radius.int");

    uint DESC= (uint)params->getValue("feature_type.int");
    //cout << dx << " " << dy << " " << dradius << endl;


    for (uint j=dradius+1;j<img->y()-dradius;j+=dy) {
        for (uint i=dradius+1;i<img->x()-dradius;i+=dx) {
            FD *fd= new FD();
            fd->setX_Y(i,j);
            fd->setScale(dradius);
            fd->setRadius(dradius);
            fd->setIntLev(0);
            features.push_back(fd);
        }
    }

    if ((DESC&DESCRIPTOR)!=0)
        computeDescriptors(img, features, params);


}

void hessian(DARY *img, vector<FeatureDescriptor*>&features, float threshold,  uint DESC, int aff, int noangle, int dradius) {

    DARY *hes=new DARY(img->y(),img->x(),FLOAT1);
    DARY *lap=new DARY(img->y(),img->x(),FLOAT1);

    hessian(img, hes, lap);

    max_detection(hes,hes,hes,hes,features,threshold,0.01/*eigenratio*/,DHESSIAN,dradius);//cout<< "haredge " << features.size()<< en

    //if(aff>1)findAffineRegion(sm,lap,sc,features,aff);
    //hes->writePNG("hes.png");
    cout << " Hessian interest points "<< features.size()<< endl;

    //if((DESC&DESCRIPTOR)!=0)computeDescriptors( sm, dx, dy, sc, features, DESC, noangle);


    delete hes;
    delete lap;

}

void hessian(DARY *img, vector<FeatureDescriptor*>&features, Params *params) {

    float threshold=(int)params->getValue("hessian_threshold.float");
    float step=params->getValue("scale_space_step.float");
    uint DESC= (uint)params->getValue("feature_type.int");
    int aff=(int)params->getValue("affine_iterations.int");
    int noangle=(int)params->getValue("descriptor_noangle.int");
    int dradius=(int)params->getValue("descriptor_radius.int");
    hessian(img, features,  threshold,  DESC, aff, noangle, dradius);
    if ((DESC&DESCRIPTOR)!=0)
        computeDescriptors(img, features, params);
}



void harris(DARY *img, vector<FeatureDescriptor*>&features, float threshold,  uint DESC, int aff, int noangle, int dradius) {

    DARY *har=new DARY(img->y(),img->x(),FLOAT1);
    DARY *har11=new DARY(img->y(),img->x(),FLOAT1);
    DARY *har12=new DARY(img->y(),img->x(),FLOAT1);
    DARY *har22=new DARY(img->y(),img->x(),FLOAT1);
    DARY *dx=new DARY(img->y(),img->x(),FLOAT1);
    DARY *dy=new DARY(img->y(),img->x(),FLOAT1);

    dX2(img, dx);
    dY2(img, dy);
    harris(dx,dy,har,har11,har12,har22);

    max_detection(har,har11,har12,har22,features,threshold,0.01/*eigenratio*/,DHARRIS,dradius);//cout<< "haredge " << features.size()<< en

    //if(aff>1)findAffineRegion(sm,lap,sc,features,aff);
    //har->writePNG("har.png");
    cout << " Harris interest points "<< features.size()<< endl;

    //if((DESC&DESCRIPTOR)!=0)computeDescriptors( sm, dx, dy, sc, features, DESC, noangle);


    delete har;
    delete har11;
    delete har12;
    delete har22;
    delete dx;
    delete dy;

}

void harris(DARY *img, vector<FeatureDescriptor*>&features, Params *params) {

    float threshold=(int)params->getValue("harris_threshold.float");
    float step=params->getValue("scale_space_step.float");
    uint DESC= (uint)params->getValue("feature_type.int");
    int aff=(int)params->getValue("affine_iterations.int");
    int noangle=(int)params->getValue("descriptor_noangle.int");
    int dradius=(int)params->getValue("descriptor_radius.int");
    harris(img, features,  threshold,  DESC, aff, noangle, dradius);
    if ((DESC&DESCRIPTOR)!=0)
        computeDescriptors(img, features, params);
}

/**********************************AFFINE**********************************************/


void getMi(DARY *img, float &a, float &b, float &c, float x, float y) {
    int row_nb= img->y();
    int col_nb= img->y();
    float  t1,t2;
    DARY  *fx= new DARY(row_nb, col_nb,FLOAT1);
    DARY  *fy= new DARY(row_nb, col_nb,FLOAT1);
    DARY  *fxy  = new DARY(row_nb, col_nb,FLOAT1);
    dX6(img, fx);
    dY6(img, fy);
    for (int row = 0; row < row_nb; row++)
        for (int col = 0; col < col_nb; col++) {
            t1 = fx->fel[row][col];
            t2 = fy->fel[row][col];
            fx->fel[row][col] = t1*t1;
            fy->fel[row][col] = t2*t2;
            fxy->fel[row][col] = t1*t2;
        }
    a= smoothf((int)x,(int)y, fx, 3);
    c= smoothf((int)x,(int)y, fy, 3);
    b= smoothf((int)x,(int)y, fxy, 3);

    delete fx;
    delete fy;
    delete fxy;
}


int fastfindAffineRegion(vector<DARY *> image,vector<DARY *> laplacian,vector<float>scale, FeatureDescriptor * cor, int lmax) {

    int level = cor->getDerLev();
    float pointx=cor->getX()/scale[level];
    float pointy=cor->getY()/scale[level];
    int sizex=19;//2*(1.44*GAUSS_CUTOFF+3)+1;
    int l,go_on=1;
    float l1=1,l2=1,ea=0;
    float eigen_ratio_act=0.1,eigen_ratio_bef=0.1,deigen=0.02;
    //static float nba=0;
    DARY *img=new DARY(sizex,sizex,FLOAT1);
    DARY *cimg=image[level];
    float a,b,c,u11=cor->getMi11(),u12=cor->getMi12(),u21=cor->getMi12(),u22=cor->getMi22(),u11t,u12t;
    getEigen(u11,u12,u21,u22,l1,l2);
    eigen_ratio_act=1-l2/l1;

    //cout << pointx << " " << pointy << " level " << level << " " << eigen_ratio_act <<endl;getchar();
    for (l=0;l<lmax && go_on;l++) { //estimate affine structure
        img->interpolate(cimg,pointx,pointy,u11,u12,u21,u22);
        //img->writePNG("img.png");
        getMi(img, a,b,c, sizex>>1, sizex>>1);
        //cout <<l1 <<"  " << l2 <<" a " << a << " b " <<b << "  c " << c << endl;
        invSqrRoot(a,b,c,l2,l1,ea);

        eigen_ratio_bef=eigen_ratio_act;
        eigen_ratio_act=1-l2/l1;

        u11t=u11;
        u12t=u12;
        u11=a*u11t+b*u21;
        u12=a*u12t+b*u22;
        u21=b*u11t+c*u21;
        u22=b*u12t+c*u22;

        //cout << u11 << " "<< u12 << endl;
        //cout << u21 << " "<< u22 << endl;
        //cout << endl << l1 << " "<< l2 << endl;
        getEigen(u11,u12,u21,u22,l1,l2);

        if (l>15 || (l1/l2>6) || (l2/l1>6)) {
            delete img;
            return 0;
        }

        if (eigen_ratio_act<deigen && eigen_ratio_bef<deigen)go_on=0;
    }
    delete img;
    //cout <<"eigen_ratio_act "<<eigen_ratio_act<<"  "<< l1<< " "<< l2 << " "<< l<< endl;getchar();

    cor->setMi(u11,u12,u21,u22,l1,l2,ea);
    return l;
}


void findAffineRegion(vector<DARY *> image,vector<DARY *> laplacian, vector<float>scale,vector<FeatureDescriptor*> &cor, int lmax) {
    for (int i=0;i<(int)cor.size();i++) {
        int l=fastfindAffineRegion(image,laplacian,scale,cor[i],lmax);
        if (l!=0) {
            //cout<<"\r  cor  "<<i<<" of "<< size << "  "<<cor[i]->getDerLev()<< "  " << cor[i]->getX() << "  " << cor[i]->getY()<<"  yes  "<< flush;
        } else {
            //cout<<"\r  cor  "<<i<<" of "<< size << "  "<<cor[i]->getDerLev()<< "  " << cor[i]->getX() << "  " << cor[i]->getY()<<"  no  "<< flush;
            cor.erase(cor.begin()+i);
            i--;
        }
    }
}


void extractFeats(DARY *image, vector<FD*> &features, Params *params) {
    unsigned long detector= (unsigned long)params->getValue("feature_type.int");
    if ((detector&DETECTOR)!=0) {
        //medianFilter(image, 2);
        if ((detector&DETECTOR)==DHARRIS) {
            fastHarris(image, features, params);
        } else if ((detector&DETECTOR)==DUNIFORM) {	    
            fastUniform(image, features, params);
        } else if ((detector&DETECTOR)==DHESSIAN) {
            fastHessian(image, features, params);
        } else if ((detector&DETECTOR)==DHES) {
            hessian(image, features, params);
        } else if ((detector&DETECTOR)==DHAR) {
            harris(image, features, params);
        } else if ((detector&DETECTOR)==DMSER) {
            fastMSER(image, features, params);
        } else if ((detector&DETECTOR)==OSIFT) {
            fastSIFT(image, features, params);
            if ((detector&DESCRIPTOR)!=0)
                computeDescriptors(image, features, params);
        } else if ((detector&DETECTOR)==DHARHES) {
            fastHarrisHessian(image, features, params);
        } else if ((detector&DETECTOR)== DSEDGE) {
            //fastSEdge(image, features,  params);
            fastDense(image, features, params);
        } else if ((detector&DETECTOR)== DSEGM) {
            //fastSegm(image, features, params);
            fastEGSF(image, features, params);
        }
    }
}


void extractFeatures(DARY *image, vector<FD*> &features, Params *params) {
    unsigned long detector= (unsigned long)params->getValue("feature_type.int");
    //cout << "extract "<< endl;
    if (image->getType()!=FLOAT1 && image->getType()!=UCHAR1FLOAT1 && image->getType()!=UCHAR3FLOAT1) {
        if (image->getType()==UCHAR1)image->convert(UCHAR3);

        if (image->getType()==UCHAR3) {
            image->RGB2opp();
            image->convert(UCHAR3FLOAT1);
        }
    }
    //cout << "OK 2 "<< image->getType()<< endl;getchar();
    //image->writePNG("test.png");
    //cout << "extract 2"<< endl;
    //medianFilter(image, 2);
    //gammaFilter(image);
    //cout << "det "<< detector << " " << (detector&DETECTOR) << " "  << getName((detector&DETECTOR))<< " " << (DSEDGE|DSIFT)<< endl;

    extractFeats(image, features, params);

    return;
    float desc_nb=params->getValue("descriptor_density_per_image.int");
    desc_nb=(image->x()*image->y())/desc_nb;
    uint prop=desc_nb/features.size();
    if (prop>2) {
        cout << " desc " <<  desc_nb << "  " << features.size() << "prop "<< prop << endl;
        //image->convert(OPP3FLOAT1);
        float sc=params->getValue("image_rescale.float");
        DARY *img=new DARY(sc*image->y(),sc*image->x(),UCHAR3FLOAT1);
        img->scale(image,1.0/sc,1.0/sc);
        //img->writePNG("im.png");
        //img->writeF("F.png");
        // img->writeR("R.png");
        //img->writeG("G.png");
        //img->writeB("B.png");
        vector<FD*> feats;
        extractFeats(img, feats, params);
        for (uint i=0;i<feats.size();i++) {
            feats[i]->setX_Y(feats[i]->getX()/sc,feats[i]->getY()/sc);
            feats[i]->setScale(feats[i]->getScale()/sc);
        }
        features.insert(features.begin(),feats.begin(),feats.end());
        feats.clear();
        float reg_thes=params->getValue("region_filter.float");
        if (reg_thes!=0)regionFilter(features, reg_thes);

        //increase image and do it again.
        delete img;
    }



}


/**********************DETECT FEATURES FOR RECOGNITION*****************************/
/**********************DETECT FEATURES FOR RECOGNITION*****************************/


void getBoundingBox(DARY *mask,int &xc, int &yc,int &width, int &height) {

    uint xs=mask->x(),ys=mask->y(),xe=0,ye=0;
    float fx=0,fy=0,sum=0;
    for (uint j=0;j<mask->y();j++) {
        for (uint i=0;i<mask->x();i++) {
            if (mask->bel[j][i]>0) {
                fx+=i;
                fy+=j;
                sum++;
                if (xs>i)xs=i;
                if (xe<i)xe=i;
                if (ys>j)ys=j;
                if (ye<j)ye=j;
            }
        }
    }
    //xc=(xs+xe)/2;
    //yc=(ys+ye)/2;
    if (xe>xs && ye>ys) {
        width=(xe-xs)/2;
        height=(ye-ys)/2;
    } else {
        width=0;
        height=0;
    }
    if (sum<10) {
        xc=0;
        yc=0;
    } else {
        xc=(int)(fx/sum);
        yc=(int)(fy/sum);
    }


}


void getObjectCenter(DARY *mask,int &x, int &y) {
    float fx=0,fy=0,sum=0;
    for (uint j=0;j<mask->y();j++) {
        for (uint i=0;i<mask->x();i++) {
            if (mask->bel[j][i]>0) {
                fx+=i;
                fy+=j;
                sum++;
            }
        }
    }
    if (sum<10) {
        x=0;
        y=0;
    } else {
        x=(int)(fx/sum);
        y=(int)(fy/sum);
    }

}

float checkOverlap(DARY *mask, FD *cor) {/*MASK is in bytes!!!*/

    int rad=(int)(0.5+cor->getRadius());
    int x=(int)cor->getX();
    int y=(int)cor->getY();
    float fig=0,score;
    //cout << x << " "<< y << " "<< cor->getScale() << " " << rad << endl;getchar();
    for (int j=y-rad;j<=y+rad;j++) {
        for (int i=x-rad;i<=x+rad;i++) {
            if (i>=0 && j>=0 && i<(int)mask->x() && j<(int)mask->y()) {
                if (mask->bel[j][i]>0)
                    fig++;
            }
        }
    }
    score=fig/square(2*rad+1);
    if (mask->bel[y][x]>0)score=1;
    cor->setArea(score);
    return score;
}

void selectWithMask(DARY *mask, vector<FD *> &desc, int &width, int &height, float overlap) {
    int x,y;
    //getObjectCenter(mask,x, y);
    getBoundingBox(mask,x, y, width, height);
    //cout << " IN " << desc.size()<< endl;
    for (uint i=0;i<desc.size();i++) {

        if (checkOverlap(mask, desc[i])>overlap) {
            desc[i]->setX_Y(x-desc[i]->getX(),y-desc[i]->getY());
            //        if((desc[i]->getType()&DCLBP)!=0){
            //   desc[i]->setX_Y(0,0);//texture descriptor has unreliable angle, can't vote for centre
            // }
        } else {
            delete desc[i];
            desc.erase(desc.begin()+i);
            i--;
        }
    }
    cout << "Selected with mask nb "<< desc.size()<< endl;
}

void extractColorDescriptor(const char *name, vector<FD *> &features,  vector<FD *> &desc, Params *params) {

    DARY *cimg = new DARY(name);
    if (cimg->getType()!=UCHAR3) {
        delete cimg;
        return;
    }

    cout << "doing color in " << name<< " for " << features.size()<< endl;

    float maxsize=0,minsize=1000;
    int nb=0;
    for (uint i=0;i<features.size();i++) {
        if ( features[i]->getType()==DHESSIAN || features[i]->getType()==DHARHES ||  features[i]->getScale()>20) {
            desc.push_back(new FD());
            desc[nb]->copy(features[i]);
            //features[i]->Cout();
            if (maxsize<desc[nb]->getScale())maxsize=desc[nb]->getScale();
            if (minsize>desc[nb]->getScale())minsize=desc[nb]->getScale();
            //computeColorDesc(cimg,desc[nb]);
            computeColorSift(cimg,desc[nb]);
            nb++;
        }
    }
    cout << "maxsize " << maxsize<< " minsize " << minsize << " nb " << nb << endl;
}

void extractMotion(const char *filein, vector<FD *> &features,  Params *params) {
    vector<FD *> desc;
    unsigned long type=(unsigned long)params->getValue("feature_type.int");
    if (( type & DMOTION)==0)return;

    cout << "Extracting motion from " << filein << "  ..."<<flush;
    char tname[512];
    char name[512];
    strcpy(tname,filein);
    char *dir=strrchr(tname,'/');
    if (!dir)dir=tname;
    else dir++;
    strcpy(name,dir);
    char *suf=strrchr(name,'.');
    sprintf(suf,"-vel.png");

    sprintf(dir,"motion/%s",name);
    //cout <<tname<< endl;
    DARY *vel=NULL;
    FILE *ident;
    if ((  ident = fopen(tname,"r") ) != NULL) {
        fclose(ident);
        vel = new DARY(tname);
    } else {
        cout<< "ERROR: no motion maps "<< tname <<endl;
        exit(0);
    }

    cout <<endl<< "Extracting motion for " << tname<<  endl;

    DARY *velx=new DARY(vel->y(),vel->x(),FLOAT1,0.0);
    DARY *vely=new DARY(vel->y(),vel->x(),FLOAT1,0.0);

    for (uint i=0;i<vely->size();i++) {
        velx->fel[0][i] = (vel->belg[0][i]-128);
        vely->fel[0][i] = (vel->belr[0][i]-128);
    }

    DARY*  grad = new DARY(vely->y(),vely->x(),FLOAT1,0.0);
    DARY*  angle = new DARY(vely->y(),vely->x(),FLOAT1,0.0);

    gradAngle(velx, vely, grad, angle);


    float ang=0;
    for (uint i=0;i<features.size();i++) {
        int xi = (int)(0.5+features[i]->getX());
        int yi = (int)(0.5+features[i]->getY());
        int size = (int)(0.5+features[i]->getScale());
        computeHistAngle(grad, angle, xi, yi, size, ang);
        features[i]->setMotion(ang);
        //cout << i << " features of " << features.size()<< " " <<ang<<  endl;
    }

    //getchar();

    delete vel;
    delete velx;
    delete vely;
    delete grad;
    delete angle;

    cout << "done " << endl;

}

void extractColorDescriptor(const char *name, vector<FD *> &features,  Params *params) {
    vector<FD *> desc;
    unsigned long type=(unsigned long)params->getValue("feature_type.int");
    if (( type & DCOLOR)==0)
        return;
    else if (( type& DESCRIPTOR)==DCOLOR) {
        extractColorDescriptor(name, features,  desc, params);
        deleteDescriptors(features);
        features.clear();
        features=desc;
    } else {
        extractColorDescriptor(name, features,  desc, params);
        features.insert(features.begin(),desc.begin(),desc.end());
        desc.clear();
    }
}

void detectFeatures(const char *name, DARY *img, DARY *mask, vector<FD *> &features,
                    int &width, int &height, Params *params) {
    vector<FD *> desc;
    int image_shift=(int)params->getValue("image_shift.int");
    float image_scale=params->getValue("image_scale.float");
    unsigned long type=(unsigned long)params->getValue("feature_type.int");
    DARY * s_img = new DARY(img->y(),img->x(),img->getType());
    s_img->set(img);
    //gray->writePNG("img.png");getchar();
    //img->flipH();
    //cout << "OK " << endl;
    if (image_shift==0 && image_scale==1) {
        extractFeatures(s_img, desc, params);
    } else {
        DARY *img2 = new DARY((int)(img->y()*image_scale),(int)(img->x()*image_scale+2*image_shift),img->getType());
        img2->scale(s_img,1.0/image_scale,1.0/image_scale);
        extractFeatures(img2, desc, params);
        //img2->writePNG("test.png");getchar();
        delete img2;
    }
    delete s_img;
    // for(uint d=0;d<desc.size();d++)desc[d]->setImageName(name);

    //cout << "tot nb "<< desc.size()<< flush;
    if (( type & DMOTION)!=0)extractMotion(name, desc, params);
    //extractColorLBP(img,gray, desc,params);
    if (( type & DCLBP)!=0)extractSegmentSift(img, desc,  params);

    if (mask!=NULL) {
        selectWithMask(mask, desc, width, height, params->getValue("train_mask_overlap.float"));
    }
    //extractColorDescriptor(name, desc,params);
//displayFeatures(img, desc, "feature.png", 255,0);cout << "features "<< endl;getchar();//remove object center
    //cout << desc.size()<< endl;
    features.insert(features.begin(),desc.begin(),desc.end());
    desc.clear();

    if (params->getValue("symmetric.int")) {
        cout << "computing Symmeric..."<<flush;
        FD *fs;
        uint nb=0;
        uint fsize=features.size();
        for (uint i=0;i<fsize;i++) {
            if ((features[i]->getType()&DSIFT)!=0) {
                fs=new FD();
                fs->copy(features[i]);
                symmetricSift(fs,img->x());
                features.push_back(fs);
                nb++;
            }
        }
        //features.clear();

        cout <<" " << nb << " total " << features.size()<< " done"<< endl;
    }

    //cout << "done"<< endl;getchar();
}


void detectFeatures(const char *filein, vector<FD *> &features, int &width, int &height,  Params *params) {

    vector<FD *> desc;
    DARY *img = new DARY(filein);
    DARY *imask=NULL;

    params->put("image_x_size.int",(float)img->x());
    params->put("image_y_size.int",(float)img->y());
    char tname[512];
    char name[512];

    strcpy(tname,filein);
    char *dir=strrchr(tname,'/');
    if (!dir)dir=tname;
    else dir++;
    strcpy(name,dir);
    char *suf=strrchr(name,'.');
    sprintf(suf,"-map.png");
    sprintf(dir,"maps/%s",name);
    //cout <<"TESTING "<< name<< " "<<  tname<< endl;
    FILE *ident;
    if ((  ident = fopen(tname,"r") ) != NULL) {
        fclose(ident);
        imask = new DARY(tname);
        imask->convert(UCHAR1);
    }
    detectFeatures(filein, img, imask, features,  width, height, params);

    delete img;
    if (imask!=NULL)delete imask;
    
}



void separateFeatures(vector<FD *> clusters, vector<FD *> &outfeat) {
    if (clusters.size()==0) {
        cout << "ERROR:problem in retrieval.separateFeatures, 0 features.\n";
    }
    vector<FD *> feat;
    FD *cl;

//  feat.push_back(new FD());feat[feat.size()-1]->setType(OSIFT | DCOLOR);
//  feat.push_back(new FD());feat[feat.size()-1]->setType((OSIFT | (OSIFT>>1)) | DCOLOR);
// feat.push_back(new FD());feat[feat.size()-1]->setType(DHARHES | DSIFT);
//  feat.push_back(new FD());feat[feat.size()-1]->setType((DHARHES | (DHARHES>>1)) | DSIFT);
    //feat.push_back(new FD());feat[feat.size()-1]->setType(DHARRIS | DSIFT);
    //feat.push_back(new FD());feat[feat.size()-1]->setType((DHARRIS | (DHARRIS>>1)) | DSIFT);
// unsigned long ntype=DHARHES;
//  feat.push_back(new FD());feat[feat.size()-1]->setType(((ntype>>2)| ntype) | DSIFT);
//  ntype=(DHARHES | (DHARHES>>1));
//  feat.push_back(new FD());feat[feat.size()-1]->setType(((ntype>>2)| ntype) | DSIFT);
//  ntype=DHARRIS;
//  feat.push_back(new FD());feat[feat.size()-1]->setType(((ntype>>2)| ntype) | DSIFT);
//  ntype=(DHARRIS | (DHARRIS>>1));
//  feat.push_back(new FD());feat[feat.size()-1]->setType(((ntype>>2)| ntype) | DSIFT);
// feat.push_back(new FD());feat[feat.size()-1]->setType(DSEDGE | DSIFT);
// feat.push_back(new FD());feat[feat.size()-1]->setType((DSEDGE | (DSEDGE>>1)) | DSIFT);
    //feat.push_back(new FD());feat[feat.size()-1]->setType(DMSER | DSIFT);
    //feat.push_back(new FD());feat[feat.size()-1]->setType((DMSER |   (DMSER>>1)) | DSIFT);
    //feat.push_back(new FD());feat[feat.size()-1]->setType(DSEGM | DSIFT);
    //feat.push_back(new FD());feat[feat.size()-1]->setType(DCLBP);
    /* feat.push_back(new FD());feat[feat.size()-1]->setType((DHARHES | (DHARHES>>1)) | DCLBP);
     feat.push_back(new FD());feat[feat.size()-1]->setType( DHARRIS| DCLBP);
     feat.push_back(new FD());feat[feat.size()-1]->setType((DHARRIS | (DHARRIS>>1)) | DCLBP);*/
//   long dhh=(DHARHES | (DHARHES>>1));
//   long dh=(DHARRIS | (DHARRIS>>1));
//   cout << ((DHARRIS | (DHARRIS>>1)) | DCLBP) << " " << ( DHARRIS| DCLBP) << " " << ((DHARHES | (DHARHES>>1)) | DCLBP) << " " << ( DHARHES| DCLBP) << " " << DCLBP << "  "<<  DHARHES << "  "<<  DHARRIS << " " <<dhh << "  "<< dh << "  " <<DSEDGE << endl;

    for (uint i=0;i<clusters.size();i++) {
        int flag=0;
        cl=clusters[i];
        //if(cl->getType()==134218240 || cl->getType()==134218496)
        //{cout << i << " " << cl->getType() << "  " << feat[0]->getType()<< "  " << feat[1]->getType()<< endl;}
        for (uint t=0;t<feat.size() && !flag;t++) {//t+=2){
            if (cl->getType()==feat[t]->getType()) {
                flag=1;
                //if(cl->getExtr()<0)
                feat[t]->features.push_back(cl);
                //else feat[t+1]->features.push_back(cl);
            }
        }
        if (!flag) {
            cout <<"Descriptors: separateFeatures no type "<< cl->getType() << "  "<< (cl->getType()&ADETECTOR) << "  " << (cl->getType()& DESCRIPTOR) << endl;//exit(1);
            feat.push_back(new FD());
            feat[feat.size()-1]->setType(cl->getType());
            feat[feat.size()-1]->features.push_back(cl);
        }
    }
    vector<pair<uint,uint> > type;

    for (uint i=0;i<feat.size();i++) {
        type.push_back(make_pair(feat[i]->getType(),i));
        //cout <<feat[i]->getType() << " " <<  feat[i]->features.size()<< endl;
    }
    sort(type.begin(),type.end());
    //cout << "type nb "<<  type.size()<< endl;
    for (uint i=0;i<feat.size();i++) {
        //cout <<"type "<<   getName(feat[type[i].second]->getType()&ADETECTOR)<< endl;
        if (feat[type[i].second]->features.size()>0) {
            outfeat.push_back(feat[type[i].second]);
        }
    }
    feat.clear();
    type.clear();

    cout << "Number of separated types: "<< outfeat.size()<< endl;


    /*
    for(uint t=0;t<feat.size();t++){

      getchar();
      for(uint i=0;i<feat[t]->features.size();i++){
        if(feat[t]->getType()!=feat[t]->features[i]->getType())cout << "ERROR "<<feat[t]->getType()<<  endl;
        cout << feat[t]->features.size()<< " " << feat[t]->getType()<< " " << feat[t]->features[i]->getType()<<" " <<feat[t]->features[i]->getExtr() <<endl;
      }
    }
    */

}
