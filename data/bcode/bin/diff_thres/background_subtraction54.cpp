#include "../../descriptor/feature.h"
#include "../../util/util.h"
#include <time.h>



struct BB
{
    int    w;
    int    h;
    int    x;
    int    y;
    int    sx;
    int    sy;
    int    ex;
    int    ey;
    int valid;
    float confidence,sparsity,coverage;
};


inline int getRoot(int *conn, int l){
    while(conn[l]!=l)l=conn[l];        
    return l;
    // cout << "             " << l << " = " << conn[l]<< endl;
}

# define INT_MAX 100000

void labelSegments(DARY *output, int &labels_nb){
    int lab=50000;  
    int *conn = new int[lab];
    bzero(conn, sizeof(int)*lab);
    for(int j=1; j<lab; j++){ conn[j]=j; }
    lab=255;
    int p0,p1,p2,p3,p4,pmin;
    for(uint j=0; j<output->y(); j++){
        output->fel[j][0]=0;
        output->fel[j][output->x()-1]=0;
    }
    for(uint j = 0; j<output->x(); j++){
        output->fel[0][j]=0;
        output->fel[output->y()-1][j]=0;
    }

    for(uint j=1; j<output->y()-1; j++){
        for(uint i=1; i<output->x()-1; i++){
            if(output->fel[j][i]>0){
                p0=output->fel[j][i];
                p1=INT_MAX; p2=INT_MAX; p3=INT_MAX; p4=INT_MAX;

                // cout << "test0" << endl;

                if(output->fel[j][i-1] > 0){
                    p1=output->fel[j][i-1]; 
                    //cout << " test1 " << p1;
                    p1=getRoot(conn,p1);
                    //cout << p1;
                }
                if(output->fel[j-1][i-1]>0){
                    p2=output->fel[j-1][i-1];    
                    //cout << " test2 " << p2;
                    p2=getRoot(conn,p2);
                    //cout << p2;
                } 
                if(output->fel[j-1][i]>0){
                    p3=output->fel[j-1][i];
                    //cout << " test3 " << p3;
                    p3=getRoot(conn,p3);
                    //cout << p3;
                } 
                if(output->fel[j-1][i+1]>0){
                    p4=output->fel[j-1][i+1];
                    //cout << " test4 " << p4;
                    p4=getRoot(conn,p4);
                    //cout << " test4 " << p4;
                }

                //cout << endl;

                pmin=(p1<p2)?p1:p2; pmin=(pmin<p3)?pmin:p3; pmin=(pmin<p4)?pmin:p4;
                if(pmin==INT_MAX){
                    output->fel[j][i]=lab;
                    lab++;
                } else {
                    output->fel[j][i]=pmin;
                    if(p1<INT_MAX){ conn[p1]=pmin; }
                    if(p2<INT_MAX){ conn[p2]=pmin; }
                    if(p3<INT_MAX){ conn[p3]=pmin; }
                    if(p4<INT_MAX){ conn[p4]=pmin; }
                }
            }
        }
    }

    int *unique = new int[lab];
    bzero(unique, sizeof(int)*lab);
    labels_nb=1;
    for(int j=255; j<lab; j++){
        conn[j]=getRoot(conn,j);
        if(conn[j]==j){
            unique[j]=labels_nb;
            labels_nb++;
        }
    }

    // int *uniquect = new int[labels_nb];
    // bzero(uniquect,sizeof(int)*labels_nb);

    for(uint j=1; j<output->y()-1; j++){
        for(uint i=1; i<output->x()-1; i++){
            p0=output->fel[j][i];
            if(p0>0){
                output->fel[j][i]=unique[conn[p0]];  
                // uniquect[unique[conn[p0]]]++;
            }
        }
    }

    //for(int j=0; j<labels_nb; j++)cout << j << " count "<< uniquect[j]<< endl;

    delete []conn;
    delete []unique;
    // output->writePNG("labels.png");
    // cout << "saved " << lab << endl;//getchar();
}


void coutBB(BB bb){
    cout << "BB v: "<< bb.valid << " cov: " << bb.coverage<< " conf: " << bb.confidence<< " sx " << bb.sx << " ex "<< bb.ex << " sy " << " " << bb.sy << " ey " << " " << bb.ey;
    cout  << " x " << " " << bb.x << " y " << " " << bb.y <<  " w " << " " << bb.w << " h " << " " << bb.h <<  endl;
}


void findBBs(DARY *binary, vector<BB> &bbs, int labels){
     
    float *mx = new float[labels];
    float *my = new float[labels];
    float *h = new float[labels];
    float *w = new float[labels];
    float *count = new float[labels];

    float cover_weights[11] = {1.0195, 1.2459, 1.1371, 1.1124, 0.9766, 0.8936, 0.8690, 0.8717, 0.8899, 0.9506, 1.0000};

    bzero(mx,labels*sizeof(float));
    bzero(my,labels*sizeof(float));
    bzero(h,labels*sizeof(float));
    bzero(w,labels*sizeof(float));
    bzero(count,labels*sizeof(float));

    int lab;
    for(uint j=0; j<binary->y(); j++){
        for(uint i=0; i<binary->x(); i++){
            if(binary->fel[j][i] > 0){
                lab=binary->fel[j][i];
                mx[lab]+=i;
                my[lab]+=j;
                count[lab]++;	     
            }
        }
    }
    for(int i=1;i<labels;i++){
        mx[0]+=mx[i]; 
        my[0]+=my[i]; 
        count[0]+=count[i]; 
    }

    for(int i=0;i<labels ;i++){
        if(count[i]==0)count[i]=1;
        mx[i]/=count[i];
        my[i]/=count[i];
    }

    for(uint j = 0 ;j<binary->y();j++){
        for(uint i = 0 ;i<binary->x();i++){
            if(binary->fel[j][i] > 0){
                lab=binary->fel[j][i];
                if(lab>=labels){ cout << lab << " " << labels << " i " << i << " j " << j <<  endl; }
                h[lab]+=(j-my[lab])*(j-my[lab]);
                w[lab]+=(i-mx[lab])*(i-mx[lab]);
                h[0]+=(j-my[0])*(j-my[0]);
                w[0]+=(i-mx[0])*(i-mx[0]);
            }
        }
    }

    float scale=1.5;
    int cind;
    for(int lab=0;lab<labels;lab++){
        //cout << lab << " " << h[lab] << " h c " << count[lab] << endl;
        h[lab]/=count[lab];
        w[lab]/=count[lab];
        h[lab]=scale*sqrt(h[lab]);
        w[lab]=scale*sqrt(w[lab]);
        BB bb;
        bb.valid=1;
        bb.x=mx[lab];
        bb.y=my[lab];
        bb.w=w[lab];
        bb.h=h[lab];
        bb.x=(((bb.x<0)?0:bb.x)>=binary->x())?(binary->x()-1):bb.x;
        bb.y=(((bb.y<0)?0:bb.y)>=binary->y())?(binary->y()-1):bb.y;
        bb.sx=((bb.x-bb.w)<0)?0:(bb.x-bb.w);
        bb.sy=((bb.y-bb.h)<0)?0:(bb.y-bb.h);
        bb.ex=((bb.x+bb.w)>=binary->x())?(binary->x()-1):(bb.x+bb.w);
        bb.ey=((bb.y+bb.h)>=binary->y())?(binary->y()-1):(bb.y+bb.h);
        bb.coverage=((bb.ex-bb.sx)*(bb.ey-bb.sy))/(float)(binary->x()*binary->y());
        cind=(int)(10*bb.coverage);
        bb.confidence = cover_weights[cind]*count[lab]/((bb.ex-bb.sx)*(bb.ey-bb.sy));       
        bbs.push_back(bb);
    }
    float norm=1.0;
    if(bbs.size()){ norm=sqrt(sqrt(norm/bbs.size())); }
    for(uint lab=0;lab<bbs.size();lab++){
        bbs[lab].confidence=bbs[lab].confidence*norm;
    }

    //cout <<  "labels "<< labels << endl;
    //getchar();

    delete []mx;
    delete []my;
    delete []w;
    delete []h;
    delete []count;

}


void detectDifferences(vector<DARY *> images, DARY *output, int thres, int median_s){
    DARY *binary=new DARY(output->y(), output->x(), UCHAR1FLOAT1);
    output->set(0.0);
    binary->set(0.0);

    for(uint j=0; j<2; j++){
        for(uint i=0; i<images[0]->size(); i++){
            output->fel[0][i] = output->fel[0][i] + fabs(images[j]->belr[0][i]-128) + fabs(images[j]->belg[0][i]-128) + fabs(images[j]->belb[0][i]-128);
        }
    }

    //output->convert(UCHAR1FLOAT1);
    /*    int hsize=256;
    int *hist=new int[hsize];

    bzero(hist,sizeof(int)*hsize);
    for(int i = 0 ;i<output->size();i++)hist[output->bel[0][i]]+=1;

    for(int i=0;i<hsize;i++)cout << hist[i]<< " ";
    cout << endl;
    */

    for(uint i=0; i<binary->size(); i++){
        output->fel[0][i] = (output->fel[0][i]>thres)?255:0;
        // output->fel[0][i]=binary->fel[0][i];
    }
    // output->writePNG("diff.png");cout << "median saved " << endl;//getchar();

    // binary->median2d(output,median_s);
    output->erosion(binary,1,255);
    output->set(0.0);
    binary->dilation(output,2,255);

    // binary->writePNG("bin.png");
    delete binary;       
}


#define G4    0.0110
#define G3    0.0432
#define G2    0.1146
#define G1    0.2060
#define G0    0.2504

void getHoG(DARY *output, DARY *gray, vector<DARY *> hog){ 
    float dx,dy,grad,ori,dori;
    int iori;
    DARY *tmp=new DARY(output->y(), output->x(), UCHAR1FLOAT1);
    for(uint j=1;j<output->y()-1;j++){
        for(uint i=1;i<output->x()-1;i++){
            if(output->fel[j][i-1]>0 && output->fel[j][i+1]>0 && output->fel[j-1][i]>0 && output->fel[j+1][i]>0){
                dx=gray->fel[j][i+1]-gray->fel[j][i-1];	
                dy=gray->fel[j-1][i]-gray->fel[j+1][i];
                grad=sqrt(dx*dx+dy*dy);
                ori=3.99*(1.0+atan2(dy,dx)/M_PI);
                iori=(int)floor(ori);
                dori=ori-iori;
                hog[iori]->fel[j][i]=(1.0-dori)*grad;
                hog[((iori+1)<8?(ori+1):0)]->fel[j][i]=(dori)*grad;
                //	if(ori ==8)
                //  cout  << ori << endl;
            }
        }
    }
    for(uint h=0;h<hog.size();h++){
        tmp->set(0.0);
        for(uint j=4;j<output->y()-4;j++){
            for(uint i=4;i<output->x()-4;i++){
                if(output->fel[j][i]>0){
                    tmp->fel[j][i]=G4*(hog[h]->fel[j][i-4]+hog[h]->fel[j][i+4])+G3*(hog[h]->fel[j][i-3]+hog[h]->fel[j][i+3])+
                    G2*(hog[h]->fel[j][i-2]+hog[h]->fel[j][i+2])+G1*(hog[h]->fel[j][i-1]+hog[h]->fel[j][i+1])+G0*hog[h]->fel[j][i];
                }
            }
        }
        hog[h]->set(0.0);
        for(uint j=4;j<output->y()-4;j++){
            for(uint i=4;i<output->x()-4;i++){
                if(output->fel[j][i]>0){
                    hog[h]->fel[j][i]=G4*(tmp->fel[j-4][i]+tmp->fel[j+4][i])+G3*(tmp->fel[j-3][i]+tmp->fel[j+3][i])+
                    G2*(tmp->fel[j-2][i]+tmp->fel[j+2][i])+G1*(tmp->fel[j-1][i]+tmp->fel[j+1][i])+G0*tmp->fel[j][i];
                }
            }
        }   
    }      
    delete tmp;
    float norm;
    for(uint j=4;j<output->y()-4;j++){
        for(uint i=4;i<output->x()-4;i++){
            norm=0;
            for(int h=0;h<8;h++){
                norm+=square(hog[h]->fel[j][i]);
            }

            norm=sqrt(norm);
            if(norm>0){
                for(int h=0;h<8;h++){
                    hog[h]->fel[j][i]/=norm;
                }
            }
        }
    }
}

void hogDistance(DARY *output, vector<DARY *> hog1, vector<DARY *> hog2, DARY *map){
  float dist;
  for(uint j=4;j<output->y()-4;j++){
	for(uint i=4;i<output->x()-4;i++){
	  dist=0;
	  if(output->fel[j][i]>0){
 	    for(int h=0;h<8;h++){	    
 	      dist+=square(hog1[h]->fel[j][i]-hog2[h]->fel[j][i]);
 	    }
 
 	    map->fel[j][i]=sqrt(dist);
	  }
	}
   }
}


void matchImages(vector<DARY *> &images, DARY *output){
    DARY *gray=new DARY(output->y(), output->x(), UCHAR1FLOAT1);
    vector<DARY *> hog1;
    vector<DARY *> hog2; 
    vector<DARY *> hog3;

    // create new HOGs
    for(uint i=0;i<8;i++){
        hog1.push_back(new DARY(output->y(), output->x(), UCHAR1FLOAT1));
        hog1[i]->set(0.0);
        hog2.push_back(new DARY(output->y(), output->x(), UCHAR1FLOAT1));
        hog2[i]->set(0.0);
        hog3.push_back(new DARY(output->y(), output->x(), UCHAR1FLOAT1));
        hog3[i]->set(0.0);
    }
    // create HOG for each of three images
    for(int k=2;k<5;k++){
        gray->set(0.0);
        for(uint j=0;j<output->size();j++){
            if(output->fel[0][j]>0){
                // cout << output->fel[0][j] << endl;
                gray->fel[0][j] = images[k]->belr[0][j] + images[k]->belg[0][j] + images[k]->belb[0][j];
                // cout << gray->fel[0][j] << endl;
            }
        }
        if(k==2){ getHoG(output, gray, hog1); }
        if(k==3){ getHoG(output, gray, hog2); }
        if(k==4){ getHoG(output, gray, hog3); }
    }
    // for(uint i=0;i<hog1.size();i++){
    //     cout << "hog images " << endl;char name[512]; sprintf(name,"hog%d.png",i);hog1[i]->normalize();hog1[i]->writePNG(name);
    // }

    delete gray;
    gray->set(0.0);
    DARY *map12=new DARY(output->y(), output->x(), UCHAR1FLOAT1);
    DARY *map13=new DARY(output->y(), output->x(), UCHAR1FLOAT1);
    DARY *map23=new DARY(output->y(), output->x(), UCHAR1FLOAT1);
    map12->set(0.0);
    map13->set(0.0);
    map23->set(0.0);
    hogDistance(output,hog1,hog2, map12);
    hogDistance(output,hog1,hog3, map13);
    hogDistance(output,hog2,hog3, map23);

    float min=100, max=0, dmean=0, tot=0;
    for(uint j=4; j<output->y()-4; j++){
        for(uint i=4;i<output->x()-4;i++){
            if(output->fel[j][i] > 0){
                map12->fel[j][i] += (map13->fel[j][i]+map23->fel[j][i]);
                dmean += map12->fel[j][i];
                tot++;
                //cout << i << " " << j << " " << dmean << " m12 " << map12->fel[j][i] << " m13 " << map13->fel[j][i] << " m23 " << map23->fel[j][i]<< endl;getchar();
                if(max < map12->fel[j][i]){ max = map12->fel[j][i]; }
                if(min > map12->fel[j][i]){ min = map12->fel[j][i]; }
                if(map12->fel[j][i] < 2.0){ output->fel[j][i] = 0; }
                else{ output->fel[j][i] = 255; }
            }
        }
    }
    dmean/=tot;
    //cout << "map images " << max << " " << min << " " << dmean << " " << tot << endl;char name[512];map12->writePNG("map12.png");getchar();

    delete map12;
    delete map13;
    delete map23;
    for(uint i=0;i<8;i++){
        delete hog1[i];
        delete hog2[i];
        delete hog3[i];
    }
    hog1.clear();
    hog2.clear();
    hog3.clear();
}


void normalizeDifference(uchar **bel, int x_size, int y_size){
    int cellnb = 5;
    int allcellnb = cellnb*cellnb;
    int xsize = 1 + x_size/cellnb;
    int ysize = 1 + y_size/cellnb;

    float *means= new float[allcellnb];
    float *count= new float[allcellnb];
    bzero(means, allcellnb*sizeof(float));
    bzero(count, allcellnb*sizeof(float));
    int x=0, y=0, b;

    for(int j=0; j<y_size; j++){
        y = (int)(j/ysize);
        for(int i=0; i<x_size; i++){
            x = (int)(i/xsize);
            b = cellnb*y+x;
            //if(b>24 || b<0)cout << "ERROR !!!!!!!!!!!!!!!!!!!" << i << " " << x << " " << j << " " << ysize << " " << y << endl;
            means[b] += bel[j][i];
            count[b]++;
        }
    }

    for(int i=0; i<allcellnb; i++){
        means[i] /= count[i];
        // cout << means[i] << " ";
    }

    float meddiff = median(means, allcellnb);
    meddiff = (int)(meddiff - 128);    
    int p;

    if(abs(meddiff)){
        //       cout << "meddiff " << meddiff  << endl;
        for(int j=0; j<y_size; j++){
            for(int i=0; i<x_size; i++){
                p = bel[j][i] - meddiff;
                bel[j][i] = (((p>255)?255:p)<0)?0:p;
            }
        }
    }
    delete []means;
    delete []count;
}


void normalizeDifference(vector<DARY *> &images){

    normalizeDifference(images[0]->belr, images[0]->x(), images[0]->y());
    normalizeDifference(images[0]->belg, images[0]->x(), images[0]->y());
    normalizeDifference(images[0]->belb, images[0]->x(), images[0]->y());

    normalizeDifference(images[1]->belr, images[0]->x(), images[0]->y());
    normalizeDifference(images[1]->belg, images[0]->x(), images[0]->y());
    normalizeDifference(images[1]->belb, images[0]->x(), images[0]->y());

}

// Calculates the RGB differences between frames 2,3 and 2,4. Results are returned in 0,1.
void computeDifferences(vector<DARY *> &images){
    int p;
    for(uint i = 0; i<images[0]->size(); i++){  //images[0]->size()==1958400.
        p = images[2]->belr[0][i] - images[3]->belr[0][i] + 128;
        images[0]->belr[0][i] = (((p>255)?255:p)<0)?0:p;
        p = images[2]->belg[0][i] - images[3]->belg[0][i] + 128;
        images[0]->belg[0][i] = (((p>255)?255:p)<0)?0:p;
        p = images[2]->belb[0][i] - images[3]->belb[0][i] + 128;
        images[0]->belb[0][i] = (((p>255)?255:p)<0)?0:p;

        p = images[2]->belr[0][i] - images[4]->belr[0][i] + 128;
        images[1]->belr[0][i] = (((p>255)?255:p)<0)?0:p;
        p = images[2]->belg[0][i] - images[4]->belg[0][i] + 128;
        images[1]->belg[0][i] = (((p>255)?255:p)<0)?0:p;
        p=images[2]->belb[0][i] - images[4]->belb[0][i] + 128;
        images[1]->belb[0][i] = (((p>255)?255:p)<0)?0:p;
    }

    // images[0]->writePNG("d1.png");
    //  images[1]->writePNG("d2.png");
     
}
 
void loadImages(vector<char *> filenames, int start, vector<DARY *> &images, Params *par){
    // cout << filenames[start]<< endl;
    DARY *sim,*im;
    int width = par->getValue("image_width.int");  
    int height = par->getValue("image_height.int");
    int mtop = par->getValue("header_size.int");

    if(start==-1){ // returns two images
        DARY *imin = new DARY(filenames[2]);
        int mhight = imin->y() - par->getValue("footnote_size.int") - mtop;
        
        im = new DARY(mhight, imin->x(), UCHAR3);    
        im->crop(imin, 0, mtop);

        // float scalex = width/(float)im->x(); 
        // height = (int)(im->y()/scalex);
        // sim = new DARY(height,width,UCHAR3);
        // sim->scale(im, scalex, scalex);
        // delete im;

        images.push_back(im);

        imin = new DARY(filenames[4]);
        im = new DARY(mhight, imin->x(), UCHAR3);    
        im->crop(imin, 0, mtop);

        // scalex=width/(float)im->x(); 
        // height = (int)(im->y()/scalex);
        // sim=new DARY(height,width,UCHAR3);
        // sim->scale(im, scalex, scalex);
        // delete im;
        
        images.push_back(im);
        delete imin;
    } else {
        DARY *imin =new DARY(filenames[start]);
        int mhight = imin->y() - par->getValue("footnote_size.int") - mtop;
        //cout << mhight << endl;
        im = new DARY(mhight,imin->x(),UCHAR3);
        im->crop(imin,0,mtop);

        float scalex=(float)im->x()/width;    
        height = (int)(im->y()/scalex);
        sim=new DARY(height,width,UCHAR3);
        sim->decrease(im);
        //sim->writePNG("test1.png");
        delete imin;

        //TEST
        // images.push_back(sim);        
        images.push_back(new DARY(sim));        //images[0]
        images.push_back(new DARY(sim));        //images[1]
        images.push_back(sim);  //              //images[2]
        //TESTEND

        imin=new DARY(filenames[start+1]);
        im->crop(imin,0,mtop);
        sim=new DARY(height,width,UCHAR3);
        sim->decrease(im);
        //sim->writePNG("test2.png");
        delete imin;
        images.push_back(sim);                  //images[3]

        imin=new DARY(filenames[start+2]);
        im->crop(imin,0,mtop);
        sim=new DARY(height,width,UCHAR3);
        sim->decrease(im);
        //sim->writePNG("test3.png");getchar();
        delete imin;
        delete im;
        images.push_back(sim);                  //images[4]     // returns five images
    }
}


void removeTimestamps(DARY * output, int hs, int he){ 
  for(int y = hs ;y<he;y++){
   for(uint x = 0 ;x<output->x();x++){
    output->fel[y][x]=0;
   }
  }
}

void  getSparsity(DARY *output, float &mom){
   float mx=0,my=0,tot=0;
  for(uint j=1;j<output->y()-1;j++){
     for(uint i=1;i<output->x()-1;i++){
       if(output->fel[j][i]){
	 mx+=i;
	 my+=j;
	 tot++;
       }              
     }               
  }
  mx/=tot;
  my/=tot;
   for(uint j=1;j<output->y()-1;j++){
     for(uint i=1;i<output->x()-1;i++){
       if(output->fel[j][i]){
	 mom+=fabs((i-mx)*(j-my));
       }              
     }                
  }
  mom/=tot;
  mom/=(output->x()*output->y());
 // cout << "mom " << mom << " ntot " << ntot << endl;
    
}


void drawBB(char *outname, DARY *image, vector<BB> bbs, int selected){
    char outputfile[512];
    sprintf(outputfile,"%s.bin.png",outname);

    // cout << selected << endl;

    if(selected>=0){   
        for(uint i = 0 ;i<image->size();i++){
            if(image->fel[0][i]>0){
                if((int)image->fel[0][i]==selected && selected>0){
                    image->fel[0][i]=255;
                }else{
                    image->fel[0][i]=170; 
                }
            }
        }
    }

    // //    for(int b=0;b<bbs.size();b++){
    if(selected>=0){
        BB bb = bbs[selected];

        // cout << outputfile << " " << bbs.size()<< " sel " << selected<< " " ;
        // coutBB(bbs[selected]);
    
        for(int i = bb.sx ;i<bb.ex;i++){ image->fel[bb.sy][i]=128; }
        for(int i = bb.sx ;i<bb.ex;i++){ image->fel[bb.ey][i]=128; }

        for(int j = bb.sy ;j<bb.ey;j++){ image->fel[j][bb.sx]=128; }
        for(int j = bb.sy ;j<bb.ey;j++){ image->fel[j][bb.ex]=128; }

    }
    // // }
    image->writePNG(outputfile);
}


int selectBBs(vector<BB> &bbs, Params *par){
  int largest=0;
  float minCov = par->getValue("min_coverage.float");
  float maxCov = par->getValue("max_coverage.float");
  float conf = par->getValue("confidence_thres.float");
  if(bbs[0].coverage>maxCov || bbs[0].coverage<minCov || bbs[0].confidence < conf){
    bbs[0].valid=0;
    //largest=-1;
  }
  
  for(uint i=1;i<bbs.size();i++){   //find largest bbs
	 if(bbs[0].valid==0 ||  bbs[i].coverage<minCov || bbs[i].coverage>maxCov || bbs[i].confidence < conf){
	   bbs[i].valid=0;
	 }else{
	   if(bbs[0].valid==1 && bbs[i].coverage>bbs[largest].coverage)
	    largest=i;
	    //coutBB(bbs[i]);
	 }
   //  cout << i << " ";coutBB(bbs[i]);
 }
 
 if(largest>=0){
    if(bbs[largest].coverage<minCov || bbs[largest].confidence < conf)bbs[0].valid=0;
 }
 
 return largest;
  
}


int processSet(Params *par, vector<char *> filenames, int findex, vector<BB> &bbs){
    // cout << "TEST\n";

    vector<DARY *> images;
    // cout << " OK 0 " << findex <<  " "<< filenames[findex]<< endl;
    loadImages(filenames, findex, images, par);
    // cout << " OK 1 " << endl;
    if(findex>0){ computeDifferences(images); }
    // cout << " OK 2 " << endl;
    normalizeDifference(images);

    DARY *output = new DARY(images[0]->y(),images[0]->x(),UCHAR1FLOAT1);
    output->set(0.0);
    // cout << " OK 3 " << endl;
    detectDifferences(images, output, par->getValue("diff_thres.int"), par->getValue("median_size.int"));

    // Detect using HOG (histogram of oriented gradients)
    if((int)par->getValue("match_images.int")){ matchImages(images,output); }
    for(uint i=0; i<images.size(); i++){ delete images[i];images.clear(); }
    // cout << " OK 4 " << endl;

    // Bounding box method
    int labels;
    labelSegments(output, labels);
    int sel=-1;
    if(labels>1){
        findBBs(output,bbs,labels);
        sel=selectBBs(bbs, par);
        // if(bbs.size()>0){ sel=0; }
    }
    // cout << " OK 4 " << endl;
    if(findex==-1){ findex=1; }
    // coutBB used in drawBB prints to terminal
    if((int)par->getValue("draw_images.int")){ drawBB(filenames[findex],output,bbs,sel); }

    delete output;
    return sel;
}

int main(int argc, char **argv){
    Params *par = new Params();
    // par->put("input_type.char","diffs");
    par->put("input_type.char","fulls");
    par->put("header_size.int",30.0);
    par->put("footnote_size.int",30.0);  
    par->put("camera_type.char","reconyx"); // can change camera type
    par->put("diff_thres.int",54.0);
    par->put("median_size.int",9.0);
    par->put("output_image.char","out.png");
    par->put("input_images.char","images.txt");
    par->put("output_file.char","output.txt");
    par->put("draw_images.int",1.0);
    par->put("match_images.int",0.0);
    // par->put("match_images.int",1.0);
    par->put("image_width.int",340);
    par->put("image_height.int",280);
    par->put("min_coverage.float",0.007);
    par->put("max_coverage.float",0.80);
    par->put("confidence_thres.float",0.24);
    par->put("conf_factor.float",10);

    for(int i=1;i<argc;i++){ 
        if(!strcmp(argv[i],"-p")){
            par->load(argv[i+1]);
        }else if(!strcmp(argv[i],"-i")){
            par->put("input_images.char",argv[i+1]);
        }else if(!strcmp(argv[i],"-o")){
            par->put("output_file.char",argv[i+1]);
        }
    }
    vector<char *> filenames;
    vector<BB> bbs;
    loadFileNames(par->getString("input_images.char"),filenames);
    if(filenames.size()<1){ return 1; }

    float conf_factor = par->getValue("conf_factor.float");
    ofstream textoutput(par->getString("output_file.char"));

    // cout << "TEST01" << endl;
    // cout << par->getString("input_type.char") << endl;
    // cout << filenames[10] << endl;

    if(!strcmp(par->getString("input_type.char"),"fulls")){
        for(uint i=0; i<filenames.size(); i+=3){
            int sel = processSet(par, filenames, i, bbs);

            if(bbs.size()>0){
                textoutput <<  filenames[i] << " " << bbs[0].valid << " " << bbs[sel].sx << " " << bbs[sel].sy << " " <<  bbs[sel].ex-bbs[sel].sx << " " << bbs[sel].ey-bbs[sel].sy  << " " << bbs[sel].confidence   << " " << bbs[sel].coverage  <<  "\n";
                // textoutput <<  filenames[i] << " " << bbs.size() << " " << bbs[sel].ex << " " << bbs[sel].ey << " " << bbs[sel].w << " " << bbs[sel].h << " " << bbs[sel].x << " " << bbs[sel].y << "\n";
                
            }else{
                textoutput <<  filenames[i] << " " << 0 << " " << 0 << " " << 0 << " " <<  0 << " " << 0  << " " << 0   << " " << 0  << "\n";
            }
            cout << i << " " << filenames[i] << " of " << filenames.size() <<" bbs.size " <<  bbs.size()<< endl;

            bbs.clear();
        }
    }else{
        // (input_type.char == "diffs") --> computeDifferences(images)
        int sel = processSet(par, filenames, -1, bbs);
        //textoutput << 0 << endl;
        if(bbs.size()>0){
            if(bbs[0].valid){
                //currently outputting this to output file
                textoutput <<  filenames[1] << " " << 1 << " " << round(conf_factor*bbs[sel].confidence) << " " << bbs[sel].sx << " " << bbs[sel].sy << " " <<  bbs[sel].ex-bbs[sel].sx << " " << bbs[sel].ey-bbs[sel].sy <<endl;
            }
        }
    }

    // cout << "TESTOUT" << endl;


    textoutput.close();

    //par->save("params.par");
    delete par;
    return 0;
}

