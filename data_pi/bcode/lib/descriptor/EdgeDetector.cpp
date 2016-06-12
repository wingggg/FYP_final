#include "feature.h"
#include "../gauss_iir/gauss_iir.h"



void cannyEdges(DARY *img, DARY *edge,  float scale, float lower_threshold, float higher_threshold){
  
  DARY *dx = new DARY(img->y(),img->x(),FLOAT1);
  DARY *dy = new DARY(img->y(),img->x(),FLOAT1);
  DARY *grad = new DARY(img->y(),img->x(),FLOAT1);

  dX(img,dx,scale);
  dY(img,dy,scale);
  for(uint j=0;j<grad->y();j++){
    for(uint i=0;i<grad->x();i++){
      grad->fel[j][i]=sqrt(dx->fel[j][i]*dx->fel[j][i]+dy->fel[j][i]*dy->fel[j][i]); 	
    }
  } 
  cannyEdges(dx, dy, grad, edge, lower_threshold, higher_threshold); 

  delete dx;delete dy;delete grad;
}

void cannyEdges(DARY *dx,DARY *dy, DARY *edge,  float lower_threshold, float higher_threshold){
  
  DARY *grad = new DARY(dx->y(),dx->x(),FLOAT1);

  for(uint j=0;j<grad->y();j++){
    for(uint i=0;i<grad->x();i++){
      grad->fel[j][i]=sqrt(dx->fel[j][i]*dx->fel[j][i]+dy->fel[j][i]*dy->fel[j][i]); 	
    }
  } 
  cannyEdges(dx, dy, grad, edge, lower_threshold, higher_threshold); 

  delete grad;
}



void cannyEdges(DARY *img, DARY *edge,  float lower_threshold, float higher_threshold){
  
  DARY *dx = new DARY(img->y(),img->x(),FLOAT1);
  DARY *dy = new DARY(img->y(),img->x(),FLOAT1);
  DARY *grad = new DARY(img->y(),img->x(),FLOAT1);

  dX2(img,dx);
  dY2(img,dy);
  for(uint j=0;j<grad->y();j++){
    for(uint i=0;i<grad->x();i++){
      grad->fel[j][i]=sqrt(dx->fel[j][i]*dx->fel[j][i]+dy->fel[j][i]*dy->fel[j][i]); 	
    }
  } 
  cannyEdges(dx, dy, grad, edge, lower_threshold, higher_threshold); 

  delete dx;delete dy;delete grad;
}


void updateLabels(vector<int> &labels, vector<int> &lvalid, int &label, int oldlabel){

  int odl=labels[oldlabel];
  int ndl=labels[label];

  if(odl<ndl){
    ndl=odl;
    odl=labels[label];
  }
  label=ndl;

  if(lvalid[ndl]==1 || lvalid[odl]==1){
    for(uint i=1;i<labels.size();i++){
      if(labels[i]==odl || labels[i]==ndl){
	labels[i]=ndl;
	lvalid[i]=1;
      }
    }    
  }else {
    for(uint i=1;i<labels.size();i++){
      if(labels[i]==odl){
	labels[i]=ndl;
      }
    }        
  }
}


void cannyEdges(DARY *dx, DARY *dy, DARY *grad, DARY *edge, DARY *tmp_edge, float lower_threshold,
		float higher_threshold, int &edge_nb){

  vector<int> valid_edge;valid_edge.push_back(0);
  vector<int> lab_edge;lab_edge.push_back(0);
  int color=1,maxcol=255;;
  float x1,y1,x2,y2,ux,uy,g,g1,g2;
  for(uint j=1;j<grad->y()-2;j++){
    for(uint i=1;i<grad->x()-2;i++){
      g=grad->fel[j][i];
      if(g<=lower_threshold)continue;
      ux=dx->fel[j][i];
      uy=dy->fel[j][i];
      x1=i+ux/g;
      y1=j+uy/g;
      x2=i-ux/g;
      y2=j-uy/g;
      g1=grad->getValue(x1,y1);
      g2=grad->getValue(x2,y2);  
      if(g<g1 || g<=g2)continue;
      if(g>higher_threshold){
	edge->fel[j][i]=maxcol;
      }
      if(g>lower_threshold){
	if(tmp_edge->fel[j-1][i-1]){
	  tmp_edge->fel[j][i]=tmp_edge->fel[j-1][i-1];	  
	}else if(tmp_edge->fel[j-1][i]){
	  tmp_edge->fel[j][i]=tmp_edge->fel[j-1][i];
	}else if(tmp_edge->fel[j-1][i+1]){
	  tmp_edge->fel[j][i]=tmp_edge->fel[j-1][i+1];
	}else if(tmp_edge->fel[j][i-1]){
	  tmp_edge->fel[j][i]=tmp_edge->fel[j][i-1];
	}else if(g>higher_threshold){
	  tmp_edge->fel[j][i]=color;valid_edge.push_back(1);lab_edge.push_back(color);color++;
	}else {
	  tmp_edge->fel[j][i]=color;valid_edge.push_back(2);lab_edge.push_back(color);color++;
	}
	if(g>higher_threshold)valid_edge[(int)tmp_edge->fel[j][i]]=1;
      }      
    }
  }

  //for(uint i=1;i<lab_edge.size();i++){if(lab_edge[i]!=i)cout << "error"<< endl;}
  int lab;
  for(uint j=1;j<tmp_edge->y()-2;j++){
    for(uint i=1;i<tmp_edge->x()-2;i++){
      if(tmp_edge->fel[j][i]!=0){
	lab=lab_edge[(int)tmp_edge->fel[j][i]];
	if(tmp_edge->fel[j][i+1]!=0 && lab!=lab_edge[(int)tmp_edge->fel[j][i+1]])
	  updateLabels(lab_edge,valid_edge, lab, (int)tmp_edge->fel[j][i+1]);	
	
	if(tmp_edge->fel[j+1][i]!=0 && lab!=lab_edge[(int)tmp_edge->fel[j+1][i]])
	  updateLabels(lab_edge, valid_edge, lab, lab_edge[(int)tmp_edge->fel[j+1][i]]);
	
	if(tmp_edge->fel[j+1][i+1]!=0 && lab!=lab_edge[(int)tmp_edge->fel[j+1][i+1]])
	  updateLabels(lab_edge, valid_edge, lab, lab_edge[(int)tmp_edge->fel[j+1][i+1]]);

	if(tmp_edge->fel[j+1][i-1]!=0 && lab!=lab_edge[(int)tmp_edge->fel[j+1][i-1]])
	  updateLabels(lab_edge, valid_edge, lab, lab_edge[(int)tmp_edge->fel[j+1][i-1]]);
	
      }
    }
  }


  for(uint j=1;j<tmp_edge->y()-2;j++){
    for(uint i=1;i<tmp_edge->x()-2;i++){
      if(tmp_edge->fel[j][i]>0){
	if(valid_edge[lab_edge[(int)tmp_edge->fel[j][i]]]==1 || edge->fel[j][i])edge->fel[j][i]=maxcol;
	else edge->fel[j][i]=1;
	tmp_edge->fel[j][i]=lab_edge[(int)tmp_edge->fel[j][i]];
	if(edge_nb<tmp_edge->fel[j][i])edge_nb=(int)tmp_edge->fel[j][i];
      }
    }
  }  
  edge_nb++;
  
  for(uint j=0;j<grad->y();j++){
    for(uint i=0;i<grad->x();i++){
      //if(edge->fel[j][i]>0)edge->fel[j][i]=grad->fel[j][i];
      //else edge->fel[j][i]=0;
    }
  }
  
  lab_edge.clear();
  valid_edge.clear();
  
}

void histAngle(Segment *segs, float min_edge){
  for (int i = 0; i < SegOriBins; i++)
    segs->hist[i] = 0.0;
 
  float gval,fbin,dbin;
  float max_angle=M_2PI;//must be M_2PI - for dominant orientation
  int bin;
  for(uint i=0;i<segs->x.size();i++){
    gval = segs->grad[i];  
    while(segs->angle[i]<0)segs->angle[i]+=max_angle;
    while(segs->angle[i]>max_angle)segs->angle[i]-=max_angle;
   
    fbin =  (SegOriBins * segs->angle[i] / max_angle);
    bin=(int)floor(fbin);
    dbin=fbin-bin;
    //assert(bin >= 0 && bin <= OriBins);
    //cout << bin << " " << SegOriBins << " " << segs->hist[bin] << endl;
    bin = (bin < SegOriBins)?bin:(0);
    segs->hist[bin] +=  (1-dbin)*gval;
    bin = (bin+1 < SegOriBins)?bin+1:(0);
    segs->hist[bin] +=  (dbin) * gval;    
  }
  int maxi=-1; 
  float maxval=0, interp=0;
  for (int i = 0; i < SegOriBins; i++){
    if (segs->hist[i] > maxval){
      maxval = segs->hist[i];
      maxi=i;
    }
  }

  int prev,next;
  prev = (maxi == 0 ? SegOriBins - 1 : maxi - 1);
  next = (maxi == SegOriBins - 1 ? 0 : maxi + 1);
  interp = interpPeak(segs->hist[prev], segs->hist[maxi], segs->hist[next]);

  segs->ang_mean=(max_angle * (maxi + interp) / SegOriBins);
  while(segs->ang_mean<0)segs->ang_mean+=max_angle;
  while(segs->ang_mean>M_PI)segs->ang_mean-=max_angle;
  float meanx=0,meany=0;
  for(uint i=0;i<segs->x.size();i++){
    meanx+=segs->x[i];
    meany+=segs->y[i];
  }
  segs->x_mean=meanx/segs->x.size();
  segs->y_mean=meany/segs->x.size();
  float x2=0,y2=0,xy=0;
  for(uint i=0;i<segs->x.size();i++){
    segs->x[i]-=segs->x_mean;
    segs->y[i]-=segs->y_mean;
    x2+=(segs->x[i])*(segs->x[i]);
    y2+=(segs->y[i])*(segs->y[i]);
    xy+=(segs->y[i])*(segs->x[i]);   
  }  
  x2/=segs->x.size();
  y2/=segs->x.size();
  xy/=segs->x.size();
  float l2=0,l1=0,ea=0;
  //cout << "OK "<< endl;
  getEigen(x2, xy, xy, y2, l1, l2);
  //cout << " x2 " << sqrt(x2) << " y2 " << sqrt(y2) << " l1 "  <<  sqrt(l1) << " l2 " << sqrt(l2) << " ea " << ea << " mx " <<segs->x_mean << " my " << segs->y_mean  <<" ang " <<  360*segs->ang_mean/max_angle<< endl;
  segs->scale=2*sqrt(l1);
  //getchar();
  
  return;
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //cout << endl;
  
  vector<float> angles;
  for (int i = 0; i < SegOriBins; i++){
    prev = (i == 0 ? SegOriBins - 1 : i - 1);
    next = (i == SegOriBins - 1 ? 0 : i + 1);
	  //cout << i  << " "<< hist[i];
    //if (segs->hist[i] > segs->hist[prev] && segs->hist[i] > segs->hist[next] && segs->hist[i]>min_edge){
    if (segs->hist[i]>min_edge){
      //interp = interpPeak(segs->hist[prev], segs->hist[i], segs->hist[next]);
      //angles.push_back(M_PI * (i + 0.5 + interp) / OriBins); 
      angles.push_back(M_PI * (i + interp) / SegOriBins); 
      //cout <<i <<" " <<  angles[angles.size()-1];      
    } 
    cout << segs->hist[i]<< " ";
  }   
 cout << endl;
  float minangle,atmp;
  int aind=0;
  cout << "size " <<segs->x.size()<<  endl;
  for(uint i=0;i<angles.size();i++){
    segs->segs.push_back(new Segment());
    cout << i << " angle " << angles[i] << endl;

  }

  for(uint a=0;a<segs->angle.size();a++){
    minangle=M_PI;
    for(uint i=0;i<angles.size();i++){
      atmp=fabs(segs->angle[a]-angles[i]);
      if(atmp<minangle){
	minangle=atmp;
	aind=i;
      }      
      segs->segs[aind]->x.push_back(segs->x[a]);;
      segs->segs[aind]->y.push_back(segs->y[a]);;
      segs->segs[aind]->grad.push_back(segs->grad[a]);;
      segs->segs[aind]->angle.push_back(segs->angle[a]);;
    } 
  }
  
  /* Apply smoothing 6 times for accurate Gaussian approximation. */
  //for(int i = 0; i < 6; i++)smoothHistogram(segs->hist, SegOriBins)


}





void cannyEdges(DARY *dx, DARY *dy, DARY *grad, DARY *edge,  float lower_threshold,
		float higher_threshold){
  DARY *lab_edge = new DARY(dx->y(),dx->x(),FLOAT1,0.0);
  int edge_nb=0;
  cannyEdges(dx, dy, grad, edge, lab_edge, lower_threshold, higher_threshold, edge_nb);
  delete lab_edge;
}

int breakPoint(DARY *lab_edge ,DARY *har, int i, int j){
  float h=0,h0=0,h1=0,h2=0,h3=0,h4=0,h5=0,h6=0,h7=0;
  int nb=0;
  if(lab_edge->fel[j-1][i-1]){h0=har->fel[j-1][i-1];nb++;}
  if(lab_edge->fel[j-1][i]){h1=har->fel[j-1][i];nb++;}
  if(lab_edge->fel[j-1][i+1]){h2=har->fel[j-1][i+1];nb++;}
  if(lab_edge->fel[j][i-1]){h3=har->fel[j][i-1];nb++;}
  if(lab_edge->fel[j][i])h=har->fel[j][i];
  if(lab_edge->fel[j][i+1]){h4=har->fel[j][i+1];nb++;}
  if(lab_edge->fel[j+1][i-1]){h5=har->fel[j+1][i-1];nb++;}
  if(lab_edge->fel[j+1][i]){h6=har->fel[j+1][i];nb++;}
  if(lab_edge->fel[j+1][i+1]){h7=har->fel[j+1][i+1];nb++;}


  if(h>100 && nb>2 && h>h0 && h>h1 && h>h2 && h>h3 && h>h4 && h>h5 && h>h6 && h>h7)return 1;
  else return 0;    
}


void getAngles(DARY *dx, DARY *dy,  DARY *ori, DARY *edge, float ang_mean, float &anglemin,  float &fadiffmin,  int x, int y, int &xn, int &yn, float angle_thres){

  if(edge->fel[y][x]==0)return;
  
  float angle=ori->fel[y][x];
  float fadiff=fabs(ang_mean-angle);
  //cout << angle << " fadiff " << fadiff<< endl;
  if(fadiff>M_PI){
    fadiff=M_2PI-fadiff;
    if(fadiff<angle_thres){
      if(ang_mean<0)
	angle=ang_mean-fadiff;
      else angle=ang_mean+fadiff;    
    }
  }  
  if(fadiff<fadiffmin){
    fadiffmin=fadiff;
    anglemin=angle;
    xn=x;
    yn=y;
  }

}

void followContour(DARY *dx, DARY *dy, DARY *edge, DARY *grad, DARY *ori, int x, int y, float angle, Segment *seg){   


  float minangle, fadiffmin=1000, angle_thres=0.5;
  int minx,miny;  
  

  edge->fel[y][x]=0;
  seg->x.push_back(x);
  seg->y.push_back(y);
  seg->angle.push_back(angle);
  seg->ang_mean=((seg->angle.size()-1)*seg->ang_mean + angle)/seg->angle.size();
  seg->grad.push_back(grad->fel[y][x]);
  

  //edge->writePNG("edge.png");cout << "OK "<< endl;
  //find point with min angle difference 
  getAngles(dx,dy,ori,edge, seg->ang_mean, minangle,  fadiffmin, x, y-1, minx, miny, angle_thres);
  getAngles(dx,dy,ori,edge, seg->ang_mean, minangle,  fadiffmin, x-1, y, minx, miny, angle_thres);
  getAngles(dx,dy,ori,edge, seg->ang_mean, minangle,  fadiffmin, x+1, y, minx, miny, angle_thres);
  getAngles(dx,dy,ori,edge, seg->ang_mean, minangle,  fadiffmin, x, y+1, minx, miny, angle_thres);
  
  //cout << x << " y "  << y << " a "  << seg->ang_mean << " x "  << minx << " y "  << miny << " a "  << minangle << " fad " << fadiffmin << " thres " << angle_thres << endl;getchar();
  //follow the contour if gradient angle is similar
  if(fadiffmin<angle_thres){
    followContour(dx, dy, edge, grad, ori, minx, miny, minangle, seg);
  }else {//otherwise do the same neighbouring pixels in corners
    getAngles(dx,dy,ori,edge, seg->ang_mean, minangle,  fadiffmin, x-1, y-1, minx, miny, angle_thres);
    getAngles(dx,dy,ori,edge, seg->ang_mean, minangle,  fadiffmin, x+1, y-1, minx, miny, angle_thres);
    getAngles(dx,dy,ori,edge, seg->ang_mean, minangle,  fadiffmin, x-1, y+1, minx, miny, angle_thres);
    getAngles(dx,dy,ori,edge, seg->ang_mean, minangle,  fadiffmin, x+1, y+1, minx, miny, angle_thres);
    //cout << x << " y "  << y << " a "  << seg->ang_mean << " x "  << minx << " y "  << miny << " a "  << minangle << " fad " << fadiffmin << " thres " << angle_thres << endl;getchar();
    if(fadiffmin<angle_thres)
      followContour(dx, dy, edge, grad, ori, minx, miny, minangle, seg);
  }
}


void drawEdge(DARY *edge, Segment *seg, int color){
  
  for(uint j=0;j<seg->x.size();j++){
    //cout << (seg->y_mean+seg->y[j]) << " " << (int)(seg->x_mean+seg->x[j])<< endl;
    edge->fel[(int)(seg->y_mean+seg->y[j])][(int)(seg->x_mean+seg->x[j])]=color;
  }  
}

void merge(Segment *nseg,Segment *seg2, int how){
  //cout << " how "<< how << " " << nseg->x.size() << " " << seg2->x.size() << endl;
  if(how==1){
    for(uint k=0;k<seg2->x.size();k++){      
      nseg->x.insert(nseg->x.begin(),seg2->x[k]);
      nseg->y.insert(nseg->y.begin(),seg2->y[k]);
      nseg->angle.insert(nseg->angle.begin(),seg2->angle[k]);
      nseg->grad.insert(nseg->grad.begin(),seg2->grad[k]);
    }
  }else if(how==2){
    nseg->x.insert(nseg->x.begin(),seg2->x.begin(),seg2->x.end());
    nseg->y.insert(nseg->y.begin(),seg2->y.begin(),seg2->y.end());
    nseg->angle.insert(nseg->angle.begin(),seg2->angle.begin(),seg2->angle.end());
    nseg->grad.insert(nseg->grad.begin(),seg2->grad.begin(),seg2->grad.end());      
  }else if(how==3){
    nseg->x.insert(nseg->x.end(),seg2->x.begin(),seg2->x.end());
    nseg->y.insert(nseg->y.end(),seg2->y.begin(),seg2->y.end());
    nseg->angle.insert(nseg->angle.end(),seg2->angle.begin(),seg2->angle.end());
    nseg->grad.insert(nseg->grad.end(),seg2->grad.begin(),seg2->grad.end());      
  }else if(how==4){
    for(int k=seg2->x.size()-1;k>=0;k--){
      nseg->x.push_back(seg2->x[k]);
      nseg->y.push_back(seg2->y[k]);
      nseg->angle.push_back(seg2->angle[k]);
      nseg->grad.push_back(seg2->grad[k]);
    }
  }
  
}
  

void merging(vector<Segment *> &segs, uint minsize){

  for(uint j=0;j<segs.size();j++){
    if(segs[j]!=NULL && segs[j]->x.size()<minsize && segs[j]->x.size()>0){
      int xs=segs[j]->x[0];
      int ys=segs[j]->y[0];
      int xe=segs[j]->x[segs[j]->x.size()-1];
      int ye=segs[j]->y[segs[j]->y.size()-1];
	
      int flag=1;
      for(uint k=0;k<segs.size() && flag;k++){
	if(segs[k]!=NULL && k!=j)if(segs[k]->x.size()>0){
	  //cout << k << " of " << segs.size() << " " << xs << " " << ys << " e " << xe << " " << ye << " " << j << " jk " << k << " " << segs[j]->x.size()<< endl;
	  int dx=abs(xs-segs[k]->x[0]);
	  int dy=abs(ys-segs[k]->y[0]);
	  if(dx<2 && dy<2){
	    merge(segs[j],segs[k],1);
	    flag=0;
	  }
	  dx=abs(xs-segs[k]->x[segs[k]->x.size()-1]);
	  dy=abs(ys-segs[k]->y[segs[k]->y.size()-1]);
	  if(dx<2 && dy<2 && flag){
	    merge(segs[j],segs[k], 2);
	    flag=0;
	  }
	  dx=abs(xe-segs[k]->x[0]);
	  dy=abs(ye-segs[k]->y[0]);
	  if(dx<2 && dy<2 && flag){
	    merge(segs[j],segs[k], 3);
	    flag=0;
	  }
	  dx=abs(xe-segs[k]->x[segs[k]->x.size()-1]);
	  dy=abs(ye-segs[k]->y[segs[k]->y.size()-1]);
	  if(dx<2 && dy<2 && flag){
	    merge(segs[j],segs[k], 4);
	    flag=0;
	  }
	  if(!flag){	    
	    delete segs[k];segs[k]=NULL;
	    j--;
	  }
	}	  
      }      
    }
  } 
  for(uint j=0;j<segs.size();j++){
    if(segs[j]==NULL){
      segs.erase(segs.begin()+j); 
      j--;
    }else if(segs[j]->x.size()<minsize){
      delete segs[j];
      segs.erase(segs.begin()+j); 
      j--;      
    }    
  } 
  if(segs.size()<2)return;
  Segment *nseg = new Segment();
  for(uint j=0;j<segs.size();j++){
    nseg->x.insert(nseg->x.end(),segs[j]->x.begin(),segs[j]->x.end());
    nseg->y.insert(nseg->y.end(),segs[j]->y.begin(),segs[j]->y.end());
    nseg->angle.insert(nseg->angle.end(),segs[j]->angle.begin(),segs[j]->angle.end());
    nseg->grad.insert(nseg->grad.end(),segs[j]->grad.begin(),segs[j]->grad.end());      
  }
  segs.push_back(nseg);  
}

void cannySegments(DARY *dx, DARY *dy, DARY *grad, DARY *ori, DARY *cedge,  vector<FD*> &features,  float lower_threshold, float higher_threshold){
  DARY *lab_edge = new DARY(dx->y(),dx->x(),FLOAT1,0.0);
  int edge_nb=0;
  vector<Segment*> segs;
  //Timer times;
  cannyEdges(dx, dy, grad, cedge, lab_edge, lower_threshold, higher_threshold, edge_nb);

  DARY *edge = new DARY(cedge);
  
  for(int j=0;j<edge_nb;j++){
    segs.push_back(NULL);
  }
  int lab=0;
  float angle=0;
  for(uint j=1;j<lab_edge->y()-2;j++){
    for(uint i=1;i<lab_edge->x()-2;i++){           
      if(edge->fel[j][i]>0){
	lab=(int)lab_edge->fel[j][i];
	if(segs[lab]==NULL){
	  segs[lab]=new Segment();
	  if(edge->fel[j][i]>10)
	  segs[lab]->valid=1;
	  segs[lab]->segs.push_back(new Segment());
	}else segs[lab]->segs.push_back(new Segment());
	angle=ori->fel[j][i];
	//if(angle<0)angle+=M_PI;
	followContour(dx, dy, edge, grad, ori, i, j, angle, segs[lab]->segs[segs[lab]->segs.size()-1]);
      }
    }
  }  
  int  nb=0;
  for(uint j=0;j<segs.size();j++){
    if(segs[j]==NULL){
      segs.erase(segs.begin()+j);
      j--;
    }else if(segs[j]->valid==0){
      delete segs[j];
      segs.erase(segs.begin()+j);
       j--;     
    }else {
      nb+=segs[j]->segs.size();
    }
    /*else if(segs[j]->x.size()<10){
    segs.erase(segs.begin()+j);
      j--;
      }else {
      //histAngle(segs[j],5,edge);
      }*/
  } 
  //cout << " nb1 "<< nb << endl;

  nb=0;
  int min_edge=20;
  for(uint j=0;j<segs.size();j++){
    //cout << j << " " <<segs.size() << "  "<<  segs[j]->segs.size();
    if(segs[j]->segs.size()>0)merging(segs[j]->segs, min_edge);     
    if(segs[j]->segs.size()==0){
      delete segs[j];
      segs.erase(segs.begin()+j);
      j--;
    }
  }
  uint epoints=0;
  uint segnb=0;
  uint sxegnb=0;
  for(uint j=0;j<segs.size();j++){
    segnb+=segs[j]->segs.size();
    epoints+=segs[j]->segs[nb]->x.size();
      for(uint p=5;p<segs[j]->segs[nb]->x.size()-5;p++){
        
      }
   
  }
  cout<<"segnb "<< segnb << " " << sxegnb<< " epoints " << epoints << " " << epoints/segnb << endl;
  //times.stop();
  cedge->set(0.0);
  vector<pair<float,uint > > mangle;
 
  for(uint j=0;j<segs.size();j++){
    uint nb=segs[j]->segs.size()-1;
    uint last=segs[j]->segs[nb]->x.size()-1;
    //cout << " seg " << segs[j]->segs.size()<<" " << last <<  " " << segs[j]->segs[nb]->x[0]<< endl;
    for(uint p=0;p<last;p++){
      cedge->fel[segs[j]->segs[nb]->y[p]][segs[j]->segs[nb]->x[p]]=100;
     
    }
    int ln=5;
    cedge->fel[segs[j]->segs[nb]->y[0]][segs[j]->segs[nb]->x[0]]=255;
    cedge->fel[segs[j]->segs[nb]->y[last]][segs[j]->segs[nb]->x[last]]=255;
    mangle.clear();
      for(uint p=ln;p<last-ln;p++){
      float ma1=0,ma2=0,a1,a2;
      for(uint p1=0;p1<=ln;p1++){      
        a1=segs[j]->segs[nb]->angle[p-p1];
        if(a1<-M_PI)a1+=M_2PI;
        else if(a1>M_PI)a1-=M_2PI;
        a2=segs[j]->segs[nb]->angle[p+p1];
        if(a2<-M_PI)a2+=M_2PI;
        else if(a2>M_PI)a2-=M_2PI;
        ma1+=a1;
        ma2+=a2;
      }
      mangle.push_back(make_pair(fabs(ma1/ln-ma2/ln),p));
     }
     std::sort(mangle.begin(),mangle.end());
     int la=mangle.size()-1;
     cedge->fel[segs[j]->segs[nb]->y[mangle[la].second]][segs[j]->segs[nb]->x[mangle[la].second]]=255;
/*     for(int p=la-1;p>0;p--){
       if(abs(la-p)>15 && mangle[la].first>0.3)
         cedge->fel[segs[j]->segs[nb]->y[mangle[la].second]][segs[j]->segs[nb]->x[mangle[la].second]]=255;
       cout << p << " " << mangle[p].first << " " << mangle[p].second<< endl; 
     }*/
     cedge->writePNG("segm.png");cout << "segm " << endl;getchar();
  }
  FD *fd;
  for(uint j=0;j<segs.size();j++){
     //cout << j << " of " << segs.size()<<" size  " <<segs[j]->segs.size()  <<  endl;
     nb+=segs[j]->segs.size();
     for(uint i=0;i<segs[j]->segs.size();i++){
       //histAngle(segs[j]->segs[i], min_edge);
       //drawEdge(edge, segs[j]->segs[i],255);
       //edge->writePNG("edge.png");getchar();
       //fd = new FD();       
       //segs[j]->segs[i]->ang_mean=0;
       //computeSift(segs[j]->segs[i],fd);
       //features.push_back(fd);
        //cout <<"     "<<  i << " of " << segs[j]->segs.size() << " sizepnb " << segs[j]->segs[i]->x.size()<< endl;
     }     
     delete segs[j];
   } 
   segs.clear();
   //cout << "OK 3 "<< segs.size() << " " << nb << endl;
   //edge->writePNG("edge.png");
   //getchar();
   
  delete edge;
  delete lab_edge;

}




void gammaFilter(DARY *img){
  for(uint j=0;j<img->size();j++){
    img->fel[0][j]=16*sqrt(img->fel[0][j]);
  }
}


void medianFilter(DARY *img, int size){

  int dim2=(2*size+1);
  int dim=dim2*dim2;
  int v=0;
  float *fel;
  DARY *img_copy=new DARY(img);
  float *array = new float[dim];
  for(uint j=size;j<img->y()-size-1;j++){
    for(uint i=size;i<img->x()-size-1;i++){
      v=0;
      for(int n=-size;n<=size;n++){
	fel=img_copy->fel[j+n]+i;
	for(int m=-size;m<=size;m++){

	  array[v++]=(float)fel[m]; 
	  //cout << (int)array[v-1]<< " " ;
	}
      } 
      //cout << endl;
      img->fel[j][i]=median(array, dim);
    }
  }  
  delete img_copy;
  delete []array;
}


