#include "../ImageContent/imageContent.h"
#include "../gauss_iir/gauss_iir.h"
#include "../descriptor/feature.h"
//#include <cstdlib>
//#include <algorithm>
#include "segment.hh"



universe::universe(int elements) {
  elts = new uni_elt[elements];
  num = elements;
  for (int i = 0; i < elements; i++) {
    elts[i].rank = 0;
    elts[i].size = 1;
    elts[i].p = i;
  }
}
  
universe::~universe() {
  delete [] elts;
}

int universe::find(int x) {
  int y = x;
  while (y != elts[y].p)
    y = elts[y].p;
  elts[x].p = y;
  return y;
}

void universe::join(int x, int y, float &w) {
  if (elts[x].rank > elts[y].rank) {
    elts[y].p = x;
    elts[x].size += elts[y].size;
    elts[x].w+=w;
    w=elts[x].w;
  } else {
    elts[x].p = y;
    elts[y].size += elts[x].size;
    elts[y].w+=w;
    w=elts[y].w;
    if (elts[x].rank == elts[y].rank)
      elts[y].rank++;
  }
  num--;
}


/*
 * Segment a graph
 *
 * Returns a disjoint-set forest representing the segmentation.
 *
 * num_vertices: number of vertices in graph.
 * num_edges: number of edges in graph
 * edges: array of edges.
 * c: constant for treshold function.
 */ 
 
 
    inline bool operator==(const rgb &a, const rgb &b) {
      return ((a.r == b.r) && (a.g == b.g) && (a.b == b.b));
    }

    bool operator<(const edge &a, const edge &b) {
      return a.w < b.w;
    } 
    bool operator>(const edge &a, const edge &b) {
      return a.w > b.w;
    } 
    bool operator<=(const edge &a, const edge &b) {
      return a.w <= b.w;
    } 
    bool operator==(const edge &a, const edge &b) {
      return a.w == b.w;
    } 
    
    bool operator<(const segment &a, const segment &b) {
      return a.w < b.w;
    } 
    bool operator>(const segment &a, const segment &b) {
      return a.w > b.w;
    } 
    bool operator<=(const segment &a, const segment &b) {
      return a.w <= b.w;
    } 
    bool operator==(const segment &a, const segment &b) {
      return a.w == b.w;
    } 

universe *segment_graph(int num_vertices, int num_edges, edge *edges, float c) { 
  // sort edges by weight
  std::sort(edges, edges + num_edges);

  // make a disjoint-set forest
  universe *u = new universe(num_vertices);

  // init thresholds
  float *threshold = new float[num_vertices];
  for (int i = 0; i < num_vertices; i++)
    threshold[i] = THRESHOLD(1,c);


  // for each edge, in non-decreasing weight order...
  for (int i = 0; i < num_edges; i++) {
    edge *pedge = &edges[i];
    float w =pedge->w;
    // components conected by this edge
    int a = u->find(pedge->a);
    int b = u->find(pedge->b);
    if (a != b) {
      if ((pedge->w <= threshold[a]) &&//was before &&
           (pedge->w <= threshold[b])) {
        u->join(a, b, w);
        a = u->find(a);
        //threshold[a] = pedge->w  +0.5;//+THRESHOLD(u->size(a), c);
        threshold[a] = pedge->w + THRESHOLD(u->size(a), c);
      }
    } 
  }

  // free up
  delete threshold;
  return u;
}



// random color
rgb random_rgb(){ 
  rgb c;
  c.r = (uchar)random();
  c.g = (uchar)random();
  c.b = (uchar)random();

  return c;
}

 
//sigma - To smooth input image before segmenting it. Recommended value 0.5.
//k - Value for threshold function. Recommended value 400.
//min - Minimum component size. Recommended value, 15.
     
void graph_segmentation(DARY *image, int k, float grad_thres, float bound_thres, int min_size, int &num_ccs, DARY *output){

  int width=image->x();
  int height=image->y(); 
  DARY *r= new DARY(height,width,FLOAT1,0.0);
  DARY *g= new DARY(height,width,FLOAT1,0.0);
  DARY *b= new DARY(height,width,FLOAT1,0.0);
  DARY *sr= new DARY(height,width,FLOAT1,0.0);
  DARY *sg= new DARY(height,width,FLOAT1,0.0);
  DARY *sb= new DARY(height,width,FLOAT1,0.0);
  DARY *igrad= new DARY(height,width,FLOAT1,0.0);
  DARY *iori= new DARY(height,width,FLOAT1,0.0);
  DARY *labout= new DARY(height,width,FLOAT1,0.0);
  image->get(r,g,b);//this is taking care of gray or color
  medianFilter(r, 2);
  medianFilter(g, 2);
  medianFilter(b, 2);
  smooth5(r,sr);
  smooth5(g,sg);
  smooth5(b,sb);
  
  
  vector<segment *> epoints;

  segment_image(sr, sg, sb,   k, grad_thres, bound_thres, min_size, num_ccs, igrad, iori, labout, epoints);
  //medianFilter(labout, 2);
  //openingFilter(labout, 2);
  //openingFilter(labout, 2);
  //openingFilter(labout, 1);
  //medianFilter(labout, 1);
  //cout << "OKseg" << endl;
  if(output->getType()==FLOAT1){
    output->set(labout);
        
  }else if(output->getType()==UCHAR3){
    rgb *colors = new rgb[width*height];
    float *rave = new float[width*height];
    float *gave = new float[width*height];
    float *bave = new float[width*height];
    float *nbave = new float[width*height];
    for (uint i = 0; i < (uint)(width*height); i++){
      colors[i] = random_rgb();
      nbave[i]=0;
      rave[i]=0;
      gave[i]=0;
      bave[i]=0;
    }
    uint comp;
    for (int y = 0; y < height-1; y++) {
      for (int x = 0; x < width-1; x++) {
        comp=(uint)labout->fel[y][x];
        nbave[comp]++;
        rave[comp]+=sr->fel[y][x];
        gave[comp]+=sg->fel[y][x];
        bave[comp]+=sb->fel[y][x];

      }
    }
    for (uint y = 0; y < (uint)(height*width); y++) {
      if(nbave[y]){
        rave[y]/=nbave[y];
        gave[y]/=nbave[y];
        bave[y]/=nbave[y];
      //cout << nbave[y] << " " << rave[y]<< " " << gave[y]<< " " << bave[y]<< endl;
      }
    
    }  
    uint lab,lab1;
    uint mk=0x00FF;

    uint maxlab=0;
    for (int y = 0; y < height-1; y++) {
      for (int x = 0; x < width-1; x++) {
        lab=(uint)labout->fel[y][x];
        if(maxlab<lab)maxlab=lab;
      //random color
/*        output->belr[y][x] = colors[lab].r;
        output->belg[y][x] = colors[lab].g;
        output->belb[y][x] = colors[lab].b;*/
      //average segment color
      output->belr[y][x] = rave[lab];
        output->belg[y][x] = gave[lab];
        output->belb[y][x] = bave[lab];
      
      //label
/*      lab1=mk&lab;
        output->belr[y][x]=lab1;
        lab1=mk&(lab>>8);
        output->belg[y][x]=lab1;
        lab1=mk&(lab>>16);
        output->belb[y][x]=lab1;*/
      

      }
    }
//     lab=maxlab;
//     lab1=mk&lab;
//     output->belr[0][0]=lab1;
//     lab1=mk&(lab>>8);
//     output->belg[0][0]=lab1;
//     lab1=mk&(lab>>16);
//     output->belb[0][0]=lab1;

  
    delete []colors;
    delete []nbave;
    delete []rave;
    delete []gave;
    delete []bave;
  }
  //image->get(ImageContent *r, ImageContent *g, ImageContent *b)
 //igrad->writePNG("igrad.png");iori->normalize(0,M_PI/2);iori->writePNG("iori.png");getchar();
 //labout->writePNG("labout.png");
  delete r;
  delete g;
  delete b;
  delete sr;
  delete sg;
  delete sb;
  delete igrad;
  delete iori;delete labout;
  
  
}

void graph_segmentation(DARY *image, DARY *output, Params *params){
  int k=(int)params->getValue("segment_k.int");
  float grad_thres=(float)params->getValue("segment_grad_thres.float");
  float bound_thres=(float)params->getValue("segment_bound_thres.float");
  int min_size=(int)params->getValue("segment_min_size.int");
  int num_ccs;
  graph_segmentation(image,  k,  grad_thres,  bound_thres,  min_size, num_ccs, output);
}


static inline float lsquare(float x) { return x*x; }

// dissimilarity measure between pixels
static inline float grad(DARY *r, DARY *g, DARY *b, int x1, int y1, int x2, int y2) {
  return sqrt(lsquare(r->fel[y1][x1]-r->fel[y2][x2]) +
      lsquare(g->fel[y1][x1]-g->fel[y2][x2]) +
      lsquare(b->fel[y1][x1]-b->fel[y2][x2]));
}
static inline float grad2(DARY *r, DARY *g, DARY *b, int pa, int pb) {
  return sqrt(lsquare(r->fel[0][pa]-r->fel[0][pb]) +
      lsquare(g->fel[0][pa]-g->fel[0][pb]) +
      lsquare(b->fel[0][pa]-b->fel[0][pb]));
}


void openingFilter(DARY *segments, DARY *sizes, int size){

  float *fel,*lel;
  DARY *img_copy=new DARY(segments);
  DARY *sizes_copy=new DARY(sizes);
  for(uint j=size;j<img_copy->y()-size-1;j++){
    for(uint i=size;i<img_copy->x()-size-1;i++){
      float maxsize=0;
      float maxlab=0;
      for(int n=-size;n<=size;n++){
        fel=sizes->fel[j+n]+i;
        lel=segments->fel[j+n]+i;
        for(int m=-size;m<=size;m++){
          if(maxsize<fel[m]){
            maxsize=fel[m];
            maxlab=lel[m];
            
          }
        }
      } 
      img_copy->fel[j][i]=maxlab;
      sizes_copy->fel[j][i]=maxsize;
/*      for(int n=-size;n<=size;n++){
        fel=segments->fel[j+n]+i;
        for(int m=-size;m<=size;m++){
         fel[m]=maxlab;
	  //cout << (int)array[v-1]<< " " ;
        }
      } */
                   
    }
  }  
  
  for(uint j=size;j<img_copy->y()-size-1;j++){
    for(uint i=size;i<img_copy->x()-size-1;i++){
      float minsize=10000000;
      float minlab=0;
      for(int n=-size;n<=size;n++){
        fel=sizes_copy->fel[j+n]+i;
        lel=img_copy->fel[j+n]+i;
        for(int m=-size;m<=size;m++){
          if(minsize>fel[m]){
            minsize=fel[m];
            minlab=lel[m];
            
          }
        }
      } 
      segments->fel[j][i]=minlab;
      sizes->fel[j][i]=minsize;
/*      for(int n=-size;n<=size;n++){
      fel=segments->fel[j+n]+i;
      for(int m=-size;m<=size;m++){
      fel[m]=maxlab;
	  //cout << (int)array[v-1]<< " " ;
    }
    } */
                   
    }
  }  

  
  
  
  delete sizes_copy;
  delete img_copy;
}


void segment_image(DARY *smooth_r, DARY *smooth_g, DARY *smooth_b, float k,  float grad_thres, float bound_thres, int min_size, int &num_ccs, DARY *igrad,  DARY *iori, DARY *labout, vector<segment *> &epoints) {
  int width = smooth_r->x();
  int height = smooth_r->y();
 // build graph
  
  
  
  DARY *dxi= new DARY(height,width,FLOAT1,0.0);
  DARY *dyi= new DARY(height,width,FLOAT1,0.0);
  DARY *dxyi= new DARY(height,width,FLOAT1,0.0);
  DARY *dyxi= new DARY(height,width,FLOAT1,0.0);
  DARY *edgei= new DARY(height,width,FLOAT1,0.0);
  DARY *ssizes= new DARY(height,width,FLOAT1,0.0);
  
  edge *edges = new edge[width*height*4];
  float **igr=igrad->fel;
  float **ior=iori->fel;
  uint num = 0,pa,pb,pc,pd,pe;
  float dx,dy,dyx,dxy,gr,maxori=0;
  for (int y = 1; y < height-1; y++) {
    for (int x = 0; x < width-1; x++) {
      pa = y * width + x;
      pb = y * width + (x+1);
      pc = (y+1) * width + x;
      pd = (y+1) * width + (x+1);
      pe = (y-1) * width + (x+1);
      dx=grad2(smooth_r, smooth_g, smooth_b, pa, pb);
      dy=grad2(smooth_r, smooth_g, smooth_b, pa, pc);
      dxy=grad2(smooth_r, smooth_g, smooth_b, pa, pd);
      dyx=grad2(smooth_r, smooth_g, smooth_b, pa, pe);
       gr=dx+dy+dxy+dyx;
      if(gr<grad_thres){
        dx=0;dy=0;dxy=0;dyx=0;
      }
      dxi->fel[0][pa]=dx;
      dyi->fel[0][pa]=dy;
      dxyi->fel[0][pa]=dxy;
      dyxi->fel[0][pa]=dyx;      
      igr[0][pa]=gr;
      ior[0][pa]=atan2(dy,dx);
      if(ior[0][pa]>maxori)maxori=ior[0][pa];
/*      if(dyx>dxy){
        if(ior[0][pa]<0)ior[0][pa]=ior[0][pa]+M_PI;
        else ior[0][pa]=ior[0][pa]-M_PI;
      }
      while(ior[0][pa]<0)ior[0][pa]+=M_2PI;
      while(ior[0][pa]>M_2PI)ior[0][pa]-=M_2PI;*/
    }
  }    
      
  cannyEdges(dxi, dyi, igrad, edgei,20,40);
 // edgei->writePNG("edge.png"); cout << "edge " << endl;

  float fac;
  for (int y = 1; y < height-1; y++) {
    for (int x = 0; x < width-1; x++) {
      pa = y * width + x;
      pb = y * width + (x+1);
      pc = (y+1) * width + x;
      pd = (y+1) * width + (x+1);
      pe = (y-1) * width + (x+1);
      fac=1;
      if(edgei->fel[0][pa]>0)fac=3;
      igr[0][pa]*=fac;
        edges[num].a = pa;
        edges[num].b = pb;
        edges[num].w = fac*dxi->fel[0][pa];
        num++;
        
        edges[num].a = pa;
        edges[num].b = pc;
        edges[num].w = fac*dyi->fel[0][pa];
        num++;
        
        edges[num].a = pa;
        edges[num].b = pd;
        edges[num].w = fac*dxyi->fel[0][pa];
        num++;
        
        edges[num].a = pa;
        edges[num].b = pe;
        edges[num].w = fac*dyxi->fel[0][pa];
        num++;           
        
    }
  } 
  delete edgei;
  delete dxi;
  delete dyi;
  delete dyxi;
  delete dxyi;
  
  // segment
  universe *u = segment_graph(width*height, num, edges, k);
  //compute extent to area
  vector<edge> bedges;
  // post process small components
  for (uint i = 0; i < num; i++) {
    uint a = u->find(edges[i].a);
    uint b = u->find(edges[i].b);
    float w=edges[i].w;
    if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size))){
     u->join(a, b, w);
    }else if(a != b){
      bedges.push_back(edges[i]);
    }
  }

  uint code,code2,flag;
  vector<uint> codes;
  codes.push_back(0);
  vector<segment> boundaries;
  segment segm;
  boundaries.push_back(segm);       
  uint ori_bin=9;
  for (uint i = 0; i < bedges.size(); i++) {
    uint a = u->find(bedges[i].a);
    uint b = u->find(bedges[i].b);
    if(a != b){
      code=a*num+b;
      code2=b*num+a;
      code=(code<code2)?code:code2;
      flag=0;
      for(uint j=1;j<codes.size() && !flag;j++){// find if the code exists
        if(codes[j]==code)flag=j;
      }
      if(!flag){//code doesn't exist
        flag=codes.size();
        codes.push_back(code);
      }
      if(boundaries.size()>flag){
        boundaries[flag].bpoints.push_back(bedges[i]);
      }else{
        segment seg;
        for(uint o=0;o<ori_bin+1;o++)seg.bweights.push_back(0);
        seg.bpoints.push_back(bedges[i]);
        boundaries.push_back(seg);       
      }
     //boundaries.push_back();
    }
  }
  codes.clear();
  vector<pair<float, int> > worder;
  float ori,maxgori=0;
  for (uint j = 1; j < boundaries.size(); j++) {
    boundaries[j].w=0;
    for (uint i = 0; i < boundaries[j].bpoints.size(); i++) {
      pa=boundaries[j].bpoints[i].a;
      boundaries[j].w+=igr[0][pa];
      ori=ior[0][pa];
      boundaries[j].bweights[(int)(ori_bin*ori/maxori)]+=igr[0][pa];
      //labout->fel[0][boundaries[j].bpoints[i].a]=255;     
    }   
    boundaries[j].w/=boundaries[j].bpoints.size();
    maxgori=0;
    for (uint i = 0; i < boundaries[j].bweights.size(); i++) {
      boundaries[j].bweights[i]/=boundaries[j].bpoints.size();
      if(maxgori<boundaries[j].bweights[i])maxgori=boundaries[j].bweights[i];
    }
    //cout << j <<  " " << maxgori<< " " << boundaries[j].w << " "  << endl;getchar();
    
    worder.push_back(make_pair(maxgori,j));
    //labout->writePNG("lab.png");cout << j << " " << boundaries.size()<< " lab" << endl;getchar();
  }
  std::sort(worder.begin(),worder.end());
  
  for (uint j = 0; j < worder.size(); j++) {
    uint a = u->find(boundaries[worder[j].second].bpoints[0].a);
    uint b = u->find(boundaries[worder[j].second].bpoints[0].b);
    //if(a==b){cout << a<< " " << b << "error "<< endl;getchar();}
    if(worder[j].first<bound_thres && a!=b)u->join(a, b, dx);
    //cout << j <<  " " << worder.size()<< " " << worder[j].first << " " << boundaries[worder[j].second].w << endl;getchar();
  }
  worder.clear();
  //RELABEL THE DATA
  delete [] edges;
  num_ccs = u->num_sets();

  edge *points = new edge[width*height];

  num=0; 
  uint pn;
  for (int y = 2; y < height-2; y++) {
    for (int x = 2; x < width-2; x++) {     
      pn=y * width + x;
      points[num].a=x; 
      points[num].b=y; 
      points[num].w=(float)u->find(pn);
      ssizes->fel[0][pn]=u->size((uint)points[num].w);    
    //if(points[num].a==27 && points[num].b==11)cout <<" " <<  points[num].w << " last " << endl;  
      num++; 
    }
  }  
  //ssizes->writePNG("size.png");
  std::sort(points,points  + num);
  uint lab=1;
  uint last=(uint)points[0].w;
  points[0].w=lab;
  epoints.push_back(new segment());//zero not used
  epoints.push_back(new segment());//first lab=1
  epoints[lab]->points.push_back(points[0]);
  for (uint y = 1; y < num; y++){
    if(points[y].w==last){
      points[y].w=lab;
      epoints[lab]->points.push_back(points[y]);
    }else {
      epoints.push_back(new segment());
      lab++;
      last=(uint)points[y].w;
      points[y].w=lab;
      epoints[lab]->points.push_back(points[y]);
    }
    labout->fel[points[y].b][points[y].a]=lab;
  }
  labout->fel[points[0].b][points[0].a]=points[0].w;//lab;
  
 //openingFilter(labout, ssizes, 1);

  //medianFilter(labout, 1);
  

  
/*  vector<int> neighbor_labels;
  float meanx,meany,varx,vary,covxy;
  for (int y = 1; y < epoints.size(); y++){
    meanx=0;meany=0;
    for(uint p=0;p<epoints[y]->points.size();p++){
      meanx+=epoints[y]->points[p].a;
      meany+=epoints[y]->points[p].b;      
    }   
    meanx/= epoints[y]->points.size();
    meany/= epoints[y]->points.size();
    varx=0;vary=0,covxy=0;
    for(uint p=0;p<epoints[y]->points.size();p++){
      varx+=lsquare(epoints[y]->points[p].a-meanx);
      vary+=lsquare(epoints[y]->points[p].b-meany);
      covxy+=((epoints[y]->points[p].a-meanx)*(epoints[y]->points[p].b-meany));      
    }    
    varx/= epoints[y]->points.size();
    vary/= epoints[y]->points.size();
    covxy/= epoints[y]->points.size();
   
    
  }*/
  
  
  
  //labout->writePNG("labout1.png");getchar();
  float ** labs=labout->fel;
  int labc,labl,labr,labt,labb;
  for (int j = 7; j < height-7; j++){
    for (int i = 7; i < width-7; i++){
      labc=(int)labs[j][i];
      labt=(int)labs[j-5][i];
      labb=(int)labs[j+5][i];
      labl=(int)labs[j][i-5];
      labr=(int)labs[j][i+5];
      if(labt!=labc){
        edge point;
        point.a=i;point.b=j-5;point.w=labt;
        epoints[labt]->points.push_back(point);      
      }    
      if(labl!=labc && labl!=labt){
        edge point;
        point.a=i-5;point.b=j;point.w=labl;
        epoints[labl]->points.push_back(point);             
      }
      if(labr!=labc && labr!=labt && labr!=labl){
        edge point;
        point.a=i+5;point.b=j;point.w=labr;
        epoints[labr]->points.push_back(point);             
      }
      if(labb!=labc && labb!=labt && labb!=labl && labb!=labr){
        edge point;
        point.a=i;point.b=j+5;point.w=labb;
        epoints[labb]->points.push_back(point);             
      }
      
    }
  } 
  epoints.erase(epoints.begin());//remove zero

  
  delete []points;
  delete ssizes;
  //cout << "num " << num_ccs << "  lab " << lab << " epoints " << epoints.size()<< endl;
  delete u;

  
}
