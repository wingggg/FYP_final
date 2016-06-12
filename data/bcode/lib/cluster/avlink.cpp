 //
// C++ Implementation: agglo
//
// Description: 
//
//
// Author: Krystian Mikolajczyk,-,-,- <kma@maui>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "clustering.h"
#include "../util/util.h"



void swapIndexes(int ind, vector<ClStep * > &trace){
  int cf1=trace[ind]->cf1;
  trace[ind]->cf1=trace[ind]->f2;
  trace[ind]->f2=cf1;
  for(uint j=ind+1;j<trace.size();j++){
    if(cf1==trace[j]->cf1){
      trace[j]->cf1=trace[ind]->cf1;
    }else if(cf1==trace[j]->f2){
      trace[j]->f2=trace[ind]->cf1;     
    }
  } 
}

void sortTrace(vector<ClStep *> &trace){
  vector<pair<float,int> > tmp;
  for(uint i=0;i<trace.size();i++){
    tmp.push_back(make_pair(trace[i]->sim,i));
  }
  sort(tmp.begin(),tmp.end());
  
  vector<ClStep *> cltmp;
  for(uint i=0;i<tmp.size();i++){
    cltmp.push_back(trace[tmp[i].second]);
  }
  int flag=1;
  for(uint i=0;i<cltmp.size();i++){
    flag=1;
    for(uint j=i+1;j<cltmp.size() && flag;j++){
      if(cltmp[i]->f2==cltmp[j]->cf1){
	swapIndexes(i, cltmp);
	flag=0;
      }
    }
  }  
  trace.clear();
  trace=cltmp;
}

    

inline float avlinkEuc(float* vec1,float* vec2,int size, float thres){
  float dist=0,d;    
  for(int i=0;i<9;i++){
    d=vec1[i]-vec2[i];
    dist+=d*d;
  }
  if(dist>thres)return dist+thres;
  for(int i=9;i<18;i++){
    d=vec1[i]-vec2[i];
    dist+=d*d;
  }
  if(dist>thres)return dist+thres;
  for(int i=18;i<size;i++){
    d=vec1[i]-vec2[i];
    dist+=d*d;
  }
  return sqrt(dist);
}


void agglomerate( FD* f1, FD *f2)
  /*******************************************************************/
  /* Merge 2 neighboring clusters if compactness is guaranteed.      */
  /*******************************************************************/
{
  /*--------------------------------*/
  /* Compute the new cluster center */
  /*--------------------------------*/  
  /* The combined cluster center can be computed as follows:       */
  /*   c_new = (1/N+M)*(N*c_x + M*c_y)                             */
 
  /*MEAN VAR*/
  int dim=f1->getSize();
  float *vec1=f1->getVec();
  float *vec2=f2->getVec();
  float s1=f1->getNbf();
  float s2=f2->getNbf();
  //cout << f1->getSize() << " " << f2->getSize()<< " " << s1 << " "<< s2 << endl;getchar();
  float *vec=new float[dim];
  for(int i=0;i<dim;i++){
    vec[i]= (s1*vec1[i]+s2*vec2[i])/(s1+s2);
  }
  /*--------------------------*/
  /* Compute the new variance */
  /*--------------------------*/
  /* The new variance can be computed as                           */
  /*   sig_new^2 = (1/N+M)*(N*sig_x^2 + M*sig_y^2 +                */
  /*                        NM/(N+M)*(mu_x - mu_y)^2)              */


  float var=(s1*f1->getVar()+s2*f2->getVar()+((s1*s2)/(s1+s2))*distEuc(vec1,vec2,dim))/(s1+s2);
  //float var=(s1*f1->getVar()+s2*f2->getVar()+((s1*s2)/(s1+s2))*avlinkEuc(vec1,vec2,dim,MAX_FLOAT))/(s1+s2);

  f1->setNbf((int)(s1+s2));

 
  memcpy(vec1,vec,dim*sizeof(float));

  f1->setVar(var);
  delete f2;
  f2=NULL;
  delete []vec; 
}

int getIndex(vector<FD *> features){  
  
  for(uint i=0;i<features.size();i++){
    if(features[i]!=NULL)return i;
  }
  return -1;
}


void  findNearestNeighbor(FD *f1, vector<FD*> features, int &ind, float &sim, int last, DARY *NNdist, DARY *NNvalid, int fullsearch){

  sim=MAX_FLOAT;
  float d=0;
  ind=-1;
  float *vec = f1->getVec();
  float var = f1->getVar();
  int dim = f1->getSize();
  
  if(fullsearch)
    for(uint i=0;i<features.size();i++){
      if(features[i]!=NULL){
	d=var+features[i]->getVar()+distEuc(vec,features[i]->getVec(),dim);
	//d=var+features[i]->getVar()+avlinkEuc(vec,features[i]->getVec(),dim,sim);
	NNdist->fel[last][i]=d;
	NNvalid->bel[last][i]=1;
	if(d<sim){
	  sim=d;
	  ind=i;
	}
      }
    }  
  else{
    for(uint i=0;i<features.size();i++){
      if(NNvalid->bel[last][i]==0 && features[i]!=NULL){
	d=var+features[i]->getVar()+distEuc(vec,features[i]->getVec(),dim);
	//d=var+features[i]->getVar()+avlinkEuc(vec,features[i]->getVec(),dim,sim);
	NNdist->fel[last][i]=d;
	NNvalid->bel[last][i]=1;
      }
    }
    for(uint i=0;i<features.size();i++){
      if(NNdist->fel[last][i]<sim && NNvalid->bel[last][i]==1 && features[i]!=NULL){
	sim=NNdist->fel[last][i];
	ind=i;
      }    
    }
  }  
}


void avlink(vector<FD *> infeatures, vector<ClStep *> &trace, float sim){
  vector<FD *> features;
  for(uint i=0;i<infeatures.size();i++){
    features.push_back(new FD());
    if(infeatures[i]!=NULL){
      features[i]->copy(infeatures[i]);
    }else{
      delete features[i];
      features[i]=NULL;
    }
  }
  int   chsize=100;
  DARY *NNdist= new DARY(chsize,features.size(),FLOAT1);
  DARY *NNvalid= new DARY(chsize,features.size(),UCHAR1);
  bzero(NNvalid->bel[0],chsize*features.size()*sizeof(uchar));

  vector<ClStep * > chain;  
  int indNN;
  float simNN;
  int last=0;
  vector<FD *> fchain;
  int start;
  start=getIndex(features);//get first non null feature
  vector<int> chindex;
  int fullsearch=1;
 
  while(start>=0){ 
    last=0;
    chain.push_back(new ClStep(MAX_FLOAT,start,start));
    fchain.push_back(features[chain[last]->cf1]);//put to the chain and 
    features[chain[last]->cf1]=NULL;//take out of the feature list.
    //int max_size=0;
    fullsearch=1;
    //cout << "building chain " << endl;
    while(chain.size()>0){
      findNearestNeighbor(fchain[last], features, indNN, simNN, last, NNdist, NNvalid, fullsearch);
      if(simNN<chain[last]->sim){// continue if new distance is closer than for pervious pair
        //cout << "building chain " << chain.size()<< endl;
        chain.push_back(new ClStep(simNN,indNN,indNN));
	fchain.push_back(features[indNN]);
	features[indNN]=NULL;
	last++;
	fullsearch=1;
        if(last>=chsize){cout << " ERROR in avlink last>=chsize " << last << endl;exit(0);}
      }else if(chain[last]->sim<sim){//agglomerate if last distance is smaller than
        //cout << "agglomerating cs " << NNvalid->y() << "  fs " <<NNvalid->x()<< "  l-1 " << last-1 << " l " << chain[last-1]->cf1<< endl;
 	//if(chain.size()>max_size)max_size=chain.size();
	agglomerate(fchain[last-1],fchain[last]);
       for(int n=0;n<=last-2;n++)NNvalid->bel[n][chain[last-1]->cf1]=0;      
        fullsearch=0;
	trace.push_back(new ClStep(chain[last]->sim,chain[last-1]->cf1,chain[last]->cf1));
 	features[chain[last-1]->cf1]=fchain[last-1];
        last-=2;  
	fchain.pop_back();fchain.pop_back();
	chain.pop_back();chain.pop_back();
	if(!(trace.size()%100))cout << "\r  "<< trace.size()<< " of " << features.size()<< " current_sim "<< trace[trace.size()-1]->sim<< " of top_sim " << sim <<  " " <<flush;
      } else {
	//deleteDescriptors(fchain); 
	fchain.clear();chain.clear();
      }      
    }
    //if(max_size<chsize)shist[max_size]++;
    start=getIndex(features);    
  }
 
  deleteDescriptors(features);
  sortTrace(trace);
  //cout << endl;for(int i=0;i<chsize;i++)cout << shist[i]<< " "; cout << endl;
  delete  NNdist;
  delete  NNvalid;
}
 
float trace2ClustersComp( vector<FD *> ivDesc, vector<ClStep *> trace,
		        vector<FD *> &vDesc, uint cl_nb){

  uint c=0; 
  uint max_iter_nb=ivDesc.size()-cl_nb;
  cout << "max_iter_nb " << max_iter_nb << " feats " << ivDesc.size() << " cl_nb "<< cl_nb<< endl;
  
  for(uint j=0;j<ivDesc.size();j++){//ADD NEXT TREE LEVEL
    vDesc.push_back(new FD());
    vDesc[j]->features.push_back(ivDesc[j]);
    vDesc[j]->setTreeLevel(1);
  }  
  
  do{
    //cout <<c <<  " "<< trace.size()<< "  " << trace[c]->sim <<  "  " << trace[c]->cf1 << " " <<trace[c]->f2<<  "  "  <<  vDesc[trace[c]->cf1]->features.size()  << "  "  <<  vDesc[trace[c]->f2]->features.size()  << endl;
   vDesc[trace[c]->cf1]->features.insert(vDesc[trace[c]->cf1]->features.end(),
					 vDesc[trace[c]->f2]->features.begin(),
					 vDesc[trace[c]->f2]->features.end());
   delete vDesc[trace[c]->f2];
   vDesc[trace[c]->f2]=NULL;
  }while((++c)<trace.size() && c<max_iter_nb);

  for(uint i=0;i<vDesc.size();i++){
   if(vDesc[i]==NULL){ 
      vDesc.erase(vDesc.begin()+i);
      i--;
   }
  }
  return trace[c]->sim;
}
 
void trace2ClustersSim( vector<FD *> ivDesc, vector<ClStep *> trace,
		        vector<FD *> &vDesc, float dSimThresh){

  uint c=0; 
  for(uint j=0;j<ivDesc.size();j++){//ADD NEXT TREE LEVEL
    vDesc.push_back(new FD());
    vDesc[j]->features.push_back(ivDesc[j]);
    vDesc[j]->setTreeLevel(1);
  }  
  do{
    //cout <<c <<  " "<< trace.size()<< "  " << trace[c]->sim <<  "  " << trace[c]->cf1 << " " <<trace[c]->f2<<  endl;
   vDesc[trace[c]->cf1]->features.insert(vDesc[trace[c]->cf1]->features.end(),
					 vDesc[trace[c]->f2]->features.begin(),
					 vDesc[trace[c]->f2]->features.end());
   delete vDesc[trace[c]->f2];
   vDesc[trace[c]->f2]=NULL;
  }while((++c)<trace.size() && trace[c]->sim<dSimThresh);



  /* uint nb=0;
     for(uint i=0;i<vDesc.size();i++){
     if(vDesc[i]!=NULL)nb++;
     }
     cout << " sim " << c << " " << trace[c]->sim<<  " " << nb << endl;
  */
}


void traceToClusters( vector<FD *> ivDesc, vector<ClStep *> trace,
		      vector<FD *> &vClusters, float dSimThresh ){
  vector<FD *> vDesc;
 trace2ClustersSim( ivDesc, trace, vDesc, dSimThresh);
  for(uint i=0;i<vDesc.size();i++){
    if(vDesc[i]!=NULL){
      vClusters.push_back(vDesc[i]);
    }
  }
  vDesc.clear();
  updateMeanVar(vClusters,1.0);
}



void updateTrace(vector<ClStep *> ktrace, uint cf, vector<ClStep *> tmptrace, uint ksize, vector<ClStep *> &trace){

  vector<uint> index;
  //find which indices from kpart[cf] correspond to global feature indices//
  index.push_back(cf);//index of the cluster
  for(uint i=0;i<ktrace.size();i++){
    if((uint)ktrace[i]->cf1==cf){
      index.push_back(ktrace[i]->f2);
    }
  }
  if(index.size()!=ksize)cout << " problem 1 in update trace"<< endl;
  //combine kmeans trace with agglomerative trace //
  for(uint i=0;i<tmptrace.size();i++){
    int cf1=index[tmptrace[i]->cf1];
    int f2=index[tmptrace[i]->f2];
    trace.push_back(new ClStep(tmptrace[i]->sim, cf1, f2));
  }  
  index.clear();
} 
 
 
void kmeansAgglo(vector<FD *> ivDesc, vector<FD *> &clusters, vector<ClStep *> &trace,  uint K, float dSim){
  cout << "Kmeans-agglomerative clustering, features = "<< ivDesc.size()<< ", sim = "<< dSim << ", K = "<< K << endl;

  vector<ClStep *> ktrace;
  kmeans(ivDesc, ktrace, K);
  vector<FD *> kpart; 
  dSim=dSim ;  
  trace2ClustersSim( ivDesc, ktrace, kpart, MAX_FLOAT);
 

  vector<ClStep *> tmpTrace;
  for(uint i=0; i<kpart.size();i++){
    if(kpart[i]!=NULL){      
      tmpTrace.clear();
      cout << "\navlink on features: "<< kpart[i]->features.size() << endl;
      avlink(kpart[i]->features, tmpTrace, dSim/2.0);
      updateTrace(ktrace, i, tmpTrace, kpart[i]->features.size(), trace);      
    }
  }
  vector<FD *> iagglo;  
  trace2ClustersSim( ivDesc, trace, iagglo, MAX_FLOAT); 
  updateMeanVar(iagglo,1.0);

  cout << endl<< "final avlink on features "<<  ivDesc.size()<< endl;
  avlink(iagglo, trace, MAX_FLOAT);

  traceToClusters( ivDesc, trace, clusters, dSim);
  cout << "done "<< ivDesc.size()<< " " << clusters.size()<< endl;
  
}
 





/*======================================================================*/
/*======================================================================*/
/*======================================================================*/
/*======================================================================*/
/*======================================================================*/
   
void saveTrace(vector<ClStep *> vClusterTrace, const char *filename){
  sortTrace(vClusterTrace);
    ofstream output(filename);
    cout << "saving trace in "<< filename<< endl; 
    if(!output){cout << "error opening " <<filename << endl;exit(0);}
    output << vClusterTrace.size() << endl;
    for(uint i=0; i<vClusterTrace.size();i++){
      output << floor(vClusterTrace[i]->sim) << " " << vClusterTrace[i]->cf1 <<" ";
      output << vClusterTrace[i]->f2  << endl;
    }
    output.close();
}

void saveTrace2(vector<ClStep *> vClusterTrace, const char *filename){
  cout << "saving trace in "<< filename<< endl; 
  //sortTrace(vClusterTrace);
  FILE *fid;
  fid=fopen(filename, "wb");
  if(!fid){cout << "error opening " <<filename << endl;return;}
  uint size=vClusterTrace.size();
  fwrite(&size, sizeof(uint), 1, fid);
  uint *buf = new uint[size*3];
  uint nb=0;
  for(uint i=0; i<vClusterTrace.size();i++){
    buf[nb++]=(uint)floor(vClusterTrace[i]->sim);
    buf[nb++]=(uint) vClusterTrace[i]->cf1;
    buf[nb++]=(uint) vClusterTrace[i]->f2;
  }
  fwrite(buf, sizeof(uint), 3*size, fid);
  delete []buf;
  
  fclose(fid);    
}
void loadTrace2(const char *filename, vector<ClStep *> &vClusterTrace){
    cout << "loading trace from "<< filename<< endl; 
    FILE *fid;
    fid=fopen(filename, "rb");
    uint  nb=0,size;
    fread(&size, sizeof(uint), 1, fid);
    uint *buf = new uint[3*size];
    fread(buf, sizeof(uint), 3*size , fid);
    for(uint i=0; i<size;i++){
      vClusterTrace.push_back(new ClStep(buf[nb],buf[nb+1],buf[nb+2]));
      nb+=3;
    }
    delete []buf;  
    fclose(fid);    

    //sortTrace(vClusterTrace);
}

void loadTrace(const char *filename, vector<ClStep *> &vClusterTrace){
    ifstream input(filename);
    cout << "loading trace from "<< filename << " ..." << endl; 
    if(!input){cout << "error opening " <<filename << endl;exit(0);}
    uint cf1,f2,nb;
    float sim;
    input>>nb;
    for(uint i=0; i<nb;i++){
      input>>sim;input>>cf1;input>>f2;
      vClusterTrace.push_back(new ClStep(sim,cf1,f2));      
    }
    input.close();
    cout << "loaded "<< endl;
    //sortTrace(vClusterTrace);
    
    //saveTrace(vClusterTrace, filename);
}
 
