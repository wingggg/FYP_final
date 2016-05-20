//
// C++ Implementation: kmeans
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

void getClusterMembers(vector<FD *> inClusters, vector<FD *> &outDesc, int level){
  for(uint i=0;i<inClusters.size();i++){
    if(inClusters[i]!=NULL)
    if(inClusters[i]->getTreeLevel() > level){
      //if(inClusters[i]->features[0]->getTreeLevel()==0) cout << "problem " << level<< " " << inClusters[i]->getTreeLevel() << endl;
      getClusterMembers(inClusters[i]->features,outDesc,level);
    }else {

      //if(inClusters[i]->getTreeLevel()!=1) cout << "problem " << level<< " " << inClusters[i]->getTreeLevel() << endl;

      outDesc.push_back(inClusters[i]);
    }
  }
}   

void updateClusterDistances(FD * cluster){
  int dim=cluster->getSize();
  uint nbf=cluster->features.size();
  if(cluster->cluster_dist==NULL){
    cluster->cluster_dist=new DARY(nbf,nbf,FLOAT1,0.0);
  }else if(cluster->cluster_dist->y()!=nbf){ 
    delete cluster->cluster_dist;
    cluster->cluster_dist=new DARY(nbf,nbf,FLOAT1,0.0);
  }

  vector<pair<float,uint> > rad;
  float dist;
  for(uint f=0;f<nbf;f++){
    dist=distEuc(cluster->getVec(),cluster->features[f]->getVec(),dim);
    //dist+=cluster->features[f]->radius;
    cluster->features[f]->setMeanDist(dist);
    rad.push_back(make_pair(dist,f));       
  } 
  sort(rad.begin(),rad.end());
  
  FD *ftmp;
  for(uint f=0;f<nbf;f++){
    ftmp=cluster->features[f];
    cluster->features[f]=cluster->features[rad[f].second];
    cluster->features[rad[f].second]=ftmp;
  } 
  rad.clear();
  for(uint i=0;i<nbf;i++){
    for(uint j=i+1;j<nbf;j++){
      cluster->cluster_dist->fel[i][j]=distEuc(cluster->features[i]->getVec(),cluster->features[j]->getVec(),dim);
      cluster->cluster_dist->fel[j][i]=cluster->cluster_dist->fel[i][j];
    }
  }

}
void updateClusterDistances(vector<FD *> &vClusters){
  if(vClusters.size()==0){
    cout << " error in updateMeanVar " << endl;return;
  }  
  for(uint i=0;i<vClusters.size();i++){
    if(vClusters[i]!=NULL)if(vClusters[i]->features.size()>1 && vClusters[i]->getTreeLevel()>1){
      updateClusterDistances(vClusters[i]);
    }
  }
}


void updateMeanVar(FD * cluster, float quant){
  float radius=0,dist,var,totnbf;
  //cluster->setImageName("cluster");
  cluster->setExtr(cluster->features[0]->getExtr());
  cluster->setType(cluster->features[0]->getType());
  vector<FD*> desc;
  getClusterMembers(cluster->features, desc,0);
  totnbf=(float)desc.size();
  if(quant<0.1 || desc.size()<2){
    cout << "WARNING:updateMeanVar  NBF " << desc.size()<< " " << cluster->features.size() << endl;
   }

  int dim=desc[0]->getSize();
  cluster->setNbf(desc.size());
  cluster->allocVec(dim);  
  float *vec=cluster->getVec();

  float *des;
  //COMPUTE CLUSTER CENTER AND RADIUS
  for(uint f=0;f<desc.size();f++){
    //if(desc[f]==NULL)cout<< f << " NUL"<< endl;
    des=desc[f]->getVec();
    for(int d=0;d<dim;d++){ 
      vec[d]+=(des[d]);
    }
  }     
  for(int d=0;d<dim;d++){
    vec[d]=vec[d]/totnbf;
  }
       
  var=0;
  radius=0;
  vector<pair<float,uint> > rad;
  for(uint f=0;f<desc.size();f++){
    dist=distEuc(vec,desc[f]->getVec(),dim);
    var+=dist;
    rad.push_back(make_pair(dist,f));
  }
  cluster->setVar(var/totnbf);
  sort(rad.begin(),rad.end());
  int ind=(int)(quant*(rad.size()-1));
  //cout << quant << " " << ind << " " << rad.size()<< endl;
  cluster->setRadius((rad[ind].first));

  desc.clear();
  rad.clear();
  
  //cout << "tot max dist "<< totmax << endl;
}
 

void updateMeanVar(vector<FD *> &vClusters, float quant){
  if(vClusters.size()==0 || quant>1.0){
    cout << "ERROR:updateMeanVar " << vClusters.size() << endl;return;
  }  
  for(uint i=0;i<vClusters.size();i++){
    //cout << "\rcomputing mean variance " << i << " of "<<vClusters.size()<< "    "<< flush
    if(vClusters[i]!=NULL)
      if(vClusters[i]->features.size()>1){
	updateMeanVar(vClusters[i],  quant);
      }else if(vClusters[i]->features.size()==1){
	vClusters[i]->copy(vClusters[i]->features[0]); 
	vClusters[i]->setTreeLevel(vClusters[i]->getTreeLevel()+1);
      }
  }
  //cout << "tot max dist "<< totmax << endl;
}
 
 


void kmeansTreeB(vector<FD *> &features, vector<ClStep *> &trace, uint min_size, int K, int max_iter){

  vector<FD *> clusters;  
  kmeans(features, clusters, trace,  K, max_iter, 0);
  //updateMeanVar(clusters,1.0);
  features=clusters; //use lower level for next iteration
  int k;
  for(uint i=0;i<features.size();i++){
    //features[i]->setTreeLevel(2);
    k=(features[i]->features.size()/min_size);
    if(k > K){//if much larger than min_size 
      kmeansTreeB(features[i]->features, trace, min_size, K, max_iter);
    }else if(k>1){//if still larger than min_size 
      kmeansTreeB(features[i]->features, trace, min_size, k, max_iter);      
    }else {//all partitions are smaller than the min_size
      //features[i]->setTreeLevel(1);
      //for(uint j=0;j<features[i]->features.size();j++){
      //	 features[i]->features[j]->setTreeLevel(0);
      // }
    }
  }  
  
}
 

void buildTrace(vector<FD *> &features, vector<ClStep *> &trace){
  
  for(uint i=0;i<features.size();i++){  
    for(uint j=0;j<features[i]->features.size();j++){
      
      if(features[i]->features[j]->getTreeLevel()>=0 && features[i]->features[j]->features.size()){
	buildTrace(features[i]->features, trace);
      }
      if(features[i]->features[j]->getTreeLevel()>=0){
	if(features[i]->getTreeLevel()!=features[i]->features[j]->getTreeLevel()){
	  trace.push_back(new ClStep(features[i]->getNbf(),features[i]->getTreeLevel(),
				     features[i]->features[j]->getTreeLevel()));     

	  //cout <<features[i]->getNbf() <<  " " << features[i]->getTreeLevel()<<" " << features[i]->features[j]->getTreeLevel()<< endl;getchar();
	}
	features[i]->features[j]->setTreeLevel(-1);
      }
    }
  }
}


void kmeansTree(vector<FD *> &features, uint min_size, int K, int max_iter, float quantile){

  //cout <<"before "<<  features.size()<< endl;
  kmeans(features, K, max_iter, 0);
  updateMeanVar(features,quantile);
  for(uint i=0;i<features.size();i++)
    features[i]->setTreeLevel(2);
  
  //  for(uint i=0;i<features.size();i++){
    //  cout <<"after " <<   features[i]->features.size()<<" level " << features[i]->features[0]->getTreeLevel()<< endl;
  // }
  int k;
   
  for(uint i=0;i<features.size();i++){
    k=(features[i]->features.size()/min_size);
    if(k > K){
      kmeansTree(features[i]->features, min_size, K, max_iter, quantile);
    }else if(k>1){
      kmeansTree(features[i]->features, min_size, k, max_iter, quantile);      
    }
  }  
  for(uint i=0;i<features.size();i++){
    if(features[i]->features.size()==1){
      features[i]=features[i]->features[0];
    }
  }
  updateClusterDistances(features);

}


void reduceLevel(FD *&cluster, int &flag){
  if(cluster!=NULL)
    if(cluster->features.size()==1 && cluster->getTreeLevel()>1){
      //cout << "reducing  " << cluster->getTreeLevel() << endl;
      FD *clust=cluster->features[0];
      delete cluster;
      cluster=clust;    
      reduceLevel(cluster, flag);
    }else if(cluster->features.size()==1 && cluster->getTreeLevel()==1){
      //cout << cluster << "killing  " << cluster->features[0]->getTreeLevel() << endl;
      //delete cluster->features[0];
      cluster->features.clear();
      delete cluster;
      cluster=NULL;
      flag=0;
    }
}


void kmeansTree(vector<FD *> &features, vector<ClStep *> &trace, uint min_size, int K, int max_iter){

  for(uint i=0;i<features.size();i++)//remember feature number
    features[i]->setTreeLevel(i);
  kmeansTreeB(features, trace, min_size,  K,  max_iter);//run hierarchical kmeans
  
  cout << "Sorting trace..."<< flush;
  trace.clear();
  buildTrace(features, trace);
  sortTrace(trace);  
  
  cout << " done "<< trace.size()<< endl;
}


//A Tree
void traceToBallTree( vector<FD *> ivDesc, vector<ClStep *> trace,
		      vector<FD *> &vClusters, float minThresh, uint levels, float maxThresh, float quant, uint min_cluster_size){

  float dsim;
  //minThresh=minThresh*minThresh; 
  //cout << "minThres " << minThresh << " maxThres " << maxThresh << " trace max " << trace[trace.size()-2]->sim<< " minSize "<< min_cluster_size<<  endl;
  maxThresh=trace[trace.size()-3]->sim;
  if(levels>0)dsim=(maxThresh-minThresh)/levels;
  else {
    dsim=minThresh;
    maxThresh=minThresh;
  }
  vector<FD *> vDesc;  
  //cout << "OK0" << endl;
  trace2ClustersSim( ivDesc, trace, vDesc, minThresh);
  updateMeanVar(vDesc,1.0);
  //cout << "WE ARE IN "<< endl;
  //vector<FD *> tvDesc;
  // getClusterMembers(vDesc,tvDesc,1);
  //cout << "ones " << vDesc.size()<< endl;
  //for(uint i=0;i<vDesc.size();i++)if(vDesc[i]!=NULL && vDesc[i]->getTreeLevel()!=1)cout << i << " before " << vDesc[i]->features.size() << endl;

  updateClusterDistances(vDesc);


  uint nb=0;
  //cout << "OK1" << endl;
  // REMOVE  N MEMBER CLUSTERS BUT KEEP THE PLACE OTHERWISE INDEXES FROM TRACE WON'T WORK
  for(uint j=0;j<vDesc.size();j++){
    if(vDesc[j]!=NULL)if(vDesc[j]->features.size()<min_cluster_size){
      //delete vDesc[j]->features[0]; //do not delete the bottom features
      vDesc[j]->features.clear();
      delete vDesc[j];
      vDesc[j]=NULL;
    }else {
      nb+=vDesc[j]->features.size();
    }
  } 
  //cout << "OK2 " << nb << endl;
  nb=0;
  uint c=0;
  int lev=1;
  while((++c)<trace.size() && trace[c]->sim<minThresh);//move c to the level of appearance clusters
  //BUILD BALL TREE
  float thres=minThresh;
  while(c<trace.size() && trace[c]->sim<maxThresh){
    
    if(trace[c]->sim>thres){
      updateMeanVar(vDesc,quant); //COMPUTE NODE CENTERS
      updateClusterDistances(vDesc);

     thres+=dsim;lev++;
      for(uint j=0;j<vDesc.size();j++){//ADD NEXT TREE LEVEL
	if(vDesc[j]!=NULL)if(vDesc[j]->features.size()>1 || lev==2){
	  FD *clust = new FD();
	  clust->features.push_back(vDesc[j]);
	  clust->setTreeLevel(lev);
	  vDesc[j]=clust;
	}//else vDesc[j]->setTreeLevel(lev);
      }
     }
    
    if(vDesc[trace[c]->f2]!=NULL && vDesc[trace[c]->cf1]!=NULL){//merge if non empty
      vDesc[trace[c]->cf1]->features.insert(vDesc[trace[c]->cf1]->features.end(),
					    vDesc[trace[c]->f2]->features.begin(),
					    vDesc[trace[c]->f2]->features.end());
      vDesc[trace[c]->f2]->features.clear();     
      delete vDesc[trace[c]->f2];
      vDesc[trace[c]->f2]=NULL;
    }else if(vDesc[trace[c]->f2]!=NULL){
      vDesc[trace[c]->cf1]=vDesc[trace[c]->f2];
      vDesc[trace[c]->f2]=NULL;
    }
    c++;        
  }


  // cout << "OK3 " << vDesc.size()<< endl;

  nb=0;
  int flag=1;
  for(uint j=0;j<vDesc.size();j++){
    reduceLevel(vDesc[j], flag);
    if(vDesc[j]!=NULL && flag){
      updateMeanVar(vDesc[j],quant);
      updateClusterDistances(vDesc[j]);
      if(vDesc[j]->features.size()==1)cout << "problem  " <<j << " " << vDesc[j]->features[0]->features.size()<< endl; 
      nb+=vDesc[j]->getNbf();
      vDesc[j]->setMeanDist(0);//top nodes have no center.
      if(vDesc[j]->getTreeLevel()>1){
	vClusters.push_back(vDesc[j]);
      }
    }
    flag=1;
  }  
  //cout << "OK4" << endl;

  if(vClusters.size()==1){
    cout <<"Tree level " <<  vClusters[0]->getTreeLevel()<< endl;
    vClusters.insert(vClusters.begin(),vClusters[0]->features.begin(),vClusters[0]->features.end());
    vClusters.pop_back();
  }
  //cout <<"Tree level " <<  vClusters[0]->getTreeLevel()<< endl;
  cout << vClusters.size()<<" " << vClusters[0]->features.size()<<  " top clusters from "<< nb << " features. " <<endl;

  for(uint i=0;i<vClusters.size();i++){
    if(vClusters[i]->getTreeLevel()<2)vClusters[i]->Cout(128);
  }
}



//K Tree
void traceToBallTree( vector<FD *> ivDesc, vector<ClStep *> trace,
		      vector<FD *> &vClusters, float minThresh,  float quant, uint min_cluster_size){

  cout << "min " << minThresh << " trace max " << trace[trace.size()-1]->sim<< endl;
  vector<FD *> vDesc;  

  for(uint i=0;i<ivDesc.size();i++){
    ivDesc[i]->setTreeLevel(0);
    vDesc.push_back(ivDesc[i]);
  }

  uint c=0;   
  FD *clust;
  clust = new FD();
  clust->copy(ivDesc[trace[0]->cf1]);
  clust->features.push_back(vDesc[trace[0]->cf1]);
  clust->setTreeLevel(vDesc[trace[c]->cf1]->getTreeLevel()+1);
  vDesc[trace[0]->cf1]=clust;
  vDesc[trace[0]->cf1]->features.push_back(vDesc[trace[0]->f2]);
  vDesc[trace[0]->f2]=NULL;
  c=1;
  //BUILD BALL TREE
  while(c<trace.size()){
    if(trace[c]->cf1!=trace[c-1]->cf1){//Add new level  
      if(vDesc[trace[c-1]->cf1]->features.size()>1){
	updateMeanVar(vDesc[trace[c-1]->cf1],quant); //COMPUTE NODE CENTERS
	updateClusterDistances(vDesc[trace[c-1]->cf1]);
	
      }
      clust = new FD();
      clust->copy(vDesc[trace[c]->cf1]);
      clust->features.push_back(vDesc[trace[c]->cf1]);
      clust->setTreeLevel(vDesc[trace[c]->cf1]->getTreeLevel()+1);      
      //if(vDesc[trace[c]->f2]!=NULL){
	clust->features.push_back(vDesc[trace[c]->f2]);
	vDesc[trace[c]->f2]=NULL;
	//}
      vDesc[trace[c]->cf1]=clust;
    }else if(vDesc[trace[c]->f2]!=NULL){//merge nodes
      vDesc[trace[c]->cf1]->features.push_back(vDesc[trace[c]->f2]);
      vDesc[trace[c]->f2]=NULL;
    }
    c++;        
  }
  uint nb=0;
  for(uint j=0;j<vDesc.size();j++){
    if(vDesc[j]!=NULL){      
      updateMeanVar(vDesc[j],quant);   //cout << "ok2"<< endl;
      updateClusterDistances(vDesc[j]);
      if(vDesc[j]->features.size()==1)cout << "problem  " <<j << " " << vDesc[j]->features[0]->features.size()<< endl; 
      //cout << vDesc[j]->getNbf()<< endl;
      nb+=vDesc[j]->getNbf();
      vDesc[j]->setMeanDist(0);//top nodes have no center.
      //cout <<"Tree mean dist " << vDesc[j]->getMeanDist()<< endl;
      vClusters.push_back(vDesc[j]);
    }
  }  
 

  cout << "top nodes "<< vClusters.size()<<" tree level " << vClusters[0]->getTreeLevel()<<  " clusters from "<< nb << " features."<< endl;

}


//A tree old
void traceToBallTree( vector<FD *> ivDesc, vector<ClStep *> trace,
		     vector<FD *> &vClusters, float minThresh, uint memNb, float quant){
 
  vector<FD *> vDesc;  
  trace2ClustersSim( ivDesc, trace, vDesc, minThresh);
  updateMeanVar(vDesc,1.0);
  

  // REMOVE  ONE MEMBER CLUSTERS BUT KEEP THE PLACE OTHERWISE INDEXES FROM TRACE WON'T WORK
  for(uint j=0;j<vDesc.size();j++){
    if(vDesc[j]!=NULL)if(vDesc[j]->features.size()<2 ){
      //delete vDesc[j]->features[0]; //do not delete the bottom features
      vDesc[j]->features.clear();
      delete vDesc[j];
      vDesc[j]=NULL;
    }else {   //ADD NEXT TREE LEVEL
      FD *clust = new FD();
      clust->features.push_back(vDesc[j]);
      clust->setTreeLevel(2);
      vDesc[j]=clust;
    }
    
  } 

  uint c=0;
  while((++c)<trace.size() && trace[c]->sim<minThresh);
  //BUILD BALL TREE
  while(c<trace.size()){
  
    if(vDesc[trace[c]->f2]!=NULL && vDesc[trace[c]->cf1]!=NULL){
      vDesc[trace[c]->cf1]->features.insert(vDesc[trace[c]->cf1]->features.end(),
					    vDesc[trace[c]->f2]->features.begin(),
					    vDesc[trace[c]->f2]->features.end());
      delete vDesc[trace[c]->f2];
      vDesc[trace[c]->f2]=NULL;
    }else if(vDesc[trace[c]->f2]!=NULL){
      vDesc[trace[c]->cf1]=vDesc[trace[c]->f2];
      vDesc[trace[c]->f2]=NULL;
      cout <<"it should not happend in traceToBallTree."<< endl;
    }
    
    if(vDesc[trace[c]->cf1]!=NULL && vDesc[trace[c]->cf1]->features.size()>=memNb){//CREATE NEW NODE
      updateMeanVar(vDesc[trace[c]->cf1],quant); 
      FD *clust = new FD();
      clust->features.push_back(vDesc[trace[c]->cf1]);
      clust->setTreeLevel(vDesc[trace[c]->cf1]->getTreeLevel()+1);
      vDesc[trace[c]->cf1]=clust;
    }
    c++;        
  }

  uint nb=0;
  for(uint j=0;j<vDesc.size();j++){
    if(vDesc[j]!=NULL){
      if(vDesc[j]->features.size()==1){//REMOVE FROM TOP IF IT HAS ONLY ONE MEMBER
	FD *clust=vDesc[j]->features[0];
	delete vDesc[j];
	vDesc[j]=clust;
      }   
      updateMeanVar(vDesc[j],quant);
      if(vDesc[j]->features.size()==1)cout << "problem  " <<j << " " << vDesc[j]->features[0]->features.size()<< endl; 
      nb+=vDesc[j]->getNbf();
      vClusters.push_back(vDesc[j]);
    }
  }  
  if(vClusters.size()==1){
    cout <<"Tree level " <<  vClusters[0]->getTreeLevel()<< endl;
    vClusters.insert(vClusters.begin(),vClusters[0]->features.begin(),vClusters[0]->features.end());
    vClusters.pop_back();
  } 
    cout <<"Tree level " <<  vClusters[0]->getTreeLevel()<< endl;

  cout << vClusters.size()<<" " << vClusters[0]->features.size()<<  " top clusters from "<< nb << " features."<< endl;
}
   


void saveBallTreeBin(vector<FD *> clusters,  FILE *fid, uint &tot, float *buf){ 
  for(uint i=0; i<clusters.size();i++){
    
    uint size=clusters[i]->pos_weight.size();
    fwrite(&size, sizeof(uint), 1, fid);
    float w;
    for(uint j=0;j<clusters[i]->pos_weight.size();j++){
      w=clusters[i]->pos_weight[j];
      fwrite(&w, sizeof(float), 1, fid);
    }
    size=clusters[i]->features.size();
    fwrite(&size, sizeof(uint), 1, fid);
    clusters[i]->writeBin(fid, buf);tot++;
    //if(clusters[i]->getTreeLevel()>1 && size<2)cout << "OK1 "<<size<< " " <<clusters[i]->getTreeLevel()<<  endl;
    if(clusters[i]->getTreeLevel()>1){
      //cout << "OK1 "<<size<< " " <<clusters[i]->getTreeLevel()<< " " << clusters[i]->cluster_dist->size()<< endl;
      fwrite(clusters[i]->cluster_dist->fel[0], sizeof(float), clusters[i]->cluster_dist->size(), fid);
    }
    //cout << "OK2 "<< endl;
    if(size>0){
      saveBallTreeBin(clusters[i]->features, fid,tot,buf);
    }
  }  
}


void saveOccurrences(vector<FD *> clusters, FILE *fid){  
}



void saveBallTreeBinary(vector<FD *> clusters,  uint fnb, const char *filename){
  cout << "saving ball tree structure in "<< filename << ", top nodes: "  <<clusters.size()<< flush;  
  vector<FD *> cl;
  getClusterMembers(clusters,cl,1);//GET POINTERS TO APPEARANCE CLUSTERS
  for(uint i=0; i<cl.size();i++){
   deleteDescriptors(cl[i]->features);    
  }  
 
  uint tot=0;
  FILE *fid;
  fid=fopen(filename, "wb");
  if(!fid){cout << "error opening " <<filename << endl;return;}

  
  uint size=clusters[0]->getSize();fwrite(&size, sizeof(uint), 1, fid);
  uint dim=size;
  size=clusters.size();fwrite(&size, sizeof(uint), 1, fid);
  size=fnb;fwrite(&size, sizeof(uint), 1, fid);
  float *buf = new float[200];
  saveBallTreeBin(clusters, fid,tot,buf);
  delete []buf;
  cout << ", all nodes: " << tot<< endl;

  //SAVE OCCURRENCES
  uint nb=0;
  uint cl_nb=cl.size();
  for(uint i=0; i<cl.size();i++){
    nb+= cl[i]->occurences.size();
  }
  fwrite(&nb, sizeof(uint), 1, fid);
  fwrite(&cl_nb, sizeof(uint), 1, fid);
  if(nb>0){
    buf = new float[nb*(OPARAMS_NB+1)];//12-number of parameters to save
    uint pos=0;
    for(uint i=0; i<cl.size();i++){
      for(uint f=0; f<cl[i]->occurences.size();f++){
	buf[pos++]=i;
	for(uint p=0; p<OPARAMS_NB;p++){
	  buf[pos++]=cl[i]->occurences[f]->par[p];
	}
      }  
    }
    fwrite(buf, sizeof(float), nb*(OPARAMS_NB+1), fid);
    delete []buf;
  }
  fclose(fid);    
  cout << ", clusters: "<< cl.size()<< ", relations: "<< nb << ",dim: " << dim<< ", type: " << cl[0]->getType()<<  endl;

}



//static float maxweight=0;
void loadBallTreeBin(vector<FD *> &clusters, FILE *fid, int size, uint &tot, float *buf){ 
  
    FD *cor = new FD();
    uint nb;
    fread(&nb, sizeof(uint), 1, fid);
    float w;
    for(uint i=0;i<nb;i++){
      fread(&w, sizeof(float), 1, fid);
      cor->pos_weight.push_back(w);
      //      if(maxweight<w)maxweight=w;
    }
    fread(&nb, sizeof(uint), 1, fid);
    cor->readBin(fid,size,buf);
    tot++;
    //cout << " OKnb "<< nb << " " << cor->getTreeLevel()<< endl;
    if(cor->getTreeLevel()>1){
      cor->cluster_dist=new DARY(nb,nb,FLOAT1);
      //cout << "OK1 "<<nb<< " " <<cor->getTreeLevel()<< " " << cor->cluster_dist->size()<< endl;
     fread(cor->cluster_dist->fel[0], sizeof(float), cor->cluster_dist->size(), fid);    
    }
    //cout << "OK "<< endl;
    for(uint i=0; i<nb;i++){
      //cout<<i<<  " nb "<< nb << " " << tot<< " " << endl;getchar();
      loadBallTreeBin(cor->features, fid, size, tot, buf);
    }
    clusters.push_back(cor);
}


void loadBallTreeBinary(const char *filename, vector<FD *> &clusters, uint &fnb){
    cout << "loading tree structure from "<< filename << flush;
    FILE *fid;
    fid=fopen(filename, "rb");
    uint  nb,size,tot=0;
    fread(&size, sizeof(uint), 1, fid);
    fread(&nb, sizeof(uint), 1, fid);
    fread(&fnb, sizeof(uint), 1, fid);
    float *buf = new float[200];
    if(nb==0){cout << "reading error in "<<filename<< endl;return;}
    for(uint i=0;i<nb;i++){
      loadBallTreeBin(clusters, fid, size, tot,buf);
    }
    delete []buf;
    cout << ", top nodes: "<< clusters.size() << ", all nodes: "<< tot << endl;

    vector<FD *> cl;
    getClusterMembers(clusters,cl,1);//GET POINTERS TO APPEARANCE CLUSTERS
    uint  cln,cl_nb;
    fread(&nb, sizeof(uint), 1, fid);
    fread(&cl_nb, sizeof(uint), 1, fid);
    if(nb==0 || cl_nb!=cl.size()){cout << "reading error in "<<filename<< endl; fclose(fid);return;}

    buf = new float[nb*(OPARAMS_NB+1)];
    fread(buf, sizeof(float), nb* (OPARAMS_NB+1), fid);
    uint pos=0;
    Occurence *occ;
    for(uint i=0;i<nb;i++){
      cln=(uint)buf[pos++];
      occ = new Occurence();
      for(uint p=0;p<OPARAMS_NB;p++){
	occ->par[p]=buf[pos++];
      }      
      cl[cln]->occurences.push_back(occ);
    }
    delete []buf;
    //cout << "  class1 " << cl1 << " class2 " << cl2<< endl;
    cout << ", clusters: "<< cl.size()<< ", relations: "<< nb << ", dim: "<< size << ", type: " << cl[0]->getType()<< endl;

    // cout << ",maxweight: "<< maxweight << endl;
    fclose(fid);

}


void searchBallTreeNN(FD *desc, vector<FD*>  nodes, float mdist, float **cl_dist, float &dThres, int level, int desc_dimension, vector<FD* > &matches, vector<float > &distances, float *n_dist, int *i_dist){  
  float dist;
  float *vec=desc->getVec();
  int  flag=1;
  //int s_dist=0;

  //float *ln_dist=new float[nodes.size()];
  //int *li_dist=new int[nodes.size()];
  for(uint i=0;i<nodes.size();i++){
    flag=0;
    //cout << mdist << " m-t  " << mdist-dThres << " <?" << nodes[i]->getMeanDist() << endl;
    if(fabs(mdist-nodes[i]->par[MEAN_DIST])-nodes[i]->par[RADIUS] < dThres)flag=1;    
    //for(int d=0;d<s_dist && flag;d++){
      //cout << "n-t "<< n_dist[d]-dThres<< " <? " << cl_dist[i][i_dist[d]]+nodes[i]->par[RADIUS]<<" i " << i << " r " <<nodes[i]->par[RADIUS]<<   endl;
      //if(ln_dist[d]-dThres>cl_dist[i][li_dist[d]]+nodes[i]->par[RADIUS])flag=0;	      	
      //if(ln_dist[d]-dThres>cl_dist[i][li_dist[d]])flag=0;	      	
    //}//getchar();
    if(flag){
      dist=distEuc(vec,nodes[i]->getVec(),desc_dimension);
      //ln_dist[s_dist]=dist;
      //li_dist[s_dist]=i; 
      //s_dist++;           
      if(dist-nodes[i]->par[RADIUS]<=dThres){
	if(nodes[i]->getTreeLevel()>level){ 
	  searchBallTreeNN(desc, nodes[i]->features, dist, nodes[i]->cluster_dist->fel, dThres, level, desc_dimension, matches,distances,n_dist,i_dist);	
	}else if(dist<dThres){
	  matches.clear();
	  distances.clear();
	  dThres=dist;
	  distances.push_back(dist);
	  matches.push_back(nodes[i]);
	}
      }    
    }   
  } 
  
  //delete []ln_dist;
  //delete []li_dist;
}

void searchBallTreeNNT(vector<FD* > desc, vector<FD* > nodes, float dThres, int level, int desc_dimension, vector<FD* > &matches, vector<float > &distances){  
  float matching_threshold,dist;
  vector<FD* >  match;
  vector<float>  nndist;

  DARY *cl_dist=new DARY(nodes.size(),nodes.size(),FLOAT1,MAX_FLOAT);
  float *n_dist=new float[200];
  int *i_dist=new int[200];
  for(uint j=0;j<desc.size();j++){
    matching_threshold=dThres;
    if(!(j%100))cout << "\r feature " << j << "   "<< flush;
    for(uint i=0;i<nodes.size();i++){
      dist=distEuc(desc[j]->getVec(),nodes[i]->getVec(),desc_dimension);
      if(dist-nodes[i]->par[RADIUS]<dThres){	
	searchBallTreeNN(desc[j], nodes[i]->features, dist, cl_dist->fel, matching_threshold, level, desc_dimension , match, nndist,n_dist,i_dist); 
      }
    }
    if(match.size()>0){
      matches.push_back(match[0]);
      distances.push_back(nndist[0]);
      nndist.clear();
      match.clear();
    }    
  } 
  delete cl_dist;
  delete []n_dist;
  delete []i_dist;

}

 //method nb 1
void searchBallTree(FD *desc, vector<FD* > nodes, float mdist,  float dThres, int level, int desc_dimension, vector<FD* > &matches1,vector<FD* > &matches2, vector<float > &distances){  
  float dist;
  float *vec=desc->getVec();
  for(uint i=0;i<nodes.size();i++){
    //cout <<i <<  " rad " << nodes[i]->par[RADIUS] << " mean_dist " << nodes[i]->mean_dist   <<  " mdist " << mdist <<endl;
    {//fabs(mdist-nodes[i]->mean_dist)-nodes[i]->par[RADIUS]<dThres){
      //cout <<fabs(mdist-nodes[i]->mean_dist)-nodes[i]->par[RADIUS]<< " < " << dThres << endl;
      //getchar();
      dist=distEuc(vec,nodes[i]->getVec(),desc_dimension);
      if(dist-nodes[i]->par[RADIUS]<=dThres){
	if(nodes[i]->getTreeLevel()>level){	
	  searchBallTree(desc, nodes[i]->features, dist, dThres, level, desc_dimension, matches1,matches2,distances);	
	}else if(dist<dThres){
	  distances.push_back(dist);
	  matches1.push_back(desc);
	  //desc->Cout(5);
	  matches2.push_back(nodes[i]);
	}
      }    
    }   
  }  
}

void searchBallTree(FD *desc, vector<FD* > nodes, float mdist,  float dThres, int level, int desc_dimension, vector<FD* > &matches, vector<float > &distances){  
  float dist;
  float *vec=desc->getVec();
  for(uint i=0;i<nodes.size();i++){
    //cout <<i << " rad " << nodes[i]->par[RADIUS] << " mean_dist " << nodes[i]->mean_dist   <<  " mdist " << mdist <<endl;
    if(1 || fabs(mdist-nodes[i]->par[MEAN_DIST])-nodes[i]->par[RADIUS]<dThres){
      //cout <<fabs(mdist-nodes[i]->par[MEAN_DIST])-nodes[i]->par[RADIUS]<< " < " << dThres << endl;
      //getchar();
      dist=distEuc(vec,nodes[i]->getVec(),desc_dimension);
      if(dist-nodes[i]->par[RADIUS]<=dThres){
	if(nodes[i]->getTreeLevel()>level){
	  searchBallTree(desc, nodes[i]->features, dist, dThres, level, desc_dimension, matches, distances);	
	}else if(dist<dThres){
	  //if(nodes[i]->getTreeLevel()<1)cout << level<< " "<< dist<< endl;
	  distances.push_back(dist);
	  matches.push_back(nodes[i]);
	}
      }    
    }   
  }  
}


void searchBallTree(FD*  desc, vector<FD* > nodes, float dThres, int level, int desc_dimension, vector<FD* > &matches, vector<float > &distances){  
  float dist;
  for(uint i=0;i<nodes.size();i++){
    dist=distEuc(desc->getVec(),nodes[i]->getVec(),desc_dimension);
    if(dist-nodes[i]->par[RADIUS]<dThres){
      searchBallTree(desc, nodes[i]->features, dist, dThres, level, desc_dimension , matches, distances); 
    }
  }  
}


//calls method nb 1 up 
void searchBallTree(vector<FD* > desc, vector<FD* > nodes, float dThres, int level, int desc_dimension, vector<FD* > &matches1, vector<FD* > &matches2, vector<float > &distances){  
  float dist;
  for(uint j=0;j<desc.size();j++){
    if(!(j%5000)){cout << "\rmatching feature " << j << " with   "<<matches1.size()<< " occurrences   "<< flush;}
    for(uint i=0;i<nodes.size();i++){
      dist=distEuc(desc[j]->getVec(),nodes[i]->getVec(),desc_dimension);
      if(dist-nodes[i]->par[RADIUS]<dThres){
	//if(nodes[i]->getTreeLevel()<2)cout << nodes[i]->features.size()<< " test " << nodes[i]->features[0]->getTreeLevel()<< endl;
	searchBallTree(desc[j], nodes[i]->features, dist, dThres, level, desc_dimension , matches1, matches2, distances); 
      }
    }
  } 
  cout  << " ...finished "<< endl;
}
 
void searchBallTreeNN2(FD *desc, vector<FD* > nodes, float &dThres, int level, int desc_dimension,  FD* match, float  &distance){  
  float dist;
  int first=0, second=0;
  float min1=MAX_FLOAT,min2=MAX_FLOAT;
  for(uint i=0;i<nodes.size();i++){
    dist=distEuc(desc->getVec(),nodes[i]->getVec(),desc_dimension);
    if(dist<min2){
      second=first;
      min2=min1;
    }else if(dist<min1){
      first=i;
      min1=dist;
    }
  }

 
    if(min1-nodes[first]->par[RADIUS]<=dThres){
      if(nodes[first]->getTreeLevel()>level){
	searchBallTreeNN2(desc, nodes[first]->features,  dThres, level, desc_dimension, match,distance);	
      }else if(min1<dThres){
	dThres=min1;
	match=nodes[first];
	distance=min1;
      }
    }  
    if(second!=first && min2-nodes[second]->par[RADIUS]<=dThres){
      if(nodes[second]->getTreeLevel()>level){
	searchBallTreeNN2(desc, nodes[second]->features,  dThres, level, desc_dimension, match,distance);	
      }else if(min2<dThres){
	dThres=min2;
	match=nodes[second];
	distance=min2;
      }
    }  
} 


void searchBallTreeNN(FD *desc, vector<FD* > nodes, float &dThres, int level, int desc_dimension,  vector<FD* > &match, vector<float > &distances){  
  float dist;
  for(uint i=0;i<nodes.size();i++){
    dist=distEuc(desc->getVec(),nodes[i]->getVec(),desc_dimension);
    if(dist-nodes[i]->par[RADIUS]<=dThres){
      if(nodes[i]->getTreeLevel()>level){
	searchBallTreeNN(desc, nodes[i]->features,  dThres, level, desc_dimension, match,distances);	
      }else if(dist<dThres){
	dThres=dist;
	match.clear();
	distances.clear();
	match.push_back(nodes[i]);
	distances.push_back(dist);
      }
    }      
  } 
} 

void searchBallTreeNN(vector<FD* > desc, vector<FD* > ballTree, float dThres, int level, int desc_dimension, vector<FD* > &match1, vector<FD* > &match2, vector<float > &distances){  
  float matching_threshold;
  vector<FD* > nnmatch;
  vector<float> dis;
  
  for(uint j=0;j<desc.size();j++){
    //if(!(j%100))cout << "\r feature " << j << "   "<< flush;
    matching_threshold=dThres;
    searchBallTreeNN(desc[j], ballTree, matching_threshold, level, desc_dimension , nnmatch, dis); 
    match1.push_back(desc[j]);
    match2.push_back(nnmatch[0]);
    distances.push_back(dis[0]);        
  } 
}


void searchBallTreekNN(FD *desc, vector<FD* > nodes, float &dThres, int k, int level, int desc_dimension, float mdist, vector<FD* > &match, vector<FD* > &kNN,vector<float > &kNNdist){  
  float dist;
  for(uint i=0;i<nodes.size();i++){
    if(fabs((mdist)-(nodes[i]->par[MEAN_DIST]))-(nodes[i]->par[RADIUS])<(dThres))
      {
      dist=distEuc(desc->getVec(),nodes[i]->getVec(),desc_dimension);
      if((dist)-(nodes[i]->par[RADIUS])<=(dThres)){
	if(nodes[i]->getTreeLevel()>level){
	  searchBallTreekNN(desc, nodes[i]->features,  dThres,k, level, desc_dimension, dist, match, kNN, kNNdist);	
	}else if(dist<kNNdist[k]){
	  dThres=kNNdist[k];
	  int f=0;
	  while(dist>kNNdist[f] && f<=k){	  
	    //cout << f << " " << dist <<  " "<< kNNdist[f]<< endl; 
	    f++;
	  }
	  if(f<k){
	    //cout << f << " " << dist <<  " "<< kNNdist[0]<<  " "<< kNNdist[1]<<  " "<< endl;
	    kNN.insert(kNN.begin()+f,nodes[i]);
	    kNNdist.insert(kNNdist.begin()+f,dist);
	    kNN.pop_back();
	    kNNdist.pop_back();
	    //cout << f << " " << dist <<  " "<< kNNdist[0]<<  " "<< kNNdist[1]<<  " "<< endl; getchar();
	  }
	}
      }      
    } 
  }
}   
    
    
void searchBallTreekNN(vector<FD* > desc, vector<FD* > ballTree, float dThres, int k, int level, int desc_dimension, vector<FD* > &match1, vector<FD* > &match2,vector<float > &distances){  

  float matching_threshold;
  vector<FD* > qmatch;
  vector<FD* > kNN;
  vector<float> kNNdist;
  
  for(int i=0;i<k+1;i++){
    kNN.push_back(new FD());
    kNNdist.push_back(MAX_FLOAT);
  }
  for(uint j=0;j<desc.size();j++){
    //if(!(j%100))cout << "\r feature " << j << "   "<< flush;
    for(int i=0;i<k+1;i++){
      kNNdist[i]=MAX_FLOAT;
    }
    matching_threshold=dThres;
    searchBallTreekNN(desc[j], ballTree, matching_threshold, k, level, desc_dimension , 0.0, qmatch, kNN, kNNdist);
    for(int i=0;i<k;i++){
      match1.push_back(desc[j]);
      match2.push_back(kNN[i]);
      distances.push_back(kNNdist[i]);        
    }
  } 
  kNN.clear();	    

  kNNdist.clear();
}





void searchBallTreeNisterT(FD *desc, vector<FD* > nodes, float mdist, float **cl_dist, float dThres, int level, int desc_dimension, vector<FD* > &matches, vector<float > &distances, float *n_dist, int *i_dist){  
  float dist,min_dist=dThres*100;
  FD *min_node=NULL;
  int  flag=1;
  //int s_dist=0;
  for(uint i=0;i<nodes.size();i++){
    flag=1;
    //cout << i << " " << mdist << " m-t  " << mdist-dThres << " >?" << nodes[i]->par[MEAN_DIST] << " " << nodes.size()<< " lev " << nodes[i]->getTreeLevel()<< endl;
    if(mdist-min_dist > nodes[i]->par[MEAN_DIST])flag=0;    
    //for(int d=0;d<s_dist && flag;d++){
      //cout << d<< " d "<< n_dist[d]<< "-"<< min_dist<< " >? " << cl_dist[i][i_dist[d]] << " r " <<nodes[i]->par[RADIUS]<<   endl;
      //if(n_dist[d]-min_dist > cl_dist[i][i_dist[d]])flag=0;	      	
      //}//getchar();
    if(flag){
      dist=distEuc(desc->getVec(),nodes[i]->getVec(),desc_dimension);
      //n_dist[s_dist]=dist;
      //i_dist[s_dist]=i; 
      //s_dist++;           
      if(dist<min_dist){
	min_dist=dist;
	min_node=nodes[i];
      }
    }
  }
  //cout << "OK "<< min_dist<< " " << min_node->cluster_dist->fel[0][2]<<endl;
  if(min_node!=NULL){
    if(min_node->getTreeLevel()>level){ 
      //cout << "Next Level "<< min_node->cluster_dist->fel[0][2]<<endl;
      searchBallTreeNisterT(desc, min_node->features, min_dist, min_node->cluster_dist->fel, dThres, level, desc_dimension, matches,distances, n_dist, i_dist);	
    }else if(min_dist<dThres){
      distances.push_back(min_dist);
      matches.push_back(min_node);
      min_node=NULL;
    }    
  }  
}

void searchBallTreeNister(FD *desc, vector<FD* > nodes, float dThres, int level, int desc_dimension, vector<FD* > &matches, vector<float > &distances){  
  float dist,min_dist=100*dThres; 
  FD *min_node=NULL;

  for(uint i=0;i<nodes.size();i++){
    dist=distEuc(desc->getVec(),nodes[i]->getVec(),desc_dimension);
    if(dist<min_dist){
      min_dist=dist;
      min_node=nodes[i];
    }
  }
  if(min_node!=NULL){
    if(min_node->getTreeLevel()>level){ 
      searchBallTreeNister(desc, min_node->features,  dThres, level, desc_dimension, matches,distances);	
    }else if(min_dist<dThres){
      distances.push_back(min_dist);
      matches.push_back(min_node);
      min_node=NULL;
    }    
  }
}

void searchBallTreeNister(vector<FD* > desc, vector<FD* > ballTree, float dThres, int level, int desc_dimension, vector<FD* > &matches, vector<float > &distances){  

  DARY *cl_dist=new DARY(ballTree.size(),ballTree.size(),FLOAT1,MAX_FLOAT);
  float *n_dist=new float[200];
  int *i_dist=new int[200];
  for(uint j=0;j<desc.size();j++){
      searchBallTreeNister(desc[j], ballTree, dThres, level, desc_dimension , matches,distances); 
      //searchBallTreeNisterT(desc[j], ballTree, 0, cl_dist->fel, dThres, level, desc_dimension , matches,distances,n_dist,i_dist);  
  } 
  delete cl_dist;
  delete []n_dist;
  delete []i_dist;
}



void searchExhaustive(vector<FD* > desc, vector<FD* > ballTree, float dThres, int desc_dimension, vector<FD* > &matches1,vector<FD* > &matches2, vector<float > &distances) {
  
  float dist;
  float *d;
  //dThres=dThres*dThres;
  for(uint j=0;j<desc.size();j++){
    d=desc[j]->getVec();
    if(!(j%100))cout << "\r feature " << j << "   "<< flush;
    for(uint i=0;i<ballTree.size();i++){
      dist=distEuc(d,ballTree[i]->getVec(),desc_dimension);
      if(dist<dThres){
	distances.push_back(dist);
	matches1.push_back(desc[j]);
	matches2.push_back(ballTree[i]);
      }      
    } 
  }  
  cout << endl;
}

void searchExhaustiveNN(vector<FD* > desc, vector<FD* > ballTree, float dThres, int desc_dimension, vector<FD* > &matches1,vector<FD* > &matches2, vector<float > &distances) {
  
  float dist,mindist;
  float *d;
  uint mini=0;
  //dThres=dThres*dThres;
  for(uint j=0;j<desc.size();j++){
    d=desc[j]->getVec();
    //if(!(j%100))cout << "\r feature " << j << "   "<< flush;
    mindist=MAX_FLOAT;
    for(uint i=0;i<ballTree.size();i++){
      dist=distEuc(d,ballTree[i]->getVec(),desc_dimension);
      if(dist<mindist){
	mini=i;
	mindist=dist;
      }      
    } 
    if(mindist<dThres){
      distances.push_back(mindist);
      matches1.push_back(desc[j]);
      matches2.push_back(ballTree[mini]);
    }          
  }  
  cout << endl;
}


void searchExhaustivekNN(vector<FD* > desc, vector<FD* > ballTree, int k, float dThres, int desc_dimension, vector<FD* > &matches1,vector<FD* > &matches2, vector<float > &distances) {
  
  float dist;//,mindist;
  float *d;
  //uint mini=0;
  //dThres=dThres*dThres;
  
  //int *kNN = new
  vector<int > kNN;  
  vector<float > kNNdist;
  for(int i=0;i<k+1;i++){
    kNN.push_back(0);
    kNNdist.push_back(MAX_FLOAT);
  }
  
  int f;
  for(uint j=0;j<desc.size();j++){
    d=desc[j]->getVec();
    if(!(j%100))cout << "\r feature " << j << "   "<< flush;
    for(int i=0;i<k+1;i++){
      kNN[i]=0;
      kNNdist[i]=MAX_FLOAT;
    }
    
    for(uint i=0;i<ballTree.size();i++){
      dist=distEuc(d,ballTree[i]->getVec(),desc_dimension);
      if(dist<kNNdist[k]){
	f=0;
	while(dist>kNNdist[f] && f<=k){	  
	  //cout << f << " " << dist <<  " "<< kNNdist[f]<< endl; 
	  f++;
	}
	//cout << f << " " << dist <<  " "<< kNNdist[0]<<  " "<< kNNdist[1]<<  " "<< kNNdist[2]<< endl; getchar();
	if(f<k){
	  kNN.insert(kNN.begin()+f,i);
	  kNNdist.insert(kNNdist.begin()+f,dist);
	  kNN.pop_back();
	  kNNdist.pop_back();
	}
      }      
    } 
    
    for(int i=0;i<k;i++){
      if(kNNdist[i]<dThres){
	distances.push_back(kNNdist[i]);
	matches1.push_back(desc[j]);
	matches2.push_back(ballTree[kNN[i]]);
      }  
    }        
  }  
  cout << endl;
}




void searchBallTrees(FD *desc, vector<FD* > nodes, int isSwap, float dThres, int level, int desc_dimension,
		    vector<FD* > &matches1,vector<FD* > &matches2, vector<float > &distances){  
  float dist;
  //cout << " searching ball "<< endl;
  for(uint i=0;i<nodes.size();i++){
    dist=distEuc(desc->getVec(),nodes[i]->getVec(),desc_dimension);
    //cout << i << " dist "<< dist << " max " <<nodes[i]->par[RADIUS]<< " d-m "<< dist-nodes[i]->par[RADIUS]<< " " <<  dThres << " "  << nodes.size()<<  endl;getchar();
    if(dist-nodes[i]->par[RADIUS]<=dThres){
      if(nodes[i]->getTreeLevel()>level){
	searchBallTrees(desc, nodes[i]->features, isSwap,  dThres, level, desc_dimension, matches1, matches2,distances);	
      }else if(dist<dThres){
	if(isSwap){
	  distances.push_back(dist);
	  matches1.push_back(nodes[i]);
	  matches2.push_back(desc);	  
	}else{
	  distances.push_back(dist);
	  matches2.push_back(nodes[i]);
	  matches1.push_back(desc);	  	  
	}
      }
    }      
  }  
}


void searchBallTrees( vector<FD* > nodes1, vector<FD* > nodes2, float dThres, int level1, int level2,  int desc_dimension,
		     vector<FD* > &matches1,vector<FD* > &matches2, vector<float > &distances){  
  float dist;
  //cout << " searching ball "<< endl;
  for(uint j=0;j<nodes1.size();j++){
    for(uint i=0;i<nodes2.size();i++){
      dist=distEuc(nodes1[j]->getVec(),nodes2[i]->getVec(),desc_dimension);
      //{cout << i << " dist "<< dist << " max " <<nodes[i]->par[RADIUS]<< " d-m "<< dist-nodes[i]->par[RADIUS]<< " " <<  dThres << " "  << nodes.size()<<  endl;getchar();}
      if(dist-nodes1[j]->par[RADIUS]-nodes2[i]->par[RADIUS] < dThres){	
	if(nodes1[j]->getTreeLevel()>level1 && nodes2[i]->getTreeLevel()>level2){
	  searchBallTrees(nodes1[j]->features, nodes2[i]->features,  dThres, level1, level2, desc_dimension, matches1, matches2,distances);		  
	}else if(nodes1[j]->getTreeLevel()==level1 && nodes2[i]->getTreeLevel()>level2){	  
	  searchBallTrees(nodes1[j], nodes2[i]->features, 0, dThres, level2, desc_dimension, matches1, matches2,distances);		  
	}else if(nodes1[j]->getTreeLevel()>level1 && nodes2[i]->getTreeLevel()==level2){	  
	  searchBallTrees(nodes2[i], nodes1[j]->features, 1,  dThres, level1, desc_dimension, matches1, matches2,distances);	  
	}else if(dist<dThres){
	  distances.push_back(dist);
	  matches1.push_back(nodes1[j]);
	  matches2.push_back(nodes2[i]);	  
	}
      }      
    }  
  }
}
