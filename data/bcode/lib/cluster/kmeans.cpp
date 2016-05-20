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

void updateMeanClusters(vector<FD *> features, vector<FD *> &clusters, uint *assign, int *changed){
  
  int dim=features[0]->getSize();
  float nbf=0;
  float *mvec=new float[dim];
  float *vec,*cc;

  for(uint i=0;i<clusters.size();i++){
    if(changed[i]){
      bzero(mvec,dim*sizeof(float));
      nbf=0;cc=clusters[i]->getVec();
      for(uint j=0;j<features.size();j++){
	if(assign[j]==i){
	  vec=features[j]->getVec();
	  nbf++;
	  for(int d=0;d<dim;d++){
	    mvec[d]+=vec[d];
	  }
	}
      }
      for(int d=0;d<dim;d++){
	cc[d]=mvec[d]/nbf;
      }  
      clusters[i]->setNbf((int)nbf);
    }
  }
  delete mvec;
}


/*greedy initialization*/
void initClusters(vector<FD *> features, vector<FD *> &clusters, int K, uint *assign, int *changed)
{
  //cout << "Furthest first initialization..." << flush;
  /*put the first center to the mean*/
  if(features.size()<5){cout <<  "Too few feature in kmeans: "<<  features.size()<< endl;exit(0);}
  FD *clust;
  clust=new FD();
  int dim=features[0]->getSize();
  clust->allocVec(dim);
  clusters.push_back(clust);  
  updateMeanClusters(features, clusters, assign, changed);  

  DARY *mindist=new DARY((int)1,(int)features.size(),FLOAT1, MAX_FLOAT);
  DARY *mincenter=new DARY(1,features.size(),FLOAT1,0.0);
  DARY *inplay=new DARY(1,features.size(),FLOAT1,1.0);
  DARY *centerdist=new DARY(K,K,FLOAT1,0.0);
  float dist, maxdist;
  int maxi;

  FD *lastclust;
  for(int k=0;k<K;k++){    
    
    /*add new cluster center as a furthest assigned point*/
    maxdist=0;
    maxi=-1;
    if(k>0){
      for(uint i=0;i<features.size();i++){
	if(mindist->fel[0][i]>maxdist){
	  maxdist=mindist->fel[0][i];
	  maxi=i;
	}      
      }
      if(maxi<0){cout << "error in greedy"<< endl;exit(0);}
      clust=new FD();
      clust->copy(features[maxi]);
      clusters.push_back(clust);  
    }
    lastclust=clusters[k];


    /*update distances between cluster centers*/
    for(int i=0;i<k;i++){      
      centerdist->fel[k][i]=0.5*distEucSqrt(clusters[i]->getVec(),lastclust->getVec(),dim);
      centerdist->fel[i][k]=centerdist->fel[k][i];
    }
    
    /*find valid points using triangular inequality compared to between class distances*/
    for(uint i=0;i<features.size() && k>0;i++){
      if(mindist->fel[0][i] > (centerdist->fel[k][(int)mincenter->fel[0][i]])){
	inplay->fel[0][i]=1;
      }else inplay->fel[0][i]=0;
    }

    /*reassign points to the new cluster*/
    for(uint i=0;i<features.size();i++){      
      if((int)inplay->fel[0][i]){
	dist=distEucSqrt(features[i]->getVec(),lastclust->getVec(),dim);
	if(mindist->fel[0][i]>dist){
	  mindist->fel[0][i]=dist; 
	  mincenter->fel[0][i]=k;
	}  
      }
    }   
  }
  delete mindist;
  delete mincenter;
  delete inplay;
  delete centerdist;
  //cout << "done" << endl;  
}


void initRandomClusters(vector<FD *> features, vector<FD *> &clusters, int K, uint *assign, int *changed)
{
  /* set the cluster centers to randomly selected points */
  if(features.size()<K){cout <<  "Too few feature in kmeans: "<<  features.size()<< endl;exit(0);}
  FD *clust;
  clust=new FD();
  int dim=features[0]->getSize();
  /* initialize random seed: */
  srand ( time(NULL) );

  uint idx = (uint) ((features.size()-1)*(rand()/((float)RAND_MAX)));  
  clust->copy(features[idx]);
  clusters.push_back(clust);  
  updateMeanClusters(features, clusters, assign, changed);  

  int flag,dist;
  uchar *fl = new uchar[features.size()];
  bzero(fl,sizeof(uchar)*features.size());
  int rounds=0;
  for ( int i=1; i < K ; i++ ) {
    if(!(clusters.size()%10000))cout << "\rClusters done: "<< clusters.size() << " of "<< K << " dist " << dist << " rounds " << rounds <<  "    "<< flush;
    flag=1;
    rounds++;
    idx = (uint)((features.size()-1)*(rand()/((float)RAND_MAX)));
    if(fl[idx])flag=0;
    float *vec = features[idx]->getVec();
    for(uint c=0;c<clusters.size() && flag;c++){
       dist = distEuc(clusters[c]->getVec(),vec,dim);
      //cout << "dist " << dist << endl;
      if(dist<0.01)flag=0;//make sure there is no similar cluster already
    }
    if(flag){
      clust=new FD();
      clust->copy(features[idx]); 
      clusters.push_back(clust); 
      fl[idx]=1;
      rounds=0;
    }else if(rounds<features.size())i--;
    else{
     cout << "Can't find " << K << " different features, found " << clusters.size() <<   endl;   
    }
    //cout << idx << "  ";
    //features[idx]->Cout();
    //clusters[i]->Cout();
  }  
  delete []fl;
  // cout << endl;
}

void randomClusters(vector<FD *> features, vector<FD *> &clusters, int K){
  uint *assign=new uint[features.size()];//table for cluster index
  bzero(assign,features.size()*sizeof(uint));
  int *changed=new int[K];
  memset(changed,1,K*sizeof(int));

  //produces  very unbalanced clusters since based on outliers, many partitions have size 1
  //initClusters(features, clusters,  K, assign, changed);
  //random is much better but results in more iterations
  initRandomClusters(features, clusters,  K, assign, changed);
  delete []assign;
  delete []changed; 
}

inline float kmeansEuc(float* vec1,float* vec2, int size1, int size2, int size, float thres){
  float dist=0,d;    
  for(int i=0;i<size1;i++){
    d=vec1[i]-vec2[i];
    dist+=(d*d);
  }
  if(dist>thres)return dist+thres;
  for(int i=size1;i<size2;i++){
    d=vec1[i]-vec2[i];
    dist+=(d*d);
  }
  if(dist>thres)return dist+thres;
  for(int i=size2;i<size;i++){
    d=vec1[i]-vec2[i];
    dist+=(d*d);
  }
  return dist;
}

void updateBetweenCustersDist(vector<FD *> clusters, int *changed, DARY *cl_dist){
  int dim=clusters[0]->getSize();
  for(uint i=0;i<cl_dist->y();i++){
    if(changed[i]){
      for(uint j=i+1;j<cl_dist->x();j++){
	cl_dist->fel[i][j]=0.5*distEucSqrt(clusters[i]->getVec(),clusters[j]->getVec(),dim);
	cl_dist->fel[j][i]=cl_dist->fel[i][j];
      }
    }
  }
}

int assignFeatures(vector<FD *> features, vector<FD *> &clusters, uint *assign, float *fassign, int *changed){
  /****************************************************/
  /* assign every point to the nearest cluster center */
  /****************************************************/
  int dim=features[0]->getSize();
  float min_dist,dist;
  int   min_idx=-1;
  int nchanged=0;
  int size1=(int)(0.3*dim);
  int size2=(int)(0.6*dim);
  
  bzero(changed,clusters.size()*sizeof(int));
  
  for ( uint i=0; i < features.size(); i++ ) { 
    float *vec=features[i]->getVec();
    /* calculate the min distance to every cluster center */
    min_dist = MAX_FLOAT;min_idx=-1;
    
    for ( uint j=0; j < clusters.size(); j++ ) {
      dist = kmeansEuc(clusters[j]->getVec(),vec,size1,size2,dim,min_dist);
      if ( dist < min_dist ) {
	min_dist = dist;
	min_idx  = j;
      }
    }       
    /* update the cluster assignment */

    /* if(min_idx==-1){cout << "ERROR: assignFeatures "<< dist << endl;exit(0);
      features[i]->Cout(features[i]->getSize());
      min_dist = MAX_FLOAT;
      for ( uint k=0; k < clusters.size(); k++ ) {
	dist = kmeansEuc(clusters[k]->getVec(),vec,size1,size2,dim,min_dist);
	cout << dist << endl;
	clusters[k]->Cout(clusters[k]->getSize());
      }
      getchar();}
    */

    if((int)assign[i]!=min_idx){
      changed[assign[i]]=1;
      assign[i]=min_idx;
      changed[min_idx]=1;
      nchanged++;
    }
    fassign[i]=min_dist;
  }
  /*
  int flag=1;
  for ( uint j=0; j < clusters.size(); j++ ) {
    flag=1;
    for ( uint n=0; n < features.size(); n++ ) { 
      if(j==assign[n])flag=0;
    }
    if(flag==1){
      cout << "nothing assigned to "<< j << endl;
      cout << j;clusters[j]->Cout();
      for ( uint k=0; k < features.size(); k++ ) { 
	if(features[k]->getX()==clusters[j]->getX()){
	  cout << k;features[k]->Cout();
	  cout << "dist "<< distEuc(clusters[j]->getVec(),features[k]->getVec(),dim)<< endl;
	  cout << "fassign[i] "<< fassign[k]<< " " << assign[k]<< endl;
	  clusters[assign[k]]->Cout();getchar();
	}
      }
    }
  }
  */
  updateMeanClusters(features, clusters, assign, changed);

  return nchanged;
} 
         
void kmeans(vector<FD *> features, vector<FD *> &clusters, vector<ClStep *> &trace, int K, int max_iter,  int verbal){ 
  if(verbal)cout << "K-means with : K=" << K << ", and "<<  max_iter<< " iterations "<< endl;
  uint *assign=new uint[features.size()];//table for cluster index
  bzero(assign,features.size()*sizeof(uint));
  float *fassign=new float[features.size()];
  bzero(fassign,features.size()*sizeof(float));
  int *changed=new int[K];
  memset(changed,1,K*sizeof(int));

  //produces  very unbalanced clusters since based on outliers, many partitions have size 1
  //initClusters(features, clusters,  K, assign, changed);
  //random is much better but results in more iterations
  initRandomClusters(features, clusters,  K, assign, changed);
  int  nchanged;
  int iter=0;
do{
    nchanged=assignFeatures(features, clusters, assign, fassign, changed);//nchanged=0;
  }while(nchanged && iter++<max_iter);

  int *cl_index=new int[clusters.size()];
  //for(uint i=0;i<clusters.size();i++)cl_index[i]=-1;
  memset(cl_index,(-1),clusters.size()*sizeof(int));

  for(uint i=0;i<features.size();i++){
    //cout << "assign  "<< i << " " << assign[i] << endl;
    clusters[assign[i]]->features.push_back(features[i]);
    if(cl_index[assign[i]]==-1){
      cl_index[assign[i]]=features[i]->getTreeLevel();
      
      //can't be assigned to itself otherwise trace2clusters problems
    }else trace.push_back(new ClStep(fassign[i],cl_index[assign[i]],features[i]->getTreeLevel()));
  }
  for(uint i=0;i<clusters.size();i++){
    if(cl_index[i]!=-1){
      clusters[i]->setTreeLevel(cl_index[i]);
    }else {
      cout << "ERROR: Problem in kmeans - zero size cluster"<< endl;
    }
  }
  /*
    int nb=0;
  for(uint i=0;i<clusters.size();i++){//remove one member partitions
    if(clusters[i]->features.size()<2){
    nb++;cout << nb<< endl;
      delete clusters[i];
      clusters.erase(clusters.begin()+i);i--;
    }
    }
  */

  //verbal=1;
  if(verbal)cout << "iter " << iter<< ": error " << nchanged <<" trace.size() " << trace.size()<<" top_sim " << trace[trace.size()-1]->sim<<  endl;
  //getchar();
  delete []assign;
  delete []fassign;
  delete []cl_index;
  delete []changed;
  
}

void kmeans(vector<FD *> &features,  int K, int max_iter,  int verbal){ 
  if(verbal)cout << "K-means with : K=" << K << ", and "<<  max_iter<< " iterations, feats: "<< features.size()<<  endl;
  uint *assign=new uint[features.size()];
  bzero(assign,features.size()*sizeof(uint));
  float *fassign=new float[features.size()];
  bzero(fassign,features.size()*sizeof(float));
  
  
  int *changed=new int[K];
  memset(changed,1,K*sizeof(int));

  vector<FD *> clusters;
  initRandomClusters(features, clusters,  K,assign, changed);
 //produces  very unbalanced clusters since based on outliers, many partitions have size 1
  //initClusters(features, clusters,  K, assign, changed);
  //random is much better but results in more iterations
  
 
  int  nchanged;
  int iter=0;
  do{
    nchanged=assignFeatures(features, clusters, assign, fassign, changed);//nchanged=0;
  }while(nchanged && iter++<max_iter);

  
  for(uint i=0;i<features.size();i++){
    clusters[assign[i]]->features.push_back(features[i]);
  }
 if(verbal)cout << "iter " << iter<< ": error " << nchanged << " clusters " << clusters.size() << " of " << features.size()<<  endl;

  features=clusters;

   //getchar();
  delete []assign;
  delete []fassign;
  delete []changed;
   
}


void kmeans(vector<FD *> features, vector<FD *> &clusters, int K, int max_iter, int verbal){
  vector<ClStep *> trace;
  for(uint i=0;i<features.size();i++)
    features[i]->setTreeLevel(i);//this numbers are used to do the clustering trace

  kmeans(features, clusters, trace,  K, max_iter, verbal);
  trace.clear(); 
  for(uint i=0;i<features.size();i++)
    features[i]->setTreeLevel(0);
}


void hkmeans(vector<FD *> &features, vector<FD *> &clusters, int K, uint tot_K, int min_cluster, int max_iter, int level, int verbal){
    
  
  level++;
  kmeans(features,   K,  max_iter,   verbal);
  
  for(uint i=0;i<features.size();i++){
      if(features[i]->features.size()>min_cluster && features.size()>1){
	hkmeans(features[i]->features, clusters, K, tot_K, min_cluster, max_iter, level, verbal);
      } else {//if(clusters.size()<tot_K){
	cout << "\rhkmeans " << level << " tot_clust  " << clusters.size() <<"  members "<< features[i]->features.size()<<  "    "<< flush;
	if(features[i]->features.size())clusters.insert(clusters.end(),features[i]->features[0]);
      }      
  }
   
}


void kmeans(vector<FD *> features, vector<ClStep *> &trace, int K, int max_iter, int verbal){
  vector<FD *> clusters;

  for(uint i=0;i<features.size();i++)
    features[i]->setTreeLevel(i);//this numbers are used to do the clustering trace
  
  kmeans(features, clusters, trace,  K, max_iter, verbal);
  
  for(uint i=0;i<clusters.size();i++){
    clusters[i]->features.clear();
    delete clusters[i];
  }
  clusters.clear();
  for(uint i=0;i<features.size();i++)
    features[i]->setTreeLevel(0);

}

void getRandomizedSets(vector<FD *> &features, int K, int verbal){
  if(verbal)cout << "Random with : K=" << K << endl;
  
  FD *clust;
  int dim=features[0]->getSize();   
  vector<FD *> clusters;
  
  for(uint i=0;i<K;i++){
    clust=new FD();
    clust->allocVec(dim);
    clusters.push_back(clust);
  }
  
  uint c=0;
  uint fsize=features.size();
  for(uint f=0;f<fsize;f++){   
      clusters[c]->features.push_back(features[f]);
      c++;
      if(c==K)c=0;      
  }        

  features=clusters;

}
