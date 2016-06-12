#ifndef _clustering_h_
#define _clustering_h_

#include "../descriptor/feature.h"

#define AVLINK 1 
#define KMEANS 2

class Tree{
 public:
  Tree(){}
  vector<FD *> features;
  uint type;
};


class ClStep{
 public:
  int cf1,f2;
  float sim; 
  ClStep(float _sim, int _cf1, int _f2){
    sim=_sim;cf1=_cf1;f2=_f2;
  }
  ~ClStep(){}
};


void getClusterMembers(vector<FD *> inClusters, vector<FD *> &outDesc, int level);
void sortTrace(vector<ClStep *> &trace);
void updateMeanVar(vector<FD *> &vClusters, float quant=1);
void updateMeanVar(FD * cluster, float quant=1);

void getRandomizedSets(vector<FD *> &features, int K,  int verbal=1);
void kmeans(vector<FD *> &features,  int K, int max_iter,  int verbal=1);
void kmeans(vector<FD *> features, vector<FD *> &clusters, int K, int max_iter=35, int verbal=1);
void kmeans(vector<FD *> features, vector<ClStep *> &trace, int K, int max_iter=35, int verbal=1);
void kmeans(vector<FD *> features, vector<FD *> &clusters, vector<ClStep *> &trace, int K, int max_iter, int verbal=1);
void kmeansTree(vector<FD *> &features,vector<ClStep *> &trace, uint min_size, int K, int max_iter);
void kmeansTreeB(vector<FD *> &features, vector<ClStep *> &trace, uint min_size, int K, int max_iter);
void kmeansTree(vector<FD *> &features, uint min_size, int K, int max_iter, float quantile);
void randomClusters(vector<FD *> features, vector<FD *> &clusters, int K);
void hkmeans(vector<FD *> &features, vector<FD *> &clusters, int K, uint tot_K, int min_cluster, int max_iter, int level, int verbal);



void avlink(vector<FD *> features, vector<ClStep *> &trace, float sim); 
void kmeansAgglo(vector<FD *> ivDesc, vector<FD *> &clusters, vector<ClStep *> &trace,  uint K, float dSim);


void traceToBallTree( vector<FD *> ivDesc, vector<ClStep *> trace,
		      vector<FD *> &vClusters, float minThresh, 
		      uint memNb, float quant);
void traceToBallTree( vector<FD *> ivDesc, vector<ClStep *> trace,
		      vector<FD *> &vClusters, float minThresh, uint levels, float maxThresh, 
		      float quant, uint min_cluster_size);

void traceToBallTree( vector<FD *> ivDesc, vector<ClStep *> trace,
		      vector<FD *> &vClusters, float minThresh,  float quant, uint min_cluster_size);

void trace2ClustersSim( vector<FD *> ivDesc, vector<ClStep *> trace,
		     vector<FD *> &vDesc, float dSimThresh);

float trace2ClustersComp( vector<FD *> ivDesc, vector<ClStep *> trace,
		        vector<FD *> &vDesc, uint cluster_nb);
void traceToClusters( vector<FD *> ivDesc, vector<ClStep *> trace,
		      vector<FD *> &vClusters, float dSimThresh );

void saveBallTreeBinary(vector<FD *> tree,  uint fnb, const char *filename);
void loadBallTreeBinary(const char *filename, vector<FD *> &tree, uint &fnb);
void saveBallTree(vector<FD *> clusters,  uint fnb, const char *filename);
void loadBallTree(const char *filename, vector<FD *> &clusters, uint &fnb);


void loadTrace(const char *filename, vector<ClStep *> &vClusterTrace);
void saveTrace(vector<ClStep *> vClusterTrace, const char *filename);


void searchBallTree(FD*  desc, vector<FD* > nodes, float dThres, int level, int desc_dimension, vector<FD* > &matches, vector<float > &distances);  
void searchBallTreeNN(vector<FD* > desc, vector<FD* > ballTree, float dThres, int level, int desc_dimension, vector<FD* > &match1, vector<FD* > &match2, vector<float > &distances);



void searchBallTree(FD *desc, vector<FD* > nodes, float dThres, int level, int desc_dimension, vector<FD* > &matches1, vector<FD* > &matches2, vector<float > &distances); 
void searchBallTree(vector<FD* > desc, vector<FD* > ballTree, float dThres, int level, int desc_dimension, vector<FD* > &matches1, vector<FD* > &matches2, vector<float > &distances);
void searchBallTreeNN(FD *desc, vector<FD* > nodes, float &dThres, int level, int desc_dimension, vector<FD* > &match, vector<float > &distances);
void searchBallTreeNister(FD *desc, vector<FD* > nodes, float dThres, int level, int desc_dimension, vector<FD* > &matches, vector<float > &distances);
void searchBallTrees(FD *desc, vector<FD* > nodes, int isSwap, float dThres, int level, int desc_dimension,
		    vector<FD* > &matches1,vector<FD* > &matches2, vector<float > &distances);
void searchBallTrees( vector<FD* > nodes1, vector<FD* > nodes2, float dThres, int level1, int level2,  int desc_dimension,
		     vector<FD* > &matches1,vector<FD* > &matches2, vector<float > &distances);  
void searchBallTreeNister(vector<FD* > desc, vector<FD* > ballTree, float dThres, int level, int desc_dimension, vector<FD* > &matches, vector<float > &distances);

void searchExhaustive(vector<FD* > desc, vector<FD* > ballTree, float dThres, int desc_dimension, vector<FD* > &matches1, vector<FD* > &matches2, vector<float > &distances) ;void searchBallTreeNNT(vector<FD* > desc, vector<FD* > ballTree, float dThres, int level, int desc_dimension, vector<FD* > &matches, vector<float > &distances); 

void searchBallTreeNN(vector<FD* > desc, vector<FD* > ballTree, float dThres, int level, int desc_dimension, vector<FD* > &match1, vector<FD* > &match2, vector<float > &distances);  

void searchExhaustiveNN(vector<FD* > desc, vector<FD* > ballTree, float dThres, int desc_dimension, vector<FD* > &matches1,vector<FD* > &matches2, vector<float > &distances);

void searchExhaustivekNN(vector<FD* > desc, vector<FD* > ballTree, int k, float dThres, int desc_dimension, vector<FD* > &matches1,vector<FD* > &matches2, vector<float > &distances);

void searchBallTreekNN(vector<FD* > desc, vector<FD* > ballTree, float dThres, int k, int level, int desc_dimension, vector<FD* > &match1, vector<FD* > &match2,vector<float > &distances);

#endif

