/*********************************************************************/
/*                                                                   */
/* FILE         clrnncagglo.cc                                       */
/* AUTHORS      Bastian Leibe                                        */
/* EMAIL        leibe@inf.ethz.ch                                    */
/*                                                                   */
/* CONTENT      Implements an average-link agglomerative clustering  */
/*              method based on the Reciprocal Nearest Neighbor      */
/*              (RNN) Chain strategy by C. de Rham and P. Benzecri.  */
/*              Details on this algorithm can be found in            */
/*                                                                   */
/*                C. de Rham,                                        */
/*                La Classification Hierarchique Ascendante Selon la */
/*                Methode des Voisins Reciproques.                   */
/*                Les Cahiers de l'Analyse des Donnees, Vol. 5,      */
/*                pp. 135-144, 1980.                                 */
/*                                                                   */
/*                J.P. Benzecri,                                     */
/*                Construction d'une Classification Ascendante Hier- */
/*                archique par la Recherche en Chaine des Voisins    */
/*                Reciproques.                                       */
/*                Les Cahiers de l'Analyse des Donnees, Vol. 7,      */
/*                pp. 209-218, 1982.                                 */
/*                                                                   */
/*              This file contains an improvement to the basic RNNC  */
/*              algorithm that does not need to compute and store an */
/*              explicit O(N^2) distance matrix. As a result, it is  */
/*              extremely efficient and has a worst-case runtime     */
/*              performance of O(N^2) while needing only O(N) space. */
/*              (A further improvement to an expected-time perfor-   */
/*              mance of only O(NlogN) can be easily achieved by in- */
/*              tegrating a fast Nearest-Neighbor search strategy).  */
/*                                                                   */
/*              The algorithm is written for average-link clustering */
/*              with correlation or Euclidean distance as a distance */
/*              measure, but it can be easily extended to any dis-   */
/*              tance measure where the group average criterion      */
/*               d(X,Y)=(1/|X||Y|)sum_{xi in X}sum_{yj in Y}d(xi,yj) */
/*              can be reformulated to only depend on the cluster    */
/*              means (and possibly the variances).                  */
/*                                                                   */
/*                                                                   */
/* BEGIN        Fri Jun 04 2004                                      */
/* LAST CHANGE  Fri Jun 04 2004                                      */
/*                                                                   */
/*********************************************************************/

/****************/
/*   Includes   */
/****************/
#include <mex.h>
#include <iostream>
#include <iomanip>

#include <math.h>
#include <stdlib.h>

#include "clrnncagglo.hh"

/****************/
/* Definitions  */
/****************/

const float MIN_SIM = -99999999999.9;




/*===================================================================*/
/*                         Class ClPoint                             */
/*===================================================================*/

/***********************************************************/
/*                      Constructors                       */
/***********************************************************/

ClPoint::ClPoint()
  /* standard constructor */
{
  m_data.clear();
}


ClPoint::ClPoint( int dim, FLOAT value )
  /* alternate constructor */
{
  m_data.resize( dim );
  for( int i=0; i<dim; i++ )
    m_data[i] = value;
}


ClPoint::ClPoint( vector<FLOAT> data )
  /* alternate constructor */
{
  m_data = data;
}


ClPoint::ClPoint( const ClPoint &other )
  /* copy constructor */
{
  copyFromOther( other );
}


ClPoint& ClPoint::operator=( ClPoint other )
  /* assignment operator */
{
  copyFromOther( other );
  return *this;
}


ClPoint& ClPoint::operator=( vector<float> data )
  /* assignment operator */
{
  m_data.resize( data.size() );
  for( int i=0; i<data.size(); i++ )
    m_data[i] = (FLOAT) data[i];
  return *this;
}


ClPoint& ClPoint::operator=( vector<double> data )
  /* assignment operator */
{
  m_data.resize( data.size() );
  for( int i=0; i<data.size(); i++ )
    m_data[i] = (FLOAT) data[i];
  return *this;
}


ClPoint::~ClPoint()
  /* standard destructor */
{
}


void ClPoint::copyFromOther( const ClPoint &other )
  /* Auxiliary function to copy from another ClPoint.                */
{
  m_data = other.m_data;
}


/***********************************************************/
/*                 Manipulation Functions                  */
/***********************************************************/

ClPoint& ClPoint::add( const ClPoint &other )
  /* add the value of another ClPoint */
{
  assert( dim() == other.dim() );

  for( int i=0; i<dim(); i++ )
    m_data[i] += other.m_data[i];
  return *this;
}


ClPoint& ClPoint::sub( const ClPoint &other )
  /* subtract the value of another ClPoint */
{
  assert( dim() == other.dim() );

  for( int i=0; i<dim(); i++ )
    m_data[i] -= other.m_data[i];
  return *this;
}


ClPoint& ClPoint::mul( const ClPoint &other )
  /* multiply the value of another ClPoint */
{
  assert( dim() == other.dim() );

  for( int i=0; i<dim(); i++ )
    m_data[i] *= other.m_data[i];
  return *this;
}


ClPoint& ClPoint::mul( const FLOAT c )
  /* multiply with a constant */
{
  for( int i=0; i<dim(); i++ )
    m_data[i] *= c;
  return *this;
}


ClPoint& ClPoint::div( const ClPoint &other )
  /* divide by the value of another ClPoint */
{
  assert( dim() == other.dim() );

  for( int i=0; i<dim(); i++ )
    m_data[i] /= other.m_data[i];
  return *this;
}


ClPoint& ClPoint::div( const FLOAT c )
  /* divide by a constant */
{
  FLOAT d = 1.0/c;
  for( int i=0; i<dim(); i++ )
    m_data[i] *= d;
  return *this;
}


FLOAT ClPoint::distTo( const ClPoint &other ) const
  /* compute the distance to another ClPoint */
{ return sqrt( distSqrTo( other ) ); }


FLOAT ClPoint::distSqrTo( const ClPoint &other ) const
  /* compute the squared distance to another ClPoint */
{
  assert( dim() == other.dim() );
  
  FLOAT sum = 0.0;
  for (int i=0; i<dim(); i++) {
    FLOAT dist = m_data[i] - other.m_data[i];
    sum += dist*dist;
  }
  return sum;
}


FLOAT ClPoint::scalarProduct( const ClPoint &other ) const
  /* compute the scalar product with another ClPoint */
{
  assert( dim() == other.dim() );
  
  FLOAT sum = 0.0;
  for (int i=0; i<dim(); i++) {
    sum += m_data[i] *other.m_data[i];
  }
  return sum;
}


/***********************************************************/
/*                  Associated Functions                   */
/***********************************************************/

ostream& operator<<( ostream& os, const ClPoint& pt )
  /* print operator */
{
  os << "(";
  for( int i=0; i<pt.dim()-1; i++ )
    os << setw(5) << setprecision(4) << pt.m_data[i] << ",";
  os << setw(5) << setprecision(4) << pt.m_data[pt.dim()-1] << ")";

  return os;
}


/*===================================================================*/
/*                         Class Cluster                             */
/*===================================================================*/

/***********************************************************/
/*                      Constructors                       */
/***********************************************************/

Cluster::Cluster()
  /* standard constructor */
{
  m_nNumPoints = 0;
  m_nDim = 0;
  m_vPoints.clear();

  init();
}


Cluster::Cluster( vector<ClPoint> &vPoints )
  /* alternate constructor */
{
  assert( vPoints.size() >= 1 );

  m_nNumPoints = vPoints.size();
  m_nDim = vPoints[0].dim();
  /* check if all the points have the same dimensionality */
  for( int i=1; i<m_nNumPoints; i++ ) {
    assert( vPoints[i].dim() == m_nDim );
  }

  m_vPoints = vPoints;

  init();
}


Cluster::~Cluster()
  /* standard destructor */
{
  m_vPoints.clear();
  m_vCenters.clear();
}


void Cluster::init()
{
  m_nNumClusters = 0;
  m_vCenters.clear();
  m_fError = 0.0;
}


/***********************************************************/
/*                    Access Functions                     */
/***********************************************************/

void Cluster::addPoint( ClPoint &pt )
  /* add a point to the list of points to be clustered. */
{
  if( m_nNumPoints == 0 )
    /* add the first point */
    m_nDim = pt.dim();
  else {
    /* check for consistency */
    assert( pt.dim() == m_nDim );
  }

  m_vPoints.push_back( pt );
  m_nNumPoints++;
}


/***********************************************************/
/*                   Clustering Functions                  */
/***********************************************************/

void Cluster::initClusters( int nNumClusters )
  /* initialize the clustering algorithm with a cluster center in every feature */
{
  cout << "    Cluster::initClusters(" << nNumClusters << ") called."
       << endl;

  m_nNumClusters = nNumClusters;;

  cout << "      Resizing the data vectors..." << endl;
  initDataVectors();

  /**********************************/
  /* Initialize the cluster centers */
  /**********************************/
  cout << "      Initializing the cluster centers..." << endl;
  initClusterCenters();

  /*******************************/
  /* Initialize the sigma values */
  /*******************************/
  cout << "      Initializing the sigma values..." << endl;
  initCovariances();

  /*********************************/
  /* Initialize the cluster priors */
  /*********************************/
  cout << "      Initializing the cluster priors..." << endl;
  initPriors();
  
  /******************************************/
  /* Initialize the posterior probabilities */
  /******************************************/
  cout << "      Initializing the posterior probabilities..." << endl;
  initPosteriors();
  
  /***********************/
  /* Calculate the error */
  /***********************/
  cout << "      Calculating the error..." << endl;
  m_fError = calculateError();

  cout << "    done." << endl;
}


void Cluster::initClusters( vector<ClPoint> vCenters )
  /* initialize the clustering algorithm with given cluster centers. */
{
  cout << "    Cluster::initClusters(" << vCenters.size() << ") called"
       << " with given cluster centers." << endl;

  m_nNumClusters = vCenters.size();

  cout << "      Resizing the data vectors..." << endl;
  initDataVectors();

  /**********************************/
  /* Initialize the cluster centers */
  /**********************************/
  cout << "      Copying the cluster centers..." << endl;
  m_vCenters = vCenters;

  /*******************************/
  /* Initialize the sigma values */
  /*******************************/
  cout << "      Initializing the sigma values..." << endl;
  initCovariances();

  /*********************************/
  /* Initialize the cluster priors */
  /*********************************/
  cout << "      Initializing the cluster priors..." << endl;
  initPriors();
  
  /******************************************/
  /* Initialize the posterior probabilities */
  /******************************************/
  cout << "      Initializing the posterior probabilities..." << endl;
  initPosteriors();
  
  /***********************/
  /* Calculate the error */
  /***********************/
  cout << "      Calculating the error..." << endl;
  m_fError = calculateError();

  cout << "    done." << endl;
}


bool Cluster::doClusterSteps( FLOAT eps, int max_iter )
  /* Iterate the clustering algorithm until either the error doesn't */
  /* decrease by more than eps anymore, or a maximum number of ite-  */
  /* rations is reached. Returns true if the convergence limit is    */
  /* reached, else false.                                            */
{
  int   iter = 0;
  FLOAT old_error = 999999999.9;
  
  cout << "    Cluster::doClusterSteps(" << eps << "," << max_iter 
       << ") called." << endl;

  /************************************/
  /* Iterate the clustering procedure */
  /************************************/
  iter = 0;
  while ( (iter < max_iter) && ((iter<1) || (old_error - m_fError > eps)) ) {
    iter++;
    old_error = m_fError;

    /****************************************************/
    /* assign every point to the nearest cluster center */
    /****************************************************/
    doReestimationStep();
    
    /******************************/
    /* update the cluster centers */
    /******************************/
    doUpdateStep();

    /*************************/
    /* recalculate the error */
    /*************************/
    m_fError = calculateError();

    cout << "      Iteration " << setw(3) << iter 
	 << ": error = " << setprecision(6) << m_fError 
	 << ", difference = " << setprecision(6)
	 << old_error - m_fError << endl;
  } /* end while */

  if (old_error - m_fError > eps)
    return false;
  else
    return true;
}


/*===================================================================*/
/*                         Class ClFastAgglo                         */
/*===================================================================*/

/***********************************************************/
/*                      Constructors                       */
/***********************************************************/

ClRNNCAgglo::ClRNNCAgglo( vector<ClPoint> &vPoints )
  : Cluster( vPoints )
{
  //m_nMetricType = METRIC_NGC;
}


/***********************************************************/
/*                     Initialization                      */
/***********************************************************/

void ClRNNCAgglo::initDataVectors()
{
  m_vCenters.resize( m_nNumClusters );
  m_vBelongsTo.resize( m_nNumPoints );
  m_vVariances.resize( m_nNumClusters, 0.0 );
  
  m_vTrace.clear();
}


void ClRNNCAgglo::initClusterCenters()
{
  /* set the cluster centers initially (one to each point) */
  for(int i=0; i < m_nNumClusters; i++ ) {
    m_vCenters[i]       = m_vPoints[i];
    m_vVariances[i]     = 0.0;
    m_vBelongsTo[i]     = i;
    m_vValid.push_back( true );

    /* initialize all Feature Vectors to each corresponding cluster center */
    vector<int> tempFV;
    tempFV.push_back( i );
    m_vvAllBelongingFV.push_back( tempFV );
  }
  
  /* free the memory used for points */
  m_vPoints.clear();

  cout << "    ClRNNCAgglo::initClusterCenters() done." << endl;      
}


void ClRNNCAgglo::initCovariances()
{}


void ClRNNCAgglo::initPriors()
{}


void ClRNNCAgglo::initPosteriors()
{
  //doReestimationStep();
}


float ClRNNCAgglo::sim( int idx1, int idx2 )
{
  float sim = MIN_SIM;
  ClPoint &p1 = m_vCenters[idx1];
  ClPoint &p2 = m_vCenters[idx2];

  /*  switch( m_nMetricType ) {
  case METRIC_NGC:
    sim = p1.scalarProduct( p2 );
    break;
 
  case METRIC_EUCLID:
  */


  sim = -(m_vVariances[idx1] + m_vVariances[idx2] + p1.distSqrTo( p2 ));
 

   /*  break;
    
  default:
      cerr << "Error in ClRNNCAgglo::sim(): "
           << "Unknown similarity measure (" << m_nMetricType << ")" << endl;
  }
    */
  return sim;
}


void ClRNNCAgglo::agglomerate( int idx1, int idx2, int newidx )
  /*******************************************************************/
  /* Merge 2 neighboring clusters if compactness is guaranteed.      */
  /*******************************************************************/
{
  /*--------------------------------*/
  /* Compute the new cluster center */
  /*--------------------------------*/
  /* The combined cluster center can be computed as follows:       */
  /*   c_new = (1/N+M)*(N*c_x + M*c_y)                             */
  ClPoint aCenter, bCenter, newCenter;
  float   N, M;
  aCenter = m_vCenters[idx1];
  bCenter = m_vCenters[idx2];
  N       = (float) m_vvAllBelongingFV[idx1].size();
  M       = (float) m_vvAllBelongingFV[idx2].size();
  aCenter.mul( N );
  bCenter.mul( M );
  newCenter = aCenter;
  newCenter.add( bCenter );
  newCenter.div( N + M );
  
  
  /*--------------------------*/
  /* Compute the new variance */
  /*--------------------------*/
  /* The new variance can be computed as                           */
  /*   sig_new^2 = (1/N+M)*(N*sig_x^2 + M*sig_y^2 +                */
  /*                        NM/(N+M)*(mu_x - mu_y)^2)              */
  double aVar, bVar, newVar;
  aVar = m_vVariances[idx1];          
  bVar = m_vVariances[idx2];
  newVar = ( ( N*aVar + M*bVar + 
               N*M/(N+M)*m_vCenters[idx1].distSqrTo( m_vCenters[idx2]) ) / 
             (N + M) );
  

  /*------------------------------------------*/
  /* Add the new center to the representation */
  /*------------------------------------------*/
  /* 2 cluster centers will be invalid */
  m_vValid[idx1] = false;
  m_vValid[idx2] = false;

  /* add the new cluster center and variance */
  m_vCenters[newidx]   = newCenter;
  m_vVariances[newidx] = newVar;
  m_vValid[newidx]     = true;

  /* adjust the corresponding Feature Vectors to the new Center */
  if( newidx != idx1 )
    for( int i=0; i<(int)m_vvAllBelongingFV[idx1].size(); i++ ) {
      //new cluster assignment
      m_vBelongsTo[ m_vvAllBelongingFV[idx1][i] ] = newidx; 
    }
  if( newidx != idx2 )
    for( int i=0; i<(int)m_vvAllBelongingFV[idx2].size(); i++ ) {
      //new cluster assignment
      m_vBelongsTo[ m_vvAllBelongingFV[idx2][i] ] = newidx; 
    }
  
  vector<int> fvecBothAB( m_vvAllBelongingFV[idx1] );
  fvecBothAB.insert( fvecBothAB.end(), m_vvAllBelongingFV[idx2].begin(),
                     m_vvAllBelongingFV[idx2].end() );
//   for(int i=0; i < m_vvAllBelongingFV[b].size(); i++) {
//     fvecBothAB.push_back( m_vvAllBelongingFV[b][i] );
//   }
  m_vvAllBelongingFV[newidx] = fvecBothAB;
}


/*************************************************************/
/*       Call update function for merging in each step       */
/*************************************************************/
bool ClRNNCAgglo::doClusterSteps( double minSimilarity )
{
  // TIMING CODE

  m_dMinSimilarity = minSimilarity;
  m_bMore = true;
  cout << "  ClRNNCAgglo::doClusterSteps() called with min similarity: " 
       << m_dMinSimilarity << endl;

  cout << "  Clustering..." << endl;

  int nNumPoints = m_vCenters.size();
  vector<int>   vNN      ( nNumPoints );
  vector<bool>  vInactive( nNumPoints );
  vector<float> vSim     ( nNumPoints );

  /******************************/
  /* Start with the first point */
  /******************************/
  int last = 0;
  vNN[last] = 0;
  vInactive[vNN[last]] = true;
  vSim     [last]  = MIN_SIM;
  int nNumActive   = nNumPoints - 1;
  int nFirstActive = 1;

  /***************/
  /* Agglomerate */
  /***************/
  while( nNumActive > 0 ) {
    /*=====================================*/
    /* Find the NN for the last chain link */
    /*=====================================*/
    int   nnidx = nFirstActive;
    float nnsim = sim( vNN[last], nnidx );
    for( int i=nFirstActive+1; i<nNumPoints; i++ )
      if( !vInactive[i] ) {
        float d = sim( vNN[last], i );
        if( d > nnsim ) {
          nnsim = d;
          nnidx = i;
        }
      }

    /*==========================*/
    /* Check for reciprocal NNs */
    /*==========================*/
    if( nnsim > vSim[last] ) {
      /*-----------------------------------------*/
      /* No RNN => Add the point to the NN chain */
      /*-----------------------------------------*/
      last++;
      vNN[last]  = nnidx;
      vSim[last] = nnsim;
      vInactive[nnidx] = true;
      nNumActive--;

      /* recompute the 'first' counter */
      if( (nnidx==nFirstActive) && (nNumActive>0) )
        for( int i=nFirstActive+1; i<nNumPoints; i++ )
          if( !vInactive[i] ) {
            nFirstActive = i;
            break;
          }

    } else
      /*----------------------------------------------------------*/
      /* RNN => check if the chain links are sufficiently similar */
      /*----------------------------------------------------------*/
      if( vSim[last] >= m_dMinSimilarity ) {
	
        /*--------------------------------------*/
        /* Agglomerate the two last chain links */
        /*--------------------------------------*/
        int newidx = vNN[last]; //min( vNN[last], vNN[last-1] );
        agglomerate( vNN[last], vNN[last-1], newidx );
        writeTrace ( vNN[last], vNN[last-1], vSim[last], newidx );
        vInactive[newidx] = false;
        nNumActive++;
        last -= 2;
        
        /* recompute the 'first' counter */
        if( newidx < nFirstActive )
          nFirstActive = newidx;

      } else {
        /*-------------------*/
        /* Discard the chain */
        /*-------------------*/
        last = -1;
      }

    /*=============================*/
    /* Check if the chain is empty */
    /*=============================*/
    if( last < 0 ){
      /* empty chain => start new chain */
      last++;
      vNN[last]  = nFirstActive;
      vSim[last] = MIN_SIM;
      vInactive[nFirstActive] = true;
      nNumActive--;

      /* recompute the 'first' counter */
      if( nNumActive > 0 )
        for( int i=nFirstActive+1; i<nNumPoints; i++ )
          if( !vInactive[i] ) {
            nFirstActive = i;
            break;
          }
    } 
  }



  /***************************/
  /* Prepare the result data */
  /***************************/
  vector<ClPoint> vNewCenters;
  vector<int> vNewIndices;
  int nNumInvalid = 0;
  for( int i=0; i<(int)m_vCenters.size(); i++ ) {
    if ( m_vValid[i] ) {
      vNewCenters.push_back( m_vCenters[i] );
      vNewIndices.push_back( i - nNumInvalid );
    }
    else {
      vNewIndices.push_back( -1 );
      nNumInvalid++;
    }
  }
  m_vCenters = vNewCenters;
  
  /* correct the following indices */
  for( int i=0; i<(int)m_vBelongsTo.size(); i++ ){
    m_vBelongsTo[i] = vNewIndices[ m_vBelongsTo[i] ];
  }
  
  /* stop agglomerative clustering */
  cout << "  Size of m_vBelongsTo after is: " << m_vBelongsTo.size() << endl;
  cout << "  Size of m_vCenters after is: " << m_vCenters.size() << endl;
  cout << "  Agglomerative clustering finished." << endl;
  cout << endl;

  return true;
}



/***********************************************************/
/*                      Reestimation                       */
/***********************************************************/

void ClRNNCAgglo::doReestimationStep()
{}



/***********************************************************/
/*                         Update                          */
/***********************************************************/

void ClRNNCAgglo::doUpdateStep()
{
}


/***********************************************************/
/*                     Error Estimation                    */
/***********************************************************/

FLOAT ClRNNCAgglo::calculateError()
{
  return -1.0;
}


/***********************************************************/
/*                     Output Functions                    */
/***********************************************************/

void ClRNNCAgglo::printResults()
{
  cout << "  ClRNNCAgglo::printResults() called." << endl;

  /* calculate the number of points in each cluster */
  vector<int> nPointsInCluster( m_nNumClusters, 0 );
  for( int i=0; i<m_nNumPoints; i++ )
    nPointsInCluster[ m_vBelongsTo[i] ]++;

  for (int j=0; j<m_nNumClusters; j++) {
    cout << "    Cluster " << setw(2) << j << ":";
    if ( m_nDim <= 3 )
      cout << " mean=" << m_vCenters[j] << ",";
    cout << " points=" << setw(5) << nPointsInCluster[j] << endl;
  }
  cout << "  done." << endl;
}



void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  // must have 3  inputs.
  if (nrhs != 3) { printf("in nrhs != 3\n"); return; }
  if (mxGetClassID(prhs[0]) != mxCHAR_CLASS) { printf("comment input must be a char array\n"); return; }
  if (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) { printf("features input must be a float array\n"); return; }
  if (mxGetClassID(prhs[2]) != mxDOUBLE_CLASS) { printf("cluster number input must be an int\n"); return; }
  // must have 3 output.
  //if (nlhs != 3) { printf("out nlhs != 3\n"); return; }
  
  // get size of x.
  int const num_dims0 = mxGetNumberOfDimensions(prhs[0]);
  int const *dims0 = mxGetDimensions(prhs[0]);
  int const num_dims1 = mxGetNumberOfDimensions(prhs[1]);
  int const *dims1 = mxGetDimensions(prhs[1]);
  int const num_dims2 = mxGetNumberOfDimensions(prhs[2]);
  int const *dims2 = mxGetDimensions(prhs[2]);
  if( dims0[0]!=1){ printf("dims != 9\n"); return; }
 
  char const *name = (char *) mxArrayToString(prhs[0]);
  double const *features = (double *)mxGetData(prhs[1]);
  double const *nb_clustx = (double *)mxGetData(prhs[2]);
  uint nb_clust=(uint)nb_clustx[0];
  float dSimThresh=nb_clustx[1];
  int dim=(int)dims1[0];
  int nb_feat=(int)dims1[1];
  if(dim>nb_feat){cout <<"nb of dimensions is less than the number of features dim="<< dim << " > nb_feat=" << nb_feat << endl; return;}
  printf("input settins:\nname: %s\ndimensions: %d\nnumber of features: %d\nnumber of output clusters: %d\ncluster distance threshold: %f\n",name, dim,nb_feat,nb_clust,dSimThresh);
  

  vector<ClPoint> vPoints;

  
  for( int i=0; i<nb_feat; i++ ) {
    vPoints.push_back(ClPoint(dim));
    for( int v=0; v<dim; v++ ) {    
      vPoints[i].set(v,features[i*dim+v]);
      //     cout << features[i*dim+v] << endl;
    }
    //cout << endl;
  }


 ClRNNCAgglo clAgglo( vPoints );
 //int nNumFeatures = (int)vFeatures.size();
 clAgglo.initClusters( vPoints.size());

 clAgglo.doClusterSteps( dSimThresh );

 vector<ClStep> vClusterTrace;
 vClusterTrace      = clAgglo.getClusterTrace();
 
 int *odim = new int[2];
 odim[0]=3;
 odim[1]=vClusterTrace.size();
 plhs[0] = mxCreateNumericArray(2, odim, mxDOUBLE_CLASS, mxREAL);
 double *trace_out = (double *)mxGetData(plhs[0]);
 odim[0]=1;
 odim[1]=nb_feat;
 plhs[1] = mxCreateNumericArray(2, odim, mxDOUBLE_CLASS, mxREAL);
 double *ass_out = (double *)mxGetData(plhs[1]);
  
 
 vector<pair<float, pair<uint, uint> > > lines;  

  for( int i=0; i<vClusterTrace.size(); i++ ) {
    lines.push_back(make_pair(vClusterTrace[i].dSim,make_pair(vClusterTrace[i].nIdx1,vClusterTrace[i].nIdx2)));
  }
  
  sort(lines.begin(),lines.end(), greater< pair<float,pair<uint,uint> > >());

  vector<pair<uint, uint> > clust;
  vector<pair<uint, uint> > nclust;
  for(int i=0;i<nb_feat;i++){
    clust.push_back(make_pair(i,i));
    nclust.push_back(make_pair(i,i));
  }

  int l_nb = (nb_feat-nb_clust);  
  if(l_nb>vClusterTrace.size())l_nb=vClusterTrace.size();
  if(nb_clust<2)nb_clust=2;
  for(int i=0;i<l_nb;i++){
    clust[lines[i].second.second].first=lines[i].second.first;
    nclust[lines[i].second.second].first=lines[i].second.first;
  }
  
  int newentry;
  for(int i=0;i<clust.size();i++){
    if(clust[i].first!=i){
      newentry=clust[i].first;
      while(clust[newentry].first!=newentry)newentry=clust[newentry].first;
      nclust[i].first=newentry;
    }
  }  
  sort(nclust.begin(),nclust.end());
  
  int cl=1,prev=nclust[0].first;
  for( int i=0; i<nclust.size(); i++ ) {
    if(nclust[i].first==prev)
      nclust[i].first=cl;
    else {
      cl++;
      prev=nclust[i].first;
      nclust[i].first=cl;     
    }
    clust[i].first=nclust[i].second;
    clust[i].second=nclust[i].first;
  }
  sort(clust.begin(),clust.end());
  for( int i=0; i<clust.size(); i++ ) {
    ass_out[i]=clust[i].second;
  }

  for( int i=0; i<vClusterTrace.size(); i++ ) {
    //  cout  << vClusterTrace[i].nIdx1 << " "<< vClusterTrace[i].nIdx2<< " " << vClusterTrace[i].nNewIdx <<  " " << vClusterTrace[i].dSim <<  endl;
    trace_out[i*3]=lines[i].second.first+1;
    trace_out[i*3+1]=lines[i].second.second+1;
    trace_out[i*3+2]=lines[i].first;
  }


  delete []odim;
  vPoints.clear();
  vClusterTrace.clear();
  lines.clear();

}
