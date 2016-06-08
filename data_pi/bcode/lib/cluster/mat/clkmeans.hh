/*********************************************************************/
/*                                                                   */
/* FILE         clkmeans.hh                                          */
/* AUTHORS      Bastian Leibe                                        */
/* EMAIL        leibe@inf.ethz.ch                                    */
/*                                                                   */
/* CONTENT      Define a general cluster class that serves as a ba-  */
/*              sis for derivative classes implementing the specific */
/*              clustering algorithms.                               */
/*                                                                   */
/* BEGIN        Tue Sep 04 2001                                      */
/* LAST CHANGE  Tue Sep 04 2001                                      */
/*                                                                   */
/*********************************************************************/

#ifndef LEIBE_CLKMEANS_HH
#define LEIBE_CLKMEANS_HH

#ifdef _USE_PERSONAL_NAMESPACES
//namespace Leibe {
#endif
  
/****************/
/*   Includes   */
/****************/
#include <vector>
#include <cassert>

#include "cluster.hh"

/*************************/
/*   Class Definitions   */
/*************************/

/*===================================================================*/
/*                         Class ClKMeans                            */
/*===================================================================*/
/* Define a k-means clustering class derived from the Cluster class. */
class ClKMeans : public Cluster
{
public:
  ClKMeans( vector<ClPoint> &vPoints ) 
    : Cluster( vPoints ) {}

public:
  /*******************************/
  /*   Content Access Functions  */
  /*******************************/
  vector<int> getClusterAssignment() { return m_vBelongsTo; }

public:
  /****************************/
  /*   Clustering Functions   */
  /****************************/

protected:
  virtual void  initDataVectors();
  virtual void  initClusterCenters();
  virtual void  initCovariances();
  virtual void  initPriors();
  virtual void  initPosteriors();

  virtual void  doReestimationStep();
  virtual void  doUpdateStep();
  virtual FLOAT calculateError();

public:
  /************************/
  /*   Output Functions   */
  /************************/
  virtual void  printResults();

private:
  vector<int> m_vBelongsTo;
};

#ifdef _USE_PERSONAL_NAMESPACES
//}
#endif

#endif // LEIBE_KMEANS_HH
