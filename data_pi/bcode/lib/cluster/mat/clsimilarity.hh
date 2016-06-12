/*********************************************************************/
/*                                                                   */
/* FILE         clsimilarity.hh                                      */
/* AUTHORS      Bastian Leibe, Andreas Remund                        */
/* EMAIL        leibe@inf.ethz.ch                                    */
/*                                                                   */
/* CONTENT      Defines a data type used for the clustering procedu- */
/*              dure when building up a priority queue.              */
/*                                                                   */
/*                                                                   */
/* BEGIN        Tue Oct 22 2002                                      */
/* LAST CHANGE  Fri Jul 18 2003                                      */
/*                                                                   */
/*********************************************************************/

#ifndef CLSIMILARITY_HH
#define CLSIMILARITY_HH


/*===================================================================*/
/*                         Class ClSimilarity                        */
/*===================================================================*/
/* Define an own data type holding the similarity between two clus-  */
/* ter centers and the corresponding indices.                        */
class ClSimilarity
{
public:
  ClSimilarity( double sim, int clA, int clB )
  {
    m_sim = sim;
    m_A = clA;
    m_B = clB;
  };

public:
  double m_sim;
  int m_A, m_B;
};


/*===================================================================*/
/*                         Class ClSimCompare                        */
/*===================================================================*/
/* Define a class (functor) which compares the similarity value used */
/* to sort the queue */
class ClSimCompare
{
public:
  ClSimCompare(){};
  bool operator() ( const ClSimilarity &A, const ClSimilarity &B ) {
    return ( A.m_sim <= B.m_sim );
  }

};

#endif
