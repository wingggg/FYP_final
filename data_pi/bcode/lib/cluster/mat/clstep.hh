/*********************************************************************/
/*                                                                   */
/* FILE         clstep.hh                                            */
/* AUTHORS      Bastian Leibe                                        */
/* EMAIL        leibe@informatik.tu-darmstadt.de                     */
/*                                                                   */
/* CONTENT      Define a helper class ClStep for keeping a record of */
/*              he clustering trace.                                 */
/*                                                                   */
/* BEGIN        Sun Mar 14 2005                                      */
/* LAST CHANGE  Sun Mar 14 2005                                      */
/*                                                                   */
/*********************************************************************/

#ifndef CLSTEP_HH
#define CLSTEP_HH

using namespace std;

/****************/
/*   Includes   */
/****************/
#include <vector>
#include <cassert>

#include "cluster.hh"

/*******************/
/*   Definitions   */
/*******************/
class ClStep 
{
public:
  ClStep( int idx1, int idx2, float sim, int newidx )
  {
    nIdx1   = idx1;
    nIdx2   = idx2;
    dSim    = sim;
    nNewIdx = newidx;
  }

  int   nIdx1, nIdx2;
  float dSim;
  int   nNewIdx;
};

/*-------------------------------------------------------------------*/
/*                         Sorting Operators                         */
/*-------------------------------------------------------------------*/
struct compClStepAsc
{
  bool operator()( ClStep x, ClStep y )
  { return (x.dSim > y.dSim); }
};

struct compClStepDesc
{
  bool operator()( ClStep x, ClStep y )
  { return (x.dSim < y.dSim); }
};

#endif
