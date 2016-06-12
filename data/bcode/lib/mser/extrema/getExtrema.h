/*--------------------------------------------------------------------------*/
/* Copyright 2006, Jiri Matas & Michal Perdoch       matas@cmp.felk.cvut.cz */
/*--------------------------------------------------------------------------*/

#ifndef __GET_EXTREMA_H__
#define __GET_EXTREMA_H__

#include "../LL/LL.h"
#include "../utls/ary.h"
#include "extremaParams.h"
#include "extremaTypes.h"
#include "suballoc.h"

namespace extrema 
{

extern t_thresh_par g_thresh_params;

void InitRegionRecycling();
void DestRegionRecycling();
t_LL GetExtrema(utls::BAry* img, t_sortpixels pixels, const ExtremaParams &par, bool invert);
void DestRegions(t_LL regions);

}
#endif
