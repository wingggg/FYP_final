/*--------------------------------------------------------------------------*/
/* Copyright 2006, Jiri Matas & Michal Perdoch       matas@cmp.felk.cvut.cz */
/*--------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include "extremaTypes.h"
#include "getExtrema.h"
#include "boundary.h"
#include "suballoc.h"


using namespace utls;

namespace extrema
{


/*----------------------------------------------------------------*/

static IAry * bck   = 0;
static BAry * image = 0;

#define is_inside(i, j, max_int) \
        (!(bck->el[j][i]<0) && image->el[j][i] <= max_int)


/*----------------------------------------------------------------*/
inline static void add_if_inside(point_vector *inside, point_vector *boundary, int max_int, int i, int j, unsigned char direct, int marker)
{ 
   if (marker == bck->el[j][i])
      return;
   t_borderpixel p;

   p.pos.x  = i; p.pos.y  = j;
   p.direct = direct;

   if (is_inside(i, j, max_int))
   { 
      inside->push_back(p);
      bck->el[j][i] = marker;
   } 
   else
      boundary->push_back(p);
}

/*----------------------------------------------------------------*/

static void Region(t_region* p_r)
{ 
   point_vector  tmp;
   point_vector  inside;
   int           marker = p_r->label; 
   int           last_thresh=0;
   t_thresh_def  *p_t;

   ForeachTyLL_M(p_r->thresholds, p_t, t_thresh_def*)
      {
         int thresh             = p_t->thresh; 
         if (last_thresh>=thresh)
         {
            printf("Threshold out of order, ignored.\n");
            continue;
         }
         point_vector *boundary = p_t->boundary = new point_vector;
        
         if (IsFirstElmLL(p_t))
         {
            /* for the first threshold of the region, insert first
               point to the inside list */
            t_borderpixel p;
            p.pos = p_r->minimum_pos;
            p.direct = 0;
            /* insert pixel into list and mark it with label */
            inside.push_back(p);
            bck->el[p.pos.y][p.pos.x] = marker;
         }
         else
         {
            /* for the other thresholds of the region, check if last
               thresholds boundary is still inside and mark background 
               pixels appropriately, leave them out otherwise */
            while (!inside.empty())
            {
               t_borderpixel &curr = inside.back();
               if (!is_inside(curr.pos.x, curr.pos.y, thresh))
                  // leave out this point
                  boundary->push_back(curr);
               else {
                  // mark pixel
                  if (bck->el[curr.pos.y][curr.pos.x] != marker)
                  {
                     tmp.push_back(curr);
                     bck->el[curr.pos.y][curr.pos.x] = marker;
                  }
               }
               inside.pop_back();
            }
            inside = tmp;
            tmp.clear();
         }

         while (!inside.empty())
         {
            t_borderpixel curr = inside.back(); inside.pop_back();
            int i      =  curr.pos.x;
            int j      =  curr.pos.y;
            add_if_inside(&inside, boundary, thresh, i , j+1, 4, marker);
            add_if_inside(&inside, boundary, thresh, i , j-1, 1, marker);
            add_if_inside(&inside, boundary, thresh, i+1,  j, 2, marker);
            add_if_inside(&inside, boundary, thresh, i-1,  j, 8, marker);
         }
         
         /* sort region's boundary */
         std::sort(boundary->begin(), boundary->end());
         
         /* if it is not last threshold, copy boundary into inside list */
         if (!IsLastElmLL(p_t))
            inside = *boundary;
      }
}

/*----------------------------------------------------------------*/
void RegionBoundaries(BAry *img, t_LL regions)
{
   int i, cols, rows;
   cols = img->cols()-2; 
   rows = img->rows()-2;

   image = img;
   /* prepare background image */
   bck = new IAry(img->lb1, img->ub1, img->lb2, img->ub2);

   memset(&bck->el[-1][-1], 0, sizeof(int)*(rows+2)*(cols+2));
   for (i=0; i<cols; i++) bck->el[-1][i] = bck->el[rows][i] = -1;
   for (i=0; i<rows; i++) bck->el[i][-1] = bck->el[i][cols] = -1;
   
   t_region * p_r;
   ForeachTyLL_M(regions, p_r, t_region*)
   {
      if (p_r->thresholds && !IsEmptyLL(p_r->thresholds))
         Region(p_r);
   }
   delete bck;
}

}
