/*--------------------------------------------------------------------------*/
/* Copyright 2006, Jiri Matas & Michal Perdoch       matas@cmp.felk.cvut.cz */
/*--------------------------------------------------------------------------*/

#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include "getExtrema.h"
#include "optThresh.h"


using namespace utls;

namespace extrema {

/* exported variables */
t_thresh_par g_thresh_params;
t_region_equiv * label_equiv = NULL;

/*----------------------------------------------------------------*/

static unsigned int  current_label;
static IAry*         labels;
static int           border_num = 0;


t_suballocator *regionSuballocator = 0;

/*----------------------------------------------------------------*/

void InitRegionRecycling()
{
   if (regionSuballocator==0)
   {
      int size = sizeof(t_region);
      regionSuballocator = (t_suballocator *)malloc(sizeof(t_suballocator));
      InitSuballocator(regionSuballocator, 100, size, 1);
   }
}

/*----------------------------------------------------------------*/

void DestRegionRecycling()
{
   if (regionSuballocator!=0)
   {
      DestSuballocator(regionSuballocator);
      free(regionSuballocator);
      regionSuballocator = 0;
   }
}

int min(int a, int b)
{
   return a<b?a:b;
}

int max(int a, int b)
{
   return a>b?a:b;
}


static unsigned int FindEquivLabel(unsigned int label)
{
   /* deal with one level case faster (no relabeling) */
   unsigned int tmp = label_equiv[label].pred;

   if (tmp == label_equiv[tmp].pred)
      return tmp;

   /* find real label (the root of the labels merging tree) */
   do {
      tmp = label_equiv[tmp].pred;
   } while (label_equiv[tmp].pred != tmp);

   unsigned int final = tmp;
   tmp                = label;

   /* flatten this branch of the tree to the first level */
   while (label_equiv[tmp].pred != tmp)
   {
      tmp                        = label_equiv[tmp].pred;
      label_equiv[label].pred    = final;
      label                      = tmp;
   }
   return final;
}

/*----------------------------------------------------------------*/

inline static unsigned int FindLabel(unsigned int label)
{
   /*  deal with the common case (label's predcessor points to itself) fast */
   if (label==0 || label_equiv[label].pred == label)
      return label;

   return FindEquivLabel(label);
}

/*----------------------------------------------------------------*/

/* add pixel on position POS with intensity INTENSITY to the region PREGION */
inline static void InsMarkPixel(t_region *pRegion, t_ipoint pos, int intensity)
{
   pRegion->maximum_int = intensity;

   /* update pixel counts */
   pRegion->pixels[intensity]  ++;
   pRegion->pixel_total        ++;

   /* update border counts */
   pRegion->borders[intensity] += border_num;	/*H*/
   pRegion->border_total       += border_num; 	/*H*/

   /* mark pixel with region's label */
   labels->el[pos.y][pos.x] = pRegion->label;
}

/*----------------------------------------------------------------*/
static void ConsRegion(t_ipoint minimum_pos, int intensity, t_LL regions)
{
   /* if there are some already created regions, reuse one of them */
   t_region *pRegion = (t_region*)SuballocatorGetItem(regionSuballocator);

   /* assign new label */
   pRegion->label       = ++current_label;
   pRegion->merge_label = -1;

   /* setup current intensity as initial */
   pRegion->minimum_int = intensity;
   pRegion->maximum_int = intensity;

   /* remember location of first pixel */
   pRegion->minimum_pos = minimum_pos;

   /* setup counters */
   pRegion->pixel_total  = 0;
   pRegion->border_total = 0;

   pRegion->pixels[intensity]  = 0;
   pRegion->borders[intensity] = 0;
   pRegion->thresholds = 0;

   /* insert region to the regions list */
   LinkInsLastLL(regions, *pRegion);

   /* store the history of labelling */
   label_equiv[current_label].pred   = current_label;

   /* tie the region with current label */
   label_equiv[current_label].region = pRegion;

   /* insert first pixel into region */
   InsMarkPixel(pRegion, minimum_pos, intensity);
}

/*----------------------------------------------------------------*/

static unsigned int labelled[4];
static int          label_num = 0;

/* returns the set of 4-connected labels in LABELLED global variable, size of
   the set is in LABEL_NUM, each label is inserted only once */

static void GetLabelled(t_ipoint pos)
/* the labels take margin into account */
/* assumes pos is not on the image edge, handled in sortedpixels */
{
   int x = pos.x, y = pos.y;

   /* image labels do not reflect merging; get the true, merged labels */
   unsigned int l1 = FindLabel(labels->el[y+1][x]);
   unsigned int l2 = FindLabel(labels->el[y][x+1]);
   unsigned int l3 = FindLabel(labels->el[y-1][x]);
   unsigned int l4 = FindLabel(labels->el[y][x-1]);

   label_num = 0;

   if (l1!=0)                                labelled[label_num++] = l1;
   if (l2!=0 && l2!=l1)                      labelled[label_num++] = l2;
   if (l3!=0 && l3!=l1 && l3!=l2)            labelled[label_num++] = l3;
   if (l4!=0 && l4!=l1 && l4!=l2 && l4!=l3)  labelled[label_num++] = l4;

   /* compute number of 4-connected neighbours */
   border_num = (l1 != 0) + (l2 != 0) + (l3 != 0) + (l4 != 0);	/*H*/
   border_num = 4 - 2*border_num;				/*H*/
}

/*----------------------------------------------------------------*/

static void MergeRegions(t_ipoint pos, int intensity)
{
   unsigned int   maxSize  = 0;
   unsigned int   maxLabel = 0;
   t_region     * pMaxRegion;
   int            i;

   /* find the region that was the largest at level: intensity-1
      this strategy makes sure that region with smallest size increase
      survives the merge */
   for(i=0; i < label_num; i++)
   {
      unsigned int label = labelled[i];
      t_region *region   = label_equiv[label].region;
      unsigned int size  = region->pixel_total - region->pixels[intensity];
      if (size >= maxSize)
      {
         maxSize  = size;
         maxLabel = label;
      }
   }

   /* finalise small regions, add their pixels totals to the largest */
   pMaxRegion = label_equiv[maxLabel].region;
   InsMarkPixel(pMaxRegion, pos, intensity);

   for(i=0; i < label_num; i++)
   {
      unsigned int label = labelled[i];
      if (label != maxLabel)
      {
         t_region *pRegion       = label_equiv[label].region;
         label_equiv[label].pred = maxLabel;

         /* merge region stats with the largest one */
         pMaxRegion->pixel_total        += pRegion->pixel_total;
         pMaxRegion->border_total       += pRegion->border_total;
         pMaxRegion->pixels[intensity]  += pRegion->pixel_total;
         pMaxRegion->borders[intensity] += pRegion->border_total;

         /* send regions smaller than MINSIZE to recycle list */
         if (pRegion->pixel_total < g_thresh_params.min_size)
            SuballocatorReturnItem(regionSuballocator, UnlinkLL(pRegion));
         else {
            /* set the final intensity and find thresholds */
            pRegion->maximum_int = intensity;
            pRegion->merge_label = maxLabel;
            FastSetOptThresholds4StableRegion(pRegion);
         }
      }
   }
}

/*----------------------------------------------------------------*/

void PrepareThresholds(BAry *img, const ExtremaParams &par, t_thresh_par &thr_par, bool invert)
{
   thr_par.min_size = par.min_size;
   thr_par.max_size = (int)((img->cols()-2) * (img->rows()-2) * par.max_area);
   thr_par.min_margin = par.min_margin;
   if (par.relative) thr_par.min_margin /= 100.0;
   thr_par.invert = invert;
   thr_par.relative_margin = par.relative;
}

/*----------------------------------------------------------------*/
t_LL GetExtrema(BAry* img, t_sortpixels pixels, const ExtremaParams &par, bool invert)
{
   int     i, rows, cols;
   t_LL    regions = ConsLL();

   cols = img->cols(); rows = img->rows();
   PrepareThresholds(img, par, g_thresh_params, invert);

   current_label   = 0;
   label_equiv     = (t_region_equiv *) 
      malloc(rows * cols * sizeof(t_region_equiv)/4);

   /* initialize label equivalency tree */
   label_equiv[current_label].pred = current_label;

   /* allocate array for labels */
   labels = new IAry(img->lb1, img->ub1, img->lb2, img->ub2);

   /* setup all pointers to zero */
   memset(&labels->el[-1][-1], 0, sizeof(int)*rows*cols);

   /* create suballocator and initialize memory */
   InitRegionRecycling();

   /* for all intensity levels */
   for(i=0; i < c_maxByte; i++)
   {
      t_ipoint *pend = pixels.data[i] + pixels.hist[i];
      /* all pixels of this intensity */
      for (t_ipoint *pelm = pixels.data[i]; pelm < pend; pelm++)
      {
         /* get the label set around pixel */
         GetLabelled(*pelm);

         switch(label_num)
         {
         case 0:
            ConsRegion(*pelm, i, regions);
            break;
         case 1:
            InsMarkPixel((t_region *)label_equiv[labelled[0]].region, *pelm, i);
            break;
         default:
            MergeRegions(*pelm, i);
         }
      }
   }

   /* process the last region (root), take any label and find root */
   t_label l = FindEquivLabel(label_equiv[labels->el[0][0]].pred);
   FastSetOptThresholds4StableRegion((t_region*)label_equiv[l].region);

   free(label_equiv);
   delete labels;
   return regions;
}

void DestRegions(t_LL regions)
{
   t_region * pRegion;

   ForeachTyLL_M(regions, pRegion, t_region *)
      {
         if (pRegion->thresholds)
         {
            t_thresh_def *p_t;
            ForeachTyLL_M(pRegion->thresholds, p_t, t_thresh_def *)
               {
                  if (p_t->boundary)
                     delete p_t->boundary;
               }
            DestLL(pRegion->thresholds);
         }
      }
   /* avoid double free in DestLL(regions); */
   SuballocatorReturnItemsLL(regionSuballocator, regions);
   free(regions);
}

}
