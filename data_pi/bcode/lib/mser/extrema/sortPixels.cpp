/*--------------------------------------------------------------------------*/
/* Copyright 2006, Jiri Matas & Michal Perdoch       matas@cmp.felk.cvut.cz */
/*--------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include "sortPixels.h"

using namespace utls;

namespace extrema 
{

void CalcHistogram(BAry * &img, t_sortpixels &pixels)
{
   int i, j, rows = img->rows()-2, cols = img->cols()-2;
   for(i=0; i < c_maxByte; i++)
      pixels.hist[i] = 0;
   unsigned char *src = &img->el[0][0];
   for (i=0; i < rows; i++)
   {
      for (j=0; j < cols; j++)
         pixels.hist[*src++]++;
      src+=2;
   }
}
void BinSortPixels(BAry * &img, t_sortpixels &pixels)
{
   int i, j;
   t_ipoint pos;
   t_ipoint *last[c_maxByte];
   /* initialize pixels counts */
   for(i=0; i < c_maxByte; i++)
      last[i] = 0;

   /* allocate memory for pixel coordinates */
   for(i=0; i < c_maxByte; i++)
   {
      if (pixels.hist[i])
         last[i] = pixels.data[i] = 
            (t_ipoint *) malloc((pixels.hist[i] + 1) * sizeof(t_ipoint));
      else
         pixels.data[i] = last[i];
   }

   /* fill pixel lists with positions */
   for (i=0; i < img->rows()-2; i++)
      for (j=0; j < img->cols()-2; j++)
      {
         pos.x = j; pos.y = i;
         *last[img->el[i][j]]++ = pos;
      }
}

/* inverts padded intensity image and it's intesity histogram */
void InvertImageAndHistogram(BAry *img, t_sortpixels &pixels)
{ 
   int i, rows = img->rows()-2, cols = img->cols();
   unsigned char *src = &img->el[0][0];
   for (i=0; i < rows*cols; i++)
   {
      *src = 255 - *src; 
      src++;
   }
   
   int       auxHist;
   t_ipoint* auxData;
   for(i=0; i < (c_maxByte / 2) ; i++)
   {
      auxHist = pixels.hist[i];
      pixels.hist[i] = pixels.hist[c_maxByte - i - 1];
      pixels.hist[c_maxByte - i - 1] = auxHist;
      
      auxData = pixels.data[i];
      pixels.data[i] = pixels.data[c_maxByte - i - 1];
      pixels.data[c_maxByte - i - 1] = auxData;
   }
}

}
