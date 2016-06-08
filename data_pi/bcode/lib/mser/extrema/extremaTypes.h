/*--------------------------------------------------------------------------*/
/* Copyright 2006, Jiri Matas & Michal Perdoch       matas@cmp.felk.cvut.cz */
/*--------------------------------------------------------------------------*/

#ifndef __EXTREMA_TYPES_H__
#define __EXTREMA_TYPES_H__

#include "../LL/LL.h"
#include <vector>
#include "extremaConfig.h"

#define c_maxByte 256

namespace extrema
{
   typedef unsigned int t_label;

   //! Internal structure, holds 2D point coordinates.
   typedef struct
   {
      int x;
      int y;
   } t_ipoint;

   //! Internal structure with intensity histogram.
   typedef struct s_sortpixels
   {
      t_ipoint * data[c_maxByte];
      int        hist[c_maxByte];
   } t_sortpixels;

   typedef unsigned int t_mregion;

   //! Internal region structure.
   typedef struct s_region
   {
      t_label   label;
      int       minimum_int;
      int       pixel_total;
      int       border_total;
      t_ipoint  minimum_pos;
      int       maximum_int;
      t_label   merge_label;
      t_LL      thresholds;
      int       pixels[c_maxByte];
      int       borders[c_maxByte];
   } t_region;

   //! Internal structure with a node of the label equivalency tree
   typedef struct s_region_equiv
   {
      unsigned int   pred;
      t_region     * region;
   } t_region_equiv;

   //! Internal structure with processed detector's parameters.
   typedef struct s_thresh_par
   {
      //! minimum size of the region in pixels
      int    min_size;

      //! maximum size of the region in pixels
      int    max_size;

      //! minimum margin and upper boundary for hystheresis thresholding
      double min_margin;

      //! margin relative to intesity level
      bool   relative_margin;

      //! do inverted margin
      int    invert;
   } t_thresh_par;

   //! Structure with pixel of the extended boundary.
   typedef struct s_borderpixel
   {
      t_ipoint      pos;
      unsigned char direct; // N = 1, E = 2, S = 4, W = 8; 0 = unknown direction
      bool operator<(const s_borderpixel &other) const 
         {
            if (pos.y > other.pos.y) return false;
            if (pos.y < other.pos.y) return true;
            if (pos.x < other.pos.x) return true;
            return false;
         }
   } t_borderpixel;

   //! Vector with extended boundary.
   typedef std::vector<t_borderpixel> point_vector;

   //! Internal structure holding threshold paramaters.
   typedef struct s_thresh_def
   {
      int            thresh;
      int            pos;
      int            margin;
      point_vector  *boundary;
   } t_thresh_def;

}
#endif
