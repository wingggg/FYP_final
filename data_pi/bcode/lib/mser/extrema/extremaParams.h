/*--------------------------------------------------------------------------*/
/* Copyright 2006, Jiri Matas & Michal Perdoch       matas@cmp.felk.cvut.cz */
/*--------------------------------------------------------------------------*/

#ifndef __EXTREMA_PARAMS_H__
#define __EXTREMA_PARAMS_H__

#define  GENERATE_MSER_PLUS    1
#define  GENERATE_MSER_MINUS   2

#include "extremaConfig.h"

namespace extrema
{

/* KEEP these constants with lowercase, otherwise correct preprocess.h */
//! An enumeration that encodes different preprocessings of an image.
enum EXTREMA_PREPROCESS
{
   PREPROCESS_CHANNEL_none          = 0x00000000,
   PREPROCESS_CHANNEL_intensity     = 0x00000001,
   PREPROCESS_CHANNEL_saturation    = 0x00000002,
   PREPROCESS_CHANNEL_hue           = 0x00000003,
   PREPROCESS_CHANNEL_redblue       = 0x00000004,
   PREPROCESS_CHANNEL_red           = 0x00000005,
   PREPROCESS_CHANNEL_green         = 0x00000006,
   PREPROCESS_CHANNEL_blue          = 0x00000007,
   PREPROCESS_CHANNEL_greenmagenta  = 0x00000008,
   PREPROCESS_CHANNEL_intensity_half= 0x00000009,
 
   PREPROCESS_CHANNEL_MASK          = 0x0000ffff,

   PREPROCESS_INTENSITY_none        = 0x00000000,
   PREPROCESS_INTENSITY_MASK        = 0xffff0000
};

//! A structure holding image parameters
struct ExtremaImage
{
   //! Width of the image.
   unsigned int    width;
   //! Height of the image.
   unsigned int    height;
   //! Number of channels of the image.
   unsigned int    channels;
   //! Pointer to image data.
   unsigned char * data;
};

//! A structure with MSER detector parameters.
struct ExtremaParams
{
   bool   relative;
   int    preprocess; /* see EXT_PREPROCESS enum */
   int    min_size;
   double max_area;
   double min_margin;
   bool   verbose;
   int    debug;

   ExtremaParams()
      {
         relative=false;
         preprocess=PREPROCESS_CHANNEL_none;
         min_size=0;
         max_area=0.01;
         min_size=30;
         min_margin=10;
         verbose=0;
         debug=0;
      }
};
}

#endif
