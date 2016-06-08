// Definition de la classe ImageContent
#ifndef _imageContent_h_
#define _imageContent_h_

//#define WITH_LIBPNG
//#define WITH_LIBJPEG
//#define WITH_LIBTIFF

#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <png.h>
#include "imageread.h"
#include "../util/util.h"


#define COL 5
#define RGB 2
#define OPP 3
#define GRAY 1
#define FLOAT1 2
#define UCHAR1 1
#define UCHAR1FLOAT1 3
#define FLOAT3 8
#define UCHAR3 5
#define UCHAR3FLOAT3 12
#define UCHAR3FLOAT1 14
#define OPP3FLOAT1 15

const int PNG_BYTES_TO_CHECK = 4;
const int ERROR = -1;

class ImageContent;
typedef ImageContent DARY;
//typedef unsigned int uint;
//typedef unsigned char uchar;


namespace imgcon
{
struct ConImage
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
}


class ImageContent {
   
  private :
    uint x_size,y_size;
    uint tsize;
    void writePGM(const char *nom,unsigned char *buff, const char* comments);
    void write(const char *nom, const char* comments);
    void init(uint ,uint, int);	   
    int buftype,coltype;
  
  public :
    float **fel;
    float **felr;
    float **felg;
    float **felb;
    char *filename;
    unsigned char **bel;
    unsigned char **belr;
    unsigned char **belg;
    unsigned char **belb;
    ImageContent(void){ y_size=0;x_size=0;buftype=0;};   
    ImageContent(const char *);	   
    ImageContent(const char *, char *type);	   
    //ImageContent(const char *, int type);	   
    ImageContent(ImageContent *im);	   
//    ImageContent(uint y_size_in,uint x_size_in){init( y_size_in, x_size_in,FLOAT1);};	   
//    ImageContent(int y_size_in ,int x_size_in){init((uint)y_size_in,(uint)x_size_in, FLOAT1);};	   
    ImageContent(uint ,uint, int);	   
    ImageContent(uint y_size_in, uint x_size_in, int buftype_in, float val){  y_size=0;x_size=0;buftype=0;init( y_size_in, x_size_in,buftype_in);buftype=buftype_in;set(val);};	   
 
    ~ImageContent();

    inline uint  const  x() const { return x_size;}
    inline uint  const  y() const { return y_size;}
    inline uint  const  size() const { return tsize;}
    int const getType() const { return buftype;}
    int const getColType() const { return coltype;}
    const char* name(){return filename;}
    void write(const char *nom);
    void writePNG(const char* name);
    void writePNG(const char* name,unsigned char **buff);
    void writeR(const char *nom);
    void writeG(const char *nom);
    void writeB(const char *nom);
    void writeF(const char* name);
    void RGB2xyY();
    void RGB2lbrg();
    void RGB2rgb();
    void RGB2opp();
    void flipH();
    void flipV();
    void convert(int type);
    void set(float);
    void set(DARY*);
    void set(unsigned char*, uint y, uint x);
    void set(float*, uint y, uint x);
    void set(const char *name){strcpy(filename,name);}
    void set(ImageContent *r, ImageContent *g, ImageContent *b);
    void get(ImageContent *r, ImageContent *g, ImageContent *b);

    void normalize(float min_in, float max_in);
    void normalize();
    void scale(DARY *im_in, float scalex, float scaley);
    void scale_float(float**im_in, int xim_in, int yim_in, float**im_out, float scalex, float scaley);
    void scale_uchar(unsigned char**im_in, int xim_in, int yim_in, unsigned char**im_out, float scalex, float scaley);
    void decrease(unsigned char **ucharArray, int xim_in, int yim_in, unsigned char**im_out);
    void decrease(DARY *im_in);
    float getValue(float x, float y);
    void interpolate(DARY *sface, float m_x, float m_y, 
                     float scalex, float scaley, float angle);
    void interpolate(DARY *im_in, float m_x, float m_y, float vec0x, float vec0y,
                     float vec1x, float vec1y);
    void median2d(ImageContent *output, int median_s);
    void erosion(ImageContent *output, int elem_rad, float val);
    void dilation(ImageContent *output, int elem_rad, float val);
    void crop(DARY *img, int x, int y, char *mode);
    void crop(DARY *img, int x, int y);
};
       


 
#endif
