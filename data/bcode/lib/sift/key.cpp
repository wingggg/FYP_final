/* Copyright (c) 2000. David G. Lowe, University of British Columbia.
   This code may be used, distributed, or modified only for research
   purposes or under license from the University of British Columbia.
   This notice must be retained in all copies.
*/
/*
 *  Find the SIFT keys in an image.  
 *     Code is by David Lowe, University of British Columbia.
 *     See paper for details at:
 *        http://www.cs.ubc.ca/spider/lowe/papers/iccv99-abs.html
 *     Call routine GetKeypoints(image) to return all keypoints for image.
 */

#include "sift.h"
#include <stdlib.h>        
                
/* ------------------------- Global variables ---------------------------- */
/* The following global variables are the parameters used to control
      the keyframe finding.
*/
int DoubleImSize = FALSE;    /* Use double size image to find keypoints. */
float KeyDisplaySize = 0.4;  /* Display output size for key squares. */

int MinImageSize = 10;    /* Smallest scale must have at least this size */
float PeakThresh = 0.1;   /* Normalized peaks must be above this threshold. */
float GradThresh = 0.1;   /* Gradients are thresholded at this maximum */
float InsigMag = 0.0;     /* Grad mags below this do not need to be indexed. */
const int MagFactor = 4;  /* Magnification from grad to index sampling. */
float ScaleWeight = 2.0;  /* Relative weight for larger scale features. */
float InitSigma = 1.4;    /* Still used in test.c */
const int NormalizeVector = TRUE;  /* Norm vecs: needed for dist thresh */

/* The default threshold on ratio by which points in circle of 3 pixels radius
   should dip lower than current point (value should be less than 1.0, but
   can be raised much above 1 to eliminate check).
*/
float Dip = 0.95;  /* Was 0.95 */

/* Number of bins in orientation histogram (10 degree spacing). */
const int OriBins = 36;

/* Size of Gaussian used to select orientations.
*/
float OriSigma = 3.0;

/* This buffer is used to hold a row or column of the image to perform
   efficient convolution.
*/
#define BufferSize 2100
float Buffer[BufferSize];

/* This is the list of keypoints that has been found in current image.
*/
Keypoint Keypoints;
int KeyCount;


/* -------------------- Local function prototypes ------------------------ */

void FindKeypoints(Image image);
Image DoubleSize(Image image);
Image ReduceSize(Image image);
void SubtractIm10(Image im1, Image im2, Image result);
void HorConvSqrt2(Image image, Image result);
void VerConvSqrt2(Image image, Image result);
void FindMaxMin(Image *gauss, Image *dog, Image *ori, int levels);
void GradOriImages(Image im, Image ori);
int LocalMaxMin(float val, Image image, int row, int col, int level);
int CheckDip(float val, Image image, int row, int col);
int CheckDipPoint(float val, Image image, float row, float col);
void InterpKeyPoint(int level, Image dog, Image grad,
		    Image ori, Image grad2, Image ori2, int r, int c, float val);
float InterpPeak(float a, float b, float c);
void MakeKeypoint(float scale, float ori, float r, float c, float row,
	  float col, Image grad, Image orim, Image grad2, Image orim2, float val);
float FindOri(Image grad, Image ori, int row, int col);
void SmoothHistogram(float *hist, int bins);
float FindOriPeaks(float *hist, int bins);


/* ----------------------------- Routines --------------------------------- */


/* Given an image, find the keypoints and return a pointer to a list of
   keypoint records.
*/

void fastSIFT(DARY *image, vector<FD *> &features, Params *params){
  Image im;
      im = CreateImage( image->y(),  image->x(), 0);
      //memcpy(image->fel[0],im->pixels[0],image->size()*sizeof(float));
      for (int r = 0; r < im->rows; r++)
	for (int c = 0; c < im->cols; c++) {
	  im->pixels[r][c] = (float)image->fel[r][c];
    }
    float radius = (float)params->getValue("descriptor_radius.int");
    float dog_threshold = (float)params->getValue("dog_threshold.float");
    unsigned long type=OSIFT;
    Keypoints = NULL;  /* Global list of keypoints. */
    KeyCount = 0;
    /* If flag is set, then double the image size prior to finding
       keypoints.  Output coordinates will need to be reduced by 2. */
    if (DoubleImSize)
        im = DoubleSize(im);
    FindKeypoints(im);
    printf( "Found %d keypoints\n", KeyCount);
    Keypoint kn=Keypoints;
    while(kn!=NULL){
      if(fabs(kn->extr)>dog_threshold){
        FD *fd = new FD();
        fd->setX_Y(kn->col,kn->row);
        fd->setScale(kn->scale*radius/2.0);
        //fd->setAngle(kn->ori);
        fd->setAngle(1000);
        fd->setFeatureness(kn->extr);
	if((kn->extr)>0)
	  fd->setType(type);
	else fd->setType(type | (type>>1));
	//cout << fd->getType()<< endl;
        features.push_back(fd);
      }
      kn=kn->next;
    }

    
}


FILE *fd;
Keypoint GetKeypoints(Image image ,char *name)
{
    Keypoints = NULL;  /* Global list of keypoints. */
    KeyCount = 0;
    fd=fopen (name, "w+");
    if(!fd)fprintf(stderr, "file open error\n");
    fprintf(fd,"10\n0\n");
    

    /* If flag is set, then double the image size prior to finding
       keypoints.  Output coordinates will need to be reduced by 2. */
    if (DoubleImSize)
        image = DoubleSize(image);

    FindKeypoints(image);
     fclose(fd);
    fprintf(stderr, "Found %d keypoints\n", KeyCount);
    return Keypoints;
}

/* Compute keypoint locations and add to the global Keypoints list.
*/
void FindKeypoints(Image image)
{
    int lev = 0, first = TRUE;
    float scale = 1.0;
    Image gauss[20], dog[20], ori[20], blur, temp, smaller;

    /* Smooth input image with Gaussian of sqrt(2). */
    temp = CreateImage(image->rows, image->cols, IMAGE_POOL);
    HorConvSqrt2(image, temp);
    VerConvSqrt2(temp, image);

    /* Create each level of pyramid.  Keep reducing image by factors of 1.5
       until one dimension is smaller than minimum size.
    */
    while (image->rows > MinImageSize &&
	   image->cols > MinImageSize && lev < 20) {

      /* Smooth image by sqrt(2) and subtract to get DOG image. */
      blur = CreateImage(image->rows, image->cols, IMAGE_POOL);
      HorConvSqrt2(image, temp);
      VerConvSqrt2(temp, blur);

      SubtractIm10(image, blur, temp);
      dog[lev] = temp;

      /* Generate next level of pyramid by reducing by factor of 1.5. */
      smaller = ReduceSize(blur); 
      //writeImage(blur,"blur.pgm");fprintf(stderr,"OK\n");getchar();
      /* Find image gradients and orientations.  Reuse image and blur.
         This does not need to be done for lowest image in pyramid. */
      gauss[lev] = image;
      ori[lev] = blur;
      if (! first)
	GradOriImages(gauss[lev], ori[lev]);
      first = FALSE;

      scale *= 1.5;
      lev++;
      image = smaller;
      temp = CreateImage(image->rows, image->cols, IMAGE_POOL);
    }

    /* Find the keypoints in the pyramid. */
    FindMaxMin(gauss, dog, ori, lev);
}


/* Double image size. Use linear interpolation between closest pixels.
   Size is two rows and columns short of double to simplify interpolation.
*/
Image DoubleSize(Image image)
{
    int rows, cols, nrows, ncols, r, c, r2, c2;
    float **im, **newe;
    Image newimage;

    rows = image->rows;
    cols = image->cols;
    nrows = 2 * rows - 2;
    ncols = 2 * cols - 2;
    newimage = CreateImage(nrows, ncols, IMAGE_POOL);
    im = image->pixels;
    newe = newimage->pixels;

    for (r = 0; r < rows - 1; r++)
      for (c = 0; c < cols - 1; c++) {
	r2 = 2 * r;
	c2 = 2 * c;
	newe[r2][c2] = im[r][c];
	newe[r2+1][c2] = 0.5 * (im[r][c] + im[r+1][c]);
	newe[r2][c2+1] = 0.5 * (im[r][c] + im[r][c+1]);
	newe[r2+1][c2+1] = 0.25 * (im[r][c] + im[r+1][c] + im[r][c+1] +
				  im[r+1][c+1]);
      }
    return newimage;
}


/* Reduce the size of the image by resampling with a pixel spacing of
   1.5 times original spacing.  Each block of 9 original pixels is
   replaced by 4 new pixels resampled with bilinear interpolation.
*/
Image ReduceSize(Image image)
{
    int nrows, ncols, r, c, r1, r2, c1, c2;
    float **im, **newe;
    Image newimage;

    nrows = (2 * image->rows) / 3;
    ncols = (2 * image->cols) / 3;
    newimage = CreateImage(nrows, ncols, IMAGE_POOL);
    im = image->pixels;
    newe = newimage->pixels;

    for (r = 0; r < nrows; r++)
      for (c = 0; c < ncols; c++) {
	if (r % 2 == 0) {
	  r1 = (r >> 1) * 3;
	  r2 = r1 + 1;
	} else {
	  r1 = (r >> 1) * 3 + 2;
	  r2 = r1 - 1;
	}      
	if (c % 2 == 0) {
	  c1 = (c >> 1) * 3;
	  c2 = c1 + 1;
	} else {
	  c1 = (c >> 1) * 3 + 2;
	  c2 = c1 - 1;
	}      
	newe[r][c] = 0.5625 * im[r1][c1] + 0.1875 * (im[r2][c1] + im[r1][c2]) +
	  0.0625 * im[r2][c2];
      }
    return newimage;
}


/* Set result to the value of 10 * (im1 - im2).  The result is multiplied
   by 10 to provide an image that can easily be displayed and visualized.
*/
void SubtractIm10(Image im1, Image im2, Image result)
{
    float **pix1, **pix2, **pixr;
    int r, c;

    pix1 = im1->pixels;
    pix2 = im2->pixels;
    pixr = result->pixels;

    for (r = 0; r < im1->rows; r++)
	for (c = 0; c < im1->cols; c++)
	    pixr[r][c] = 10.0 * (pix1[r][c] - pix2[r][c]);
}


/* Convolve image with the a 1-D Gaussian kernel vector along image rows.
   The Gaussian has a sigma of sqrt(2), which results in the following kernel:
     (.030 .105 .222 .286 .222 .105 .030)
   Pixels outside the image are set to the value of the closest image pixel.
*/
void HorConvSqrt2(Image image, Image result)
{
    int rows, cols, r, c;
    float **pix, **rpix, *prc, p1;

    rows = image->rows;
    cols = image->cols;
    pix = image->pixels;
    rpix = result->pixels;

    for (r = 0; r < rows; r++) {
      /* Handle easiest case of pixels that do not overlap the boundary. */
      for (c = 3; c < cols - 3; c++) {
	prc =  pix[r] + c;
	rpix[r][c] = 0.030 * prc[-3] + 0.105 * prc[-2] + 0.222 * prc[-1] +
	  0.286 * prc[0] + 0.222 * prc[1] + 0.105 * prc[2] + 0.030 * prc[3];
      }
      /* For pixels near boundary, use value of boundary pixel. */
      for (c = 0; c < 3; c++) {
	p1 = (c < 1) ? pix[r][0] : pix[r][c-1];
	prc =  pix[r] + c;
	rpix[r][c] = 0.135 * pix[r][0] + 0.222 * p1 +
	  0.286 * prc[0] + 0.222 * prc[1] + 0.105 * prc[2] + 0.030 * prc[3];
      }
      for (c = cols - 3; c < cols; c++) {
	p1 = (c >= cols - 1) ? pix[r][cols-1] : pix[r][c+1];
	prc =  pix[r] + c;
	rpix[r][c] = 0.030 * prc[-3] + 0.105 * prc[-2] + 0.222 * prc[-1] +
	  0.286 * prc[0] + 0.222 * p1 + 0.135 * pix[r][cols-1];
      }
    }
}

/* Same as HorConvSqrt2, but convolve along vertical direction.
*/
void VerConvSqrt2(Image image, Image result)
{
    int rows, cols, r, c;
    float **pix, **rpix, p1;

    rows = image->rows;
    cols = image->cols;
    pix = image->pixels;
    rpix = result->pixels;

    for (c = 0; c < cols; c++) {
      /* Handle easiest case of pixels that do not overlap the boundary. */
      for (r = 3; r < rows - 3; r++) {
	rpix[r][c] = 0.030 * pix[r-3][c] + 0.105 * pix[r-2][c] +
	  0.222 * pix[r-1][c] + 0.286 * pix[r][c] + 0.222 * pix[r+1][c] +
	  0.105 * pix[r+2][c] + 0.030 * pix[r+3][c];
      }
      /* For pixels near boundary, use value of boundary pixel. */
      for (r = 0; r < 3; r++) {
	p1 = (r < 1) ? pix[0][c] : pix[r-1][c];
	rpix[r][c] = 0.135 * pix[0][c] + 0.222 * p1 +
	  0.286 * pix[r][c] + 0.222 * pix[r+1][c] +
	  0.105 * pix[r+2][c] + 0.030 * pix[r+3][c];
      }
      for (r = rows - 3; r < rows; r++) {
	p1 = (r >= rows - 1) ? pix[rows-1][c] : pix[r+1][c];
	rpix[r][c] = 0.030 * pix[r-3][c] + 0.105 * pix[r-2][c] +
	  0.222 * pix[r-1][c] + 0.286 * pix[r][c] + 0.222 * p1 +
	  0.135 * pix[rows-1][c];
      }
    }
}


/* Given a smoothed image "im", return edge gradients and orientations
      in "im" and "ori".  The gradient is computed from pixel differences,
      so is offset by half a pixel from original image.
   Result is normalized so that threshold value is 1.0, for ease of
      display and to make results less sensitive to change in threshold.
*/
void GradOriImages(Image im, Image ori)
{
    float xgrad, ygrad, invthresh;
    float **pix;
    int rows, cols, r, c;

    rows = im->rows;
    cols = im->cols;
    pix = im->pixels;

    for (r = 0; r < rows-1; r++)
      for (c = 0; c < cols-1; c++) {
	xgrad = pix[r][c+1] - pix[r][c];
	ygrad = pix[r][c] - pix[r+1][c];
	pix[r][c] = sqrt(xgrad * xgrad + ygrad * ygrad);
	ori->pixels[r][c] = atan2 (ygrad, xgrad);
      }
    /* Threshold all edge magnitudes at GradThresh and scale to 1.0. */
    invthresh = 1.0 / GradThresh;
    for (r = 0; r < rows; r++)
      for (c = 0; c < cols; c++) {
	if (pix[r][c] > GradThresh)
	  pix[r][c] = 1.0;
	else pix[r][c] *= invthresh;
      }
}
    

/* Find the local maxima and minima of the DOG images in the dog pyramid.
*/
void FindMaxMin(Image *gauss, Image *dog, Image *ori, int levels)
{
    int i, r, c, rows, cols;
    float val;

    fprintf(stderr, "Searching %d image levels, from size %dx%d\n",
	    levels, dog[0]->rows, dog[0]->cols);

    /* Search through each scale, leaving 1 scale below and 2 above. */
    for (i = 1; i < levels - 2; i++) {
        rows = gauss[i]->rows;
	cols = gauss[i]->cols;

        /* Only find peaks at least 5 pixels from borders. */
	for (r = 5; r < rows - 5; r++)
	    for (c = 5; c < cols - 5; c++) {
		val = dog[i]->pixels[r][c];
		if (fabs(val) > PeakThresh  &&
		    LocalMaxMin(val, dog[i], r, c, 0) &&
		    LocalMaxMin(val, dog[i-1], r, c, -1) &&
		    LocalMaxMin(val, dog[i+1], r, c, 1) &&
		    CheckDip(val, dog[i], r, c))
		    InterpKeyPoint(i, dog[i], gauss[i], ori[i],
				 gauss[i + 2], ori[i + 2], r, c, val);
	    }
    }
}

/* Check whether this pixel is a local maximum (positive value) or
   minimum (negative value) compared to the 3x3 neighbourhood that
   is centered at (row,col).  If level is -1 (or 1), then (row,col) must
   be scaled down (or up) a level by factor 1.5.
*/
int LocalMaxMin(float val, Image image, int row, int col, int level)
{
    int r, c;
    float **pix;

    pix = image->pixels;

    /* Calculate coordinates for image one level down (or up) in pyramid. */
    if (level < 0) {
      row = (3 * row + 1) / 2;
      col = (3 * col + 1) / 2;
    } else if (level > 0) {
      row = (2 * row) / 3;
      col = (2 * col) / 3;
    }

    if (val > 0.0) {
      for (r = row - 1; r <= row + 1; r++)
	for (c = col - 1; c <= col + 1; c++)
	  if (pix[r][c] > val)
	    return FALSE;
    } else {
      for (r = row - 1; r <= row + 1; r++)
	for (c = col - 1; c <= col + 1; c++)
	  if (pix[r][c] < val)
	    return FALSE;
    }
    return TRUE;
}


/* Check whether 20 points sampled on a ring of radius 3 pixels around
   the center point are all less than Dip * "val".
*/
int CheckDip(float val, Image image, int row, int col)
{
    float x, y, nx, ny, cosine, sine;
    int i;

    /* Create a vector of length 4 pixels (2*scale was used in orig). */
    x = 3.0;
    y = 0.0;

    /* Rotate this vector around a circle in increments of PI/20 and
       check the dip of each point. */
    sine = sin(PI / 20.0);
    cosine = cos(PI / 20.0);
    for (i = 0; i < 40; i++) {
	if (! CheckDipPoint(Dip * val, image, row + 0.5 + y, col + 0.5 + x))
	    return FALSE;
	nx = cosine * x - sine * y;
	ny = sine * x + cosine * y;
	x = nx;
	y = ny;
    }
    return TRUE;
}


/* Check whether the image point closest to (row,col) is below the
   given value.  It would be possible to use bilinear interpolation to
   sample the image more accurately, but this is not necessary as we
   are looking for a worst-case peak value.
*/
int CheckDipPoint(float val, Image image, float row, float col)
{
    int r, c;
    
    r = (int) (row + 0.5);
    c = (int) (col + 0.5);
    if (r < 0  ||  c < 0  ||  r >= image->rows  ||
        	c >= image->cols)
	return TRUE;
    if (val > 0.0) {
	if (image->pixels[r][c] > val)
	    return FALSE;
    } else {
	if (image->pixels[r][c] < val)
	    return FALSE;
    }
    return TRUE;
}


/* Find the position of the key point by interpolation in position, and
   create key vector.  Output its location as a square.
     The (grad2,ori2) provide gradient and orientation sampled one octave
   above that of (grad,ori).
*/
void InterpKeyPoint(int level, Image dog, Image grad,
		    Image ori, Image grad2, Image ori2, int r, int c, float val)
{
    float center, scale, rval, cval, opeak;

    center = dog->pixels[r][c];
    /* Scale relative to input image.  Represents size of Gaussian
       for smaller of DOG filters used at this scale. */
    scale = 1.414 * pow(1.5, (double) level);

    rval = r + InterpPeak(dog->pixels[r-1][c],
			 center, dog->pixels[r+1][c]);
    cval = c + InterpPeak(dog->pixels[r][c-1],
			 center, dog->pixels[r][c+1]);

    /* Scale up center location to account for level's scale change.
       A peak at (0,0) in gradient image would be at (0.5,0,5) in
       original image.  Add another 0.5 to place at pixel center.
       Therefore, peak at (0,0) gives location on boundary from 0 to 1. */
    rval = (rval + 1.0) * scale / 1.414;
    cval = (cval + 1.0) * scale / 1.414;

    /* If image was doubled to find keypoints, then row,col must be reduced
       by factor of 2 to put in coordinates of original image. */
    if (DoubleImSize) {
      rval *= 0.5;
      cval *= 0.5;
      scale *= 0.5;
    }

    /* Find orientation(s) for this keypoint.*/
    opeak = FindOri(grad, ori, r, c);

    MakeKeypoint(scale, opeak, r, c, rval, cval, grad, ori, grad2, ori2, val);
}


/* Return a number in the range [-0.5, 0.5] that represents the location
   of the peak of a parabola passing through the 3 samples.  The center
   value is assumed to be greater than or equal to the other values if
   positive, or less than if negative.
*/
float InterpPeak(float a, float b, float c)
{
    if (b < 0.0) {
	a = -a; b = -b; c = -c;
    }
    assert(b >= a  &&  b >= c);
    return 0.5 * (a - c) / (a - 2.0 * b + c);
}


/* Create a Keypoint record and store in global Keypoints list.
*/
void MakeKeypoint(float scale, float ori, float r, float c, float row,
	  float col, Image grad, Image orim, Image grad2, Image orim2, float val)
{
    int pool;
    Keypoint k;

    /* Store in model memory pool if we are processing training images.
       Otherwise, key only needs to be kept for the current image. */
    pool = IMAGE_POOL;
    /* fprintf(stderr, "%f  %f %d %f %f\n", row, col,10, scale,ori);*/
    if(fd)fprintf(fd, "%f  %f %f %f %f 1.0 1.0 1.0 1.0 1\n", col, row,10.0, scale,ori);
    
    k = NEW(KeypointSt, pool);
    k->next = Keypoints;
    Keypoints = k;
    k->row = row;
    k->col = col;
    k->scale = scale;
    k->ori = ori;
    k->extr = val;
    KeyCount++;
}


/* Display a square marking this location.  The size of the square
   is made 10 times the sigma of the equivalent small DOG filter, which
   gives a size 10 * 1.414 = 14 pixels in pixel sampling for this level.
   Key sampling region is size 16 pixels, for image at level of smaller
   sigma blur. 
     The square is oriented in direction "ori" with a tick mark to show the
   top.  The "smalltick" flag gives a short tick for top-down matches.
*/


/* Assign an orientation to this keypoint.  This is done by
     creating a Gaussian weighted histogram of the gradient directions in
     the region.  The histogram is smoothed and the largest peak selected.
     The results are in the range of -PI to PI.
*/
float FindOri(Image grad, Image ori, int row, int col)
{
    int i, r, c, rows, cols, radius, bin;
    float hist[OriBins], distsq, gval, weight, angle;

    rows = grad->rows;
    cols = grad->cols;

    for (i = 0; i < OriBins; i++)
      hist[i] = 0.0;

    /* Look at pixels within 3 sigma around the point and put their
       Gaussian weighted values in the histogram. */
    radius = (int) (OriSigma * 3.0);
    for (r = row - radius; r <= row + radius; r++)
      for (c = col - radius; c <= col + radius; c++)
	/* Do not use last row or column, which are not valid. */
	if (r >= 0 && c >= 0 && r < rows - 2 && c < cols - 2) {
	  gval = grad->pixels[r][c];
	  distsq = (r - row) * (r - row) + (c - col) * (c - col);
	  if (gval > 0.0  &&  distsq < radius * radius + 0.5) {
	    weight = exp(- distsq / (2.0 * OriSigma * OriSigma));
	    /* Ori is in range of -PI to PI. */
	    angle = ori->pixels[r][c];
	    bin = (int) (OriBins * (angle + PI + 0.001) / (2.0 * PI));
	    assert(bin >= 0 && bin <= OriBins);
	    bin = MIN(bin, OriBins - 1);
	    hist[bin] += weight * gval;
	  }
	}
    /* Apply smoothing twice. */
    SmoothHistogram(hist, OriBins);
    SmoothHistogram(hist, OriBins);

    return FindOriPeaks(hist, OriBins);
}


/* Smooth a histogram by using the [0.25, 0.5, 0.25] kernel.  Assume
   the histogram is connected in a circular buffer.
*/
void SmoothHistogram(float *hist, int bins)
{
    int i;
    float prev, temp;

    prev = hist[bins - 1];
    for (i = 0; i < bins; i++) {
      temp = hist[i];
      hist[i] = 0.25 * prev + 0.5 * hist[i] + 0.25 *
	          hist[(i + 1 == bins) ? 0 : i + 1];
      prev = temp;
    }
}


/* Find a peak in the histogram and return corresponding angle.
*/
float FindOriPeaks(float *hist, int bins)
{
    int i, maxloc = 0;
    float maxval = 0.0;

    /* Find peak in histogram. */
    for (i = 0; i < bins; i++)
      if (hist[i] > maxval) {
	maxval = hist[i];
	maxloc = i;
      }
    /* Set angle in range -PI to PI. */
    return (2.0 * PI * (maxloc + 0.5) / bins - PI);
}

