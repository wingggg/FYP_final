Package contains statically linked binary of Maximum Stable Extremal 
Regions detector.

Usage: extrema-bin [options]
Usage: mserdet.exe [options]
    -i (null) [null] input image (png, tiff, jpg, ppm, pgm)
    -o (null) [null] output file
    -es (1.000) [1.000] ellipse scale, (output types 2 and 4)
    -pre (0) [0] image preprocessing type
    -per (0.010) [0.010] maximum relative area
    -ms (30) [30] minimum size of output region
    -mm (10) [10] minimum margin
    -rel (0) [0] use relative margins
    -pe (0) [0] speedup detection by preprocessing local extrema
    -v (0) [0] verbose output
    -d (0) [0] debug outputs
    -r (1) [1] number of runs (for timing purposes)
    -t (0) [0] output file type 0 - RLE, 1 - ExtBound., 2 - Ellipse

    -help (0) [0] print out usage info
Dependencies:
    Option -i is compulsory
Errors detected during  option parsing:
    Missing compulsory option -i

Detected MSER+ and MSER- regions are stored in output file as follows:


Extrema file format RLE:

NUM_MSER_PLUS
NUM_RLE LINE COL1 COL2 LINE COL1 COL2 ... LINE COL1 COL2
...
NUM_RLE LINE COL1 COL2 LINE COL1 COL2 ... LINE COL1 COL2
NUM_MSER_MINUS
NUM_RLE LINE COL1 COL2 LINE COL1 COL2 ... LINE COL1 COL2
...
NUM_RLE LINE COL1 COL2 LINE COL1 COL2 ... LINE COL1 COL2

where NUM_MSER_PLUS and NUM_MSER_MINUS are number of MSER+ and MSER- regions 
respectively. Each region is described as one line in output file. 
NUM_RLE specifies number of RLE triples LINE COL1 COL2.


Extrema file format Extended boundary:

NUM_MSER_PLUS
NUM_PTS X Y X Y ... X Y
...
NUM_PTS X Y X Y ... X Y
NUM_MSER_MINUS
NUM_PTS X Y X Y ... X Y
...
NUM_PTS X Y X Y ... X Y

where NUM_MSER_PLUS and NUM_MSER_MINUS are number of MSER+ and MSER- regions 
respectively. Each region is described as one line in output file. NUM_PTS 
specifies number of points in extended boundary.


Extrema file format Ellipse:

1.0
NUM_TOTAL_MSER
U V A B C 
...
U V A B C

where NUM_TOTAL_MSER is number of MSER+ and MSER- regions. Each region is 
described as one line in output file. Each affine region is described as an 
ellipse with parameters U,V,A,B,C in:
 
a(x-u)(x-u)+2b(x-u)(y-v)+c(y-v)(y-v)=1

with (0,0) at image top left corner.
