/* Utility routines for handling matrices, sorting, etc. */

#include "sift.h"


void writeImage(Image im, char *name){
  DARY *image=new DARY(im->rows,im->cols,UCHAR1);
  for (int r = 0; r < im->rows; r++)
    for (int c = 0; c < im->cols; c++) {
      image->bel[r][c] = (uchar)im->pixels[r][c];
    }
  image->write(name);

}





/*----------------------- Error messages --------------------------------*/
 
/* Print message and quit. */
void FatalError(char *msg)
{
    fprintf(stderr, "FATAL ERROR: %s\n",msg);
    abort();
}
 


/*-------------------- Pooled storage allocator ---------------------------*/

/* The following routines allow for the efficient allocation of storage in
     small chunks from a named pool.  Rather than allowing each structure to be
     freed individually, an entire pool of storage is freed at once.
   This method has two advantages over just using malloc() and free().  First,
     it is far more efficient for allocating small objects, as there is
     no overhead for remembering all the information needed to free each
     object.  Second, the decision about how long to keep an object is made
     at the time of allocation, and there is no need to track down all the
     objects to free them.
*/

/* We maintain memory alignment to word boundaries by requiring that all
   allocations be in multiples of the machine wordsize.  
*/
#define WORDSIZE 4   /* Size of machine word in bytes.  Must be power of 2. */
#define BLOCKSIZE 2048	/* Minimum number of bytes requested at a time from
			   the system.  Must be multiple of WORDSIZE. */

/* Pointers to base of current block for each storage pool (C automatically
   initializes them to NULL). */
static char *PoolBase[POOLNUM];

/* Number of bytes left in current block for each storage pool (initialized
   to 0). */
static int PoolRemain[POOLNUM];


/* Returns a pointer to a piece of new memory of the given size in bytes
   allocated from a named pool. 
*/
void *MallocPool(int size, int pool)
{
    char *m, **prev;
    int bsize;

    /* Round size up to a multiple of wordsize.  The following expression
       only works for WORDSIZE that is a power of 2, by masking last bits of
       incremented size to zero. */
    size = (size + WORDSIZE - 1) & ~(WORDSIZE - 1);

    /* Check whether new block must be allocated.  Note that first word of
       block is reserved for pointer to previous block. */
    if (size > PoolRemain[pool]) {
	bsize = (size + sizeof(char **) > BLOCKSIZE) ?
	           size + sizeof(char **) : BLOCKSIZE;
	m = (char*)malloc(bsize);
	if (! m)
	    FatalError("Failed to allocate memory.");
	PoolRemain[pool] = bsize - sizeof(void *);
	/* Fill first word of new block with pointer to previous block. */
	prev = (char **) m;
	prev[0] = PoolBase[pool];
	PoolBase[pool] = m;
    }
    /* Allocate new storage from end of the block. */
    PoolRemain[pool] -= size;
    return (PoolBase[pool] + sizeof(char **) + PoolRemain[pool]);
}

/* Free all storage that was previously allocated with MallocPool from
   a particular named pool. 
*/
void FreeStoragePool(int pool)
{
    char *prev;

    while (PoolBase[pool] != NULL) {
	prev = *((char **) PoolBase[pool]);  /* Get pointer to prev block. */
	free(PoolBase[pool]);
	PoolBase[pool] = prev;
    }
    PoolRemain[pool] = 0;
}


/*----------------- 2D matrix and image allocation ------------------------*/

/* Allocate memory for a 2D float matrix of size [row,col].  This returns
     a vector of pointers to the rows of the matrix, so that routines
     can operate on this without knowing the dimensions.
*/
float **AllocMatrix(int rows, int cols, int pool)
{
    int i;
    float **m, *v;

    m = (float **) MallocPool(rows * sizeof(float *), pool);
    v = (float *) MallocPool(rows * cols * sizeof(float), pool);
    for (i = 0; i < rows; i++) {
	m[i] = v;
	v += cols;
    }
    return (m);
}


/* Create a new image with uninitialized pixel values.
*/
Image CreateImage(int rows, int cols, int pool)
{
    Image im;

    im = NEW(ImageSt, pool);
    im->rows = rows;
    im->cols = cols;
    im->pixels = AllocMatrix(rows, cols, pool);
    return im;
}


/*--------------------- Least-squares solutions ---------------------------*/

/* Give a least-squares solution to the system of linear equations given in
   the jacobian and errvec arrays.  Return result in solution.
   This uses the method of solving the normal equations.
*/
void SolveLeastSquares(float *solution, int rows, int cols, float **jacobian,
		       float *errvec, float **sqarray)
{
    int r, c, i;
    float sum;

    assert(rows >= cols);
    /* Multiply Jacobian transpose by Jacobian, and put result in sqarray. */
    for (r = 0; r < cols; r++)
	for (c = 0; c < cols; c++) {
	    sum = 0.0;
	    for (i = 0; i < rows; i++)
		sum += jacobian[i][r] * jacobian[i][c];
	    sqarray[r][c] = sum;
	}
    /* Multiply transpose of Jacobian by errvec, and put result in solution. */
    for (c = 0; c < cols; c++) {
	sum = 0.0;
	for (i = 0; i < rows; i++)
	    sum += jacobian[i][c] * errvec[i];
	solution[c] = sum;
    }
    /* Now, solve square system of equations. */
    SolveLinearSystem(solution, sqarray, cols);
}
  

/* Solve the square system of linear equations given in sq and solution.
   Leave result in solution.   Uses Gaussian elimination with pivoting.
*/
void SolveLinearSystem(float *solution, float **sq, int size)
{
    int row, col, c, pivot = 0, i;
    float maxc, coef, temp, mult, val;

    /* Triangularize the matrix. */
    for (col = 0; col < size - 1; col++) {
	/* Pivot row with largest coefficient to top. */
	maxc = -1.0;
	for (row = col; row < size; row++) {
	    coef = sq[row][col];
	    coef = ABS(coef);
	    if (coef > maxc) {
		maxc = coef;
		pivot = row;
	    }
	}
	if (pivot != col) {
	    /* Exchange "pivot" with "col" row (this is no less efficient
	       than having to perform all array accesses indirectly). */
	    for (i = 0; i < size; i++) {
		temp = sq[pivot][i];
		sq[pivot][i] = sq[col][i];
		sq[col][i] = temp;
	    }
	    temp = solution[pivot];
	    solution[pivot] = solution[col];
	    solution[col] = temp;
	}
	/* Do reduction for this column. */
	for (row = col + 1; row < size; row++) {
	    mult = sq[row][col] / sq[col][col];
	    for (c = col; c < size; c++)	/* Could start with c=col+1. */
		sq[row][c] -= mult * sq[col][c];
	    solution[row] -= mult * solution[col];
	}
    }

    /* Do back substitution.  Pivoting does not affect solution order. */
    for (row = size - 1; row >= 0; row--) {
	val = solution[row];
	for (col = size - 1; col > row; col--)
	    val -= solution[col] * sq[row][col];
	solution[row] = val / sq[row][row];
    }
}

