/** @file Header of a minimalistic library for matrix computation.
 *  The objective of this kind of implementation is to provided
 *  elementary algorithms with the minimum overhead to be integrated
 *  in very small systems where resources are limited and dependencies
 *  must be limited. */
#ifndef __RAWMATRIX_H
#define __RAWMATRIX_H

/** Enumeration of basic status codes for matrix. */
enum {
    RAWMATRIX_OK                      =  0, /**< Matrix valid. */
    RAWMATRIX_ERR_NULL_HANDLE         =  1, /**< Null matrix pointer. */
    RAWMATRIX_ERR_INVALID_ROWS        =  2, /**< Invalid number of rows. */
    RAWMATRIX_ERR_INVALID_COLS        =  3, /**< Invalid number of columns. */
    RAWMATRIX_ERR_NO_DATA             =  4, /**< The pointer to the matrix data is null. */
    RAWMATRIX_ERR_SIZE_MISMATCH       =  5, /**< Size of two matrices are inconsistent. */
    RAWMATRIX_ERR_NOT_SQUARED         =  6, /**< The matrix should be square (n=m) but is not. */
    RAWLATRIX_ERR_NULL_INPUT          =  7, /**< An intermediate parameter is null. */
    RAWMATRIX_ERR_NOT_A_VECTOR        =  8, /**< A vector should be provided (m should be 1, but is not). */
    RAWMATRIX_MAX                     =  9  /**< Maximum value of matrix enumeration. */
};

/** Minimalistic structure for matrix. */
typedef struct struct_rawmatrix
{
    unsigned int n; /**< Number of rows. */
    unsigned int m; /**< Number of columns. */
    double *x;      /**< Pointer to the first element of the matrix.
                     *   The data of the matrix are handled as a (n x m) x array of doubles.*/
} rawmatrix_t;



/** Consistency checks for a matrix. */
extern int rawmatrix_check(const rawmatrix *const A); /**< Pointer to the matrix to be checked. */

/** LU decomposition used to solve a linear system. */
extern int rawmatrix_lu_decomposition(rawmatrix_t *A,     /**< Matrix to be decomposed. */
                                      unsigned int *idx,  /**< Vector of permutation indexes. */
                                      double *s,          /**< Vector with scale factors. */
                                      int *multiplier,    /**< +1 if the number of permuted rows is even, -1 if odd. */
                                      int *isnan);        /**< Returns the index of current row when singularity has been detected. */


/** Resolution of linear system using LU decomposition. */
extern int rawmatrix_lu_resolution(rawmatrix_t *LU,
                                   rawmatrix_t *B,
                                   unsigned int *idx);


#endif
