//
//  reorder.c
//  
//  Created by on November 2013
//

#include <stdio.h>
#include <stdlib.h>
#include "reorder.h"

int index_gc = 0;

int P[20000] = {0};

void column_intersection(double** input, double** output, int nrows, int ncols) {
    /* input: original sparse matrix - will not be modified
     output: reordered sparse matrix */
    int i, ij, j, ijk, max_pos;
    int column_length = nrows;
    int row_length = ncols;
    
    int *intersection = (int *)malloc(ncols*sizeof(int));
    
    for (j = 0; j < ncols; ++j) {
        intersection[j] = 0;
    }
    
    ijk = 0;
    P[ijk] = 0;       // choice of 1st column effects resulting reorder
    
    while (ijk < ncols-1) {
        #pragma omp parallel for default(shared) private(j, i)
        for (j = 0; j < ncols; ++j) {
            for (i = 0; i < nrows; ++i) {
                if (input[i][j] != 0) {
                    if (input[i][P[ijk]] != 0) {
                        intersection[j]++;
                    }
                }
            }
        }
        
        max_pos = findMax(intersection, ncols);
        #pragma omp parallel for default(shared) private(ij)
        for (ij = 0; ij < ncols; ++ij) {
            /* if max_pos is in permutation find next max */
            if(exists(P, ncols, max_pos)) { /* find next max */
                intersection[max_pos] = -1;
                max_pos = findMax(intersection, ncols);
            }
        }
        
        ijk++;

        P[ijk] = max_pos;
        
        /* reset intersection weights */
        #pragma omp parallel for default(shared) private(j)
        for (j = 0; j < ncols; ++j) {
            intersection[j] = 0;
        }
    }  // while-loop done
    
    free(intersection);
    
    return;
}

void gray_code(double** A, int* C, int rowIndex, int sign, int nrows, int ncols) {
    int ijk, jj, ii, i1, i2;
    int *C1 = (int*)malloc(ncols*sizeof(int));
    // C1: int array holding indices of all columns with nonzero at rowIndex
    int *C2 = (int*)malloc(ncols*sizeof(int));
    // C2: int array holding indices of all columns with zero at rowIndex
    
    int length1 = ncols;
    int length2 = ncols;
    
    i1 = 0;
    i2 = 0;
    
    #ifdef LOG
    printf("ncols %d, nrows %d, rowIndex %d, sign %d\n", ncols, nrows, rowIndex, sign);
    #endif
    
    if (rowIndex == nrows || ncols <= 1) {
        for (jj = 0; jj < ncols; ++jj) {
            P[index_gc] = C[jj];
            index_gc++;
        }
        free(C1);
        free(C2);
        return;
    }
    
    /* reset C1 and C2 (is this necessary?) */
    for (ijk = 0; ijk < ncols; ++ijk) {
        C1[ijk] = -1;
        C2[ijk] = -1;
    }
    
    /* partition C into two disjoint sets C1 and C2 such that C1 contains all column indices that have nonzeros in rowIndex and C2 contains the other indices */
    for (ii = 0; ii < ncols; ++ii) {
        if(A[rowIndex][C[ii]] != 0) {
            C1[i1] = C[ii];
            i1++;
        } else if(A[rowIndex][C[ii]] == 0) {
            C2[i2] = C[ii];
            i2++;
        }
    }
    
    rowIndex++;
    
    // update new length (ncols on next recursive call)
    length1 = i1;
    length2 = i2;
    
    if (sign == +1) {
        gray_code(A, C1, rowIndex, -1, nrows, length1);
        gray_code(A, C2, rowIndex, +1, nrows, length2);
    } else {
        gray_code(A, C2, rowIndex, -1, nrows, length2);
        gray_code(A, C1, rowIndex, +1, nrows, length1);
    }
    
}
