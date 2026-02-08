//
//  blocked.c — BCRS (block compressed row) and ELL (ELLPACK) formats and SpMV.
//  Created on April 2013
//

#include "blocked.h"

/******* BCRS (Block Compressed Row Storage) *******/

/* Count contiguous column blocks of nonzeros per row (used for BCRS layout). */
int countBlocks(MATRIX* my_m) {
    int i, j, nblocks = 0;
    int nrows = my_m->nrows;
    int ncols = my_m->ncols;
    double** mal = my_m->mel;
    
    for (i = 0; i < nrows; ++i) {
        for (j = 0; j < ncols; ++j) {
            if (mal[i][j] != 0) {
                // is this the first element in the block?
                if (j == 0) {
                    nblocks++;
                } else if ((j != 0) && (mal[i][j-1] == 0)) {
                    nblocks++;
                }
            }
        }
    }
    return nblocks;
}

/* Build BCRS from dense matrix: rowPtr, colInd, nnzPtr, value. */
BCRS* CreateBCRS(MATRIX* my_m) {
    int i, jj, k, index, block_number = 0;
    int nrows = my_m->nrows;
    int ncols = my_m->ncols;
    int nnz = my_m->nnz;
    double** mal = my_m->mel;
    double current;
    
    BCRS* cc = (BCRS*)malloc(sizeof(BCRS));
    
    int nblocks = countBlocks(my_m);
 
    int*  colInd = malloc(nblocks * sizeof(int));
    int*  rowPtr = malloc((nrows+1) * sizeof(int));
    int*  nnzPtr = malloc((nblocks+1) * sizeof(int));
    double*  value = malloc(nnz * sizeof(double));
     
    for(i=0; i<nrows+1; ++i) {
         rowPtr[i] = 0;
    }

    index = 0;
    
    for(k=0; k<nrows; k++){
        for(jj=0; jj<ncols; ++jj){
            if (mal[k][jj] != 0) {
                if ((jj == 0) || (jj != 0 &&(mal[k][jj-1] == 0))) {
                    colInd[block_number] = jj;
                    nnzPtr[block_number] = index;
                    block_number++;
                }
                value[index] = mal[k][jj];
                index++;
            }
        }
        rowPtr[k+1] = block_number;
    }
    nnzPtr[nblocks] = nnz;
    
    cc->colInd = colInd;
    cc->rowPtr = rowPtr;
    cc->value = value;
    cc->nnzPtr = nnzPtr;
    cc->nrows = nrows;
    cc->ncols = ncols;
    cc->nnz = nnz;
    cc->nblocks = nblocks;
    
    return cc;
}

void PrintBCRS(BCRS *c) {
    printf("To-Do: Not implemented");
}

/* SpMV for BCRS: r = A * x. r length nrows, x length ncols. */
void MultiplyBCRS(BCRS *c, double *x, double *r) {
    int          i, j, s, startCol, k, tt;
    double        temp;
	int *rowPtr = c->rowPtr;
    int *nnzPtr = c->nnzPtr;
	int *colInd = c->colInd;
    int nrows = c->nrows;
    int ncols = c->ncols;
    double *value = c->value;
    int nnz = c->nnz;
    double t;
    

    memset(r, (int)0, ncols*sizeof(double));

      #pragma omp parallel for default(shared) private(i, j, k, startCol, t, tt)
      for (i = 0; i < nrows; ++i) {
        t = 0.0;
        for (j = rowPtr[i]; j < rowPtr[i+1]; ++j) {
            startCol = colInd[j];
            tt = 0;
            #pragma ivdep
            for (k = nnzPtr[j]; k < nnzPtr[j+1]; ++k) {
                t += value[k]*x[startCol+tt];
                tt++;
            }
        }
        r[i] = t;
    }
}

void DestroyBCRS(BCRS *cc) {
    free(cc->colInd);
    free(cc->rowPtr);
    free(cc->value);
    free(cc->nnzPtr);
    free(cc);
}


/******* ELL (ELLPACK) *******/

/* Maximum number of nonzeros in any row (determines ELL column width). */
int compute_max_entries_per_row(MATRIX* m)
{
    int max_entries_per_row = 0;
    int current = 0;
    int i, j;
    
    for (i = 0; i < m->nrows; ++i) {
        current = 0;
        for (j = 0; j < m->ncols; ++j) {
            if(m->mel[i][j] != 0) {
                current++;
            }
        }
        max_entries_per_row = max(current, max_entries_per_row);
    }
    
    return max_entries_per_row;
}

/* Build ELL: values and indices 2D arrays; short rows padded with 0 and index -1. */
ELL* CreateELL(MATRIX *m) {
    ELL* ell = (ELL*) malloc(sizeof(ELL));
    double** values;
    int** indices;
    int nnz = m->nnz;
    int ncols = m->ncols;
    int nrows = m->nrows;
    int i, j, jj;
    double** mal = m->mel;
    
    int max_entries_per_row = compute_max_entries_per_row(m);
    
    values = (double**) malloc(m->nrows * sizeof(double*));

    for (i = 0; i < m->nrows; ++i) {
        values[i] = (double*) malloc(max_entries_per_row * sizeof(double));
    }
    
    indices = (int**) malloc(m->nrows * sizeof(int*));
    for (i = 0; i < m->nrows; ++i) {
        indices[i] = (int*) malloc(max_entries_per_row * sizeof(int));
    }
    
    for (i = 0; i < nrows; ++i) {
        for (jj = 0; jj < max_entries_per_row; ++jj) {
            values[i][jj] = 0;
            indices[i][jj] = -1;
        }
    }
        
    for (i = 0; i < nrows; ++i) {
        jj = 0;
        for (j = 0; j < ncols; ++j) {
            if(mal[i][j] != 0) {
                values[i][jj] = mal[i][j];
                indices[i][jj] = j;
                jj++;               // less than or equal to max_entries_per_row
            }
        }
    }
        
    ell->values = values;
    ell->indices = indices;
    ell->nnz = nnz;
    ell->nrows = nrows;
    ell->ncols = ncols;
    ell->max_entries_per_row = max_entries_per_row;

    return ell;
}


void PrintELL(ELL *e) {
    printf("To-Do: Not implemented");
}


/* SpMV for ELL: r = A * x. Padding entries (index -1) are skipped. */
void MultiplyELL(ELL* e, double* x, double* r) {
    int i, j;
    int **indices = e->indices;
    int nrows = e->nrows;
    int ncols = e->ncols;
    int max_entries_per_row = e->max_entries_per_row;
    double **A = e->values;
    double t;
    
    memset(r, 0, (size_t)nrows * sizeof(double));
    
    #pragma omp parallel for default(shared) private(i, j, t)
    for (i = 0; i < nrows; ++i) {
        t = 0.0;
        #pragma ivdep
        for (j = 0; j < max_entries_per_row; ++j) {
            if (indices[i][j] >= 0)  /* skip padding entries */
                t += A[i][j] * x[indices[i][j]];
        }
        r[i] = t;
    }
}

void DestroyELL(ELL* c) {
    int i;
    
    for (i = 0; i < c->nrows; ++i) {
        free(c->values[i]);
        free(c->indices[i]);
    }
    free(c->values);
    free(c->indices);
    free(c);
}
