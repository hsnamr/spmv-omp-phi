//
//  sparse.c
//  
//  Created on September 2012
//

#include "sparse.h"

/******* CCS functions *******/

CCS* CreateCCS(MATRIX *m) {
    CCS *c = malloc(sizeof(CCS));
    int i, j;
    CRS *r = CreateCRS(m);
    
    int *colInd = malloc(r->nnz * sizeof(int));
    
    for(i=0; i < r->nnz; ++i)
		colInd[i] = r->colInd[i];
    
    c->nrows = r->nrows;
	c->ncols = r->ncols;
	c->nnz = r->nnz;
    
	c->A = malloc(c->nnz * sizeof(double));
	c->colPtr = malloc((c->ncols + 1) * sizeof(int));
	c->rowInd = malloc(c->nnz * sizeof(int));
	
    for (i=0; i <= c->ncols; ++i)
		c->colPtr[i] = 0;
	
    for (i=0; i < c->nnz; ++i)
		c->colPtr[r->colInd[i] + 1] += 1;
	
    for (i=0; i < c->ncols; ++i)
		c->colPtr[i+1] += c->colPtr[i];
	
    for (i=0; i < c->nrows; ++i)
		for (j = r->rowPtr[i]; j < r->rowPtr[i+1]; ++j)
			c->rowInd[j] = i;
	
    for (i=0; i < c->nnz; ++i)
		c->A[i] = r->A[i];
    
    QuickSort(c->A, c->rowInd, colInd, 0, c->nnz);
	
    DestroyCRS(r);
    free(colInd);
    
    return c;
}


void PrintCCS(CCS *c) {
    int i;
	printf("Header\n");
	printf("%u ", (int)c->nrows);
	printf("%u ", (int)c->ncols);
	printf("%u\n", (int)c->nnz);
	printf("Values - Row\n");
	
    for (i=0; i < c->nnz; ++i) {
		printf("%f\t", (double)c->A[i]);
		printf("%u\n", (int)c->rowInd[i]);
	}
	
    printf("Column Pointer (colPtr)\n");
    
    for (i=0; i < c->ncols + 1; ++i)
        printf("%u ", (int)c->colPtr[i]);
    
    printf("\n");
}


void MultiplyCCS(CCS *c, double *v, double* r) {
	int          i, j;
    double        temp;
	int*  rowInd = c->rowInd;
    int nrows = c->nrows;
	int*  colPtr = c->colPtr;
	int ncols = c->ncols;
    double *A = c->A;
    int nnz = c->nnz;
    

    memset(r, (int)0, ncols*sizeof(double));
    
    #pragma omp parallel for default(shared) private(i, j, temp)
    for (i=0; i < ncols; ++i) {
        temp = 0.0;
        #pragma ivdep
        for (j = colPtr[i]; j < colPtr[i+1]; ++j) {
            r[rowInd[j]] += (double)(A[j] * v[i]);
              temp += (double)(A[j] * v[i]);
        }
        r[rowInd[j]] = temp;
    }
}


void DestroyCCS(CCS *c) {
	free(c->colPtr);
	free(c->rowInd);
	free(c->A);
	free(c);
}


/******* CRS functions *******/

CRS* CreateCRS(MATRIX *m) {
	int          i, k, j, index;
    int nrows = m->nrows;
    int ncols = m->ncols;
    int          nnz = m->nnz;
    double** mal = m->mel;
    
    CRS *cc = (CRS*)malloc(sizeof(CRS));
    int*  colInd = malloc(nnz * sizeof(int));
	int*  rowPtr = malloc((nrows+1) * sizeof(int));
    double*  A = malloc(nnz * sizeof(double));
        
    for(i=0; i<nrows+1; ++i) {
        rowPtr[i] = 0;
    }
    
    for (i=0; i<nnz; ++i) {
        A[i] = 0;
        colInd[i] = 0;
    }
    
    index = 0;

    for(k=0; k<nrows; k++){
        for(j=0; j<ncols; j++){
            if(mal[k][j] != 0){
                A[index] = mal[k][j];
                colInd[index] = j;
                index++;
            }
        }
        rowPtr[k+1] = index;
    }
        
    cc->colInd = colInd;
    cc->rowPtr = rowPtr;
    cc->A = A;
    cc->nrows = nrows;
    cc->ncols = ncols;
    cc->nnz = nnz;

	return cc;
}


void PrintCRS(CRS *c) {
    int i;
	printf("Header\n");
	printf("%u ", (int)c->nrows);
	printf("%u ", (int)c->ncols);
	printf("%u\n", (int)c->nnz);
	printf("Values - Column\n");
	
    for (i=0; i < c->nnz ; ++i) {
		printf("%f\t", (double)c->A[i]);
		printf("%u\n", (int)c->colInd[i]);
	}
	
    printf("Row Pointer\n");
    
    for (i=0; i < c->nrows + 1; ++i)
        printf("%u ", (int)c->rowPtr[i]);
    
    printf("\n");
}


void MultiplyCRS(CRS *c, double *v, double* r) {
	int          i, j;
    double        temp;
	int *rowPtr = c->rowPtr;
	int *colInd = c->colInd;
    int nrows = c->nrows;
    int ncols = c->ncols;
    double *A = c->A;
    int nnz = c->nnz;
    double t;
    
    memset(r, (int)0, ncols*sizeof(double));
    
    #pragma omp parallel for default(shared) private(i, j, t)
    for (i=0; i < nrows; ++i) {
        t = 0.0;
        #pragma ivdep
        for (j = rowPtr[i]; j < rowPtr[i+1]; ++j) {
            t += (double)(A[j] * v[colInd[j]]);
        }
        r[i] = t;
    }
}


void DestroyCRS(CRS *c) {
	
	free(c->colInd);
	free(c->rowPtr);
	free(c->A);
	free(c);
}
