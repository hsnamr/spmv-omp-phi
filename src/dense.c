//
//  dense.c — Dense MATRIX type, Matrix Market I/O, and dense SpMV.
//  Created on September 2012
//

#include "dense.h"

/* Read a sparse matrix in Matrix Market coordinate format into dense MATRIX (mel). */
MATRIX* ReadMatrix(FILE *file) {
    MM_typecode matcode;
    int         M, N, nz, ret_code, i, j, k;
    double      *val;
    double      **A;
    MATRIX      *m;
        
    if (mm_read_banner(file, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }    
    
    
    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */
    
    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
        mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Matrix Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }
    
    
    /* find out size of sparse matrix .... */
    
    if ((ret_code = mm_read_mtx_crd_size(file, &M, &N, &nz)) !=0)
        exit(1);
        
    /* reseve memory for matrices */
    m = (MATRIX *) malloc(sizeof(MATRIX));
    int*  I = (int *) malloc(nz * sizeof(int));
    int*  J = (int *) malloc(nz * sizeof(int));
    val = (double *) malloc(nz * sizeof(double));
            
    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */
    
    for (i=0; i<nz; i++)
    {
        fscanf(file, "%d %d %lf\n", &I[i], &J[i], &val[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }
    
    if (file !=stdin) fclose(file);
        
    A = (double**) malloc(M * sizeof(double*));
            
    for (i = 0; i < M; ++i) {
        A[i] = (double*) malloc(N * sizeof(double));
    }
                
    for(j=0; j < M; ++j) {
        for(k=0; k < N; ++k) {
            A[j][k] = 0.0;
        }
    }
    
    for(i=0; i < nz; ++i) {
        j = I[i];
        k = J[i];
        A[j][k] = val[i];
    }
        
    m->nnz = nz;
    m->nrows = M;
    m->ncols = N;
    m->mel = A;
    
	return m;
}

/* Write MATRIX to file in Matrix Market coordinate format. */
void WriteMatrix(FILE *file, MATRIX* reM) {
    MM_typecode matcode;
    int i, jjj, k;
    int nz = reM->nnz;
    int*  I = (int *) malloc(nz * sizeof(int));
    int*  J = (int *) malloc(nz * sizeof(int));
    double* val = (double *) malloc(nz * sizeof(double));
    int nrows = reM->nrows;
    int ncols = reM->ncols;
    int N = nrows;
    int M = ncols;

    k = 0;    
    while (k < nz) {
        for (i = 0; i < nrows; ++i) {
            for (jjj = 0; jjj < ncols; ++jjj) {
                if((double)reM->mel[i][jjj] != (double)0) {
                    I[k] = i;
                    J[k] = jjj;
                    val[k] = reM->mel[i][jjj];
                    k++;
                }
            }
        }
        break;
    } // end while
    
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_real(&matcode);
    
    mm_write_banner(file, matcode);
    mm_write_mtx_crd_size(file, M, N, nz);
    
    /* NOTE: matrix market files use 1-based indices, i.e. first element
     of a vector has index 1, not 0.  */
    
    for (i=0; i<nz; i++)
        fprintf(file, "%d %d %10.3g\n", I[i]+1, J[i]+1, val[i]);
    
    fflush(file);
    
    #ifdef LOG
    printf("Closing file...\n");
    #endif
    
    fclose(file);
    
    free(I);
    free(J);
    free(val);
}

/* Parallel dot product of vectors a and b of length n. */
double DotProduct(const double *a, const double *b, int n) {
     double result = 0.0f;
    int i;
    
    #pragma omp parallel for default(shared) private(i) reduction(+:result)
    for(i = 0; i < n; i++) {
        result += (a[i] * b[i]);
    }
    
    return result;
}


/* Dense SpMV: return r = m * val. Caller must free the returned vector. */
double* MultiplyMatrix(MATRIX *m, double *val) {
	int          i, j, s;
    int nrows = m->nrows;
    int ncols = m->ncols;
    double*  r = malloc((size_t)nrows * sizeof(double));
    double**  mal = m->mel;
    
    #pragma omp parallel for
    for (s = 0; s < nrows; s++) {
        r[s] = 0;
    }

    #pragma omp parallel for default(shared) private(j)
    for (i = 0; i < nrows; ++i) {
        double result = 0.0;
        for (j = 0; j < ncols; ++j) {
            result += (double)(mal[i][j] * val[j]);
        }
        r[i] = result;
    }
    
    return r;
}

/* Print dense matrix to stdout. */
void PrintMatrix(MATRIX *m) {
    int i,
        j;
    
    for (i = 0; i < m->nrows; i++) {
        for (j = 0; j < m->ncols; j++) {
            printf("%f ", m->mel[i][j]);
        }
        printf("\n");
    }
}

/* Free all rows and the MATRIX struct. */
void DestroyMatrix(MATRIX *m) {
    int i;
    
    for (i = 0; i < m->nrows; ++i) {
        free(m->mel[i]);
    }
    
    free(m->mel);
    free(m);
}
