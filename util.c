//
//  util.c — Timing, sorting, and small helpers used by sparse/reorder.
//  Created on September 2012
//

#include "util.h"
#include <sys/time.h>

/* Return current time in seconds (for simple wall-clock timing). */
double dtime() {
    double tseconds = 0.0;
    struct timeval mytime;
    gettimeofday(&mytime, (struct timezone*)0);
    tseconds = (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);
    return tseconds;
}

double Random() {    
    return (double)drand48();
}

/* Quicksort by key array c; permute arrays a and b accordingly (for CCS construction). */
void QuickSort(double *a, int *b, int *c, int beg, int end) {
    
	if (end > beg + 1) {
        int piv = c[beg];
		int l = beg + 1;
		int r = end;
        while (l < r) {
			if (c[l] <= piv)
                ++l;
			else {
				--r;
				swapI(&c[l], &c[r]);
				swapF(&a[l], &a[r]);
				swapI(&b[l], &b[r]);
			}
		}
		--l;
		swapI(&c[l], &c[beg]);
		swapF(&a[l], &a[beg]);
		swapI(&b[l], &b[beg]);
		QuickSort(a, b, c, beg, l);
		QuickSort(a, b, c, r, end);
	}
}

void swapF(double *a, double *b) {
	double temp = *a;
	*a = *b;
	*b = temp;
}

void swapI(int *a, int *b) {
	int temp = *a;
	*a = *b;
	*b = temp;
}

int findMax(int *input, int n) {
    int j;
    int ncols = n;
    
    int max = 0;
    int max_pos = 0;
    
    for (j = 0; j < ncols; ++j) {
        if (input[j] >= max) {
            max = input[j];
            max_pos = j;
        }
    }
    
    return max_pos;
}

int exists(int *input, int ncols, int n) {
    int j;
    
    for (j = 0; j < ncols; ++j) {
        if (input[j] == n) {
            return 1;
        }
    }
    
    return 0;
}
