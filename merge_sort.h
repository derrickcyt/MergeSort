#include<xmmintrin.h>

void registerSort(__m128 *m, int n);
void bitonicSort(__m128 *a, __m128 *b);
void merge2sq(float *x, float *y, float *r, int n);
void merge2sq(float *x, float *y, float *r, int n1, int n2);
void blockInnerParallelMerge(float *a, float *b, float *r, int n, int threadNum);
void blockMergeSort(float *a, int n);
void commonMerge(float *a, float *b, float *r, int n1, int n2);