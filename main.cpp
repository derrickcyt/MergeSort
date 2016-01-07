#include<stdlib.h>
#include<time.h>
#include<iostream>
#include"merge_sort.h"

using namespace std;

void generateArray(float *a, int n)
{
	srand(time(NULL));
	for (int i = 0; i < n; i++)
		a[i] = rand()%10000/pow(10,rand()%5);
}

int cmp(const void *a, const void *b)
{
	if (*(float *)a < *(float *)b)
		return -1;
	else if (*(float *)a >*(float *)b)
		return 1;
	else return 0;
}

float calError(float *a, float *b, int n)
{
	float e = 0;
	for (int i = 0; i < n; i++)
	{
		float t = a[i] - b[i];
		if (t < 0) t = -t;
		e += t;
	}
	return e;
}


int main()
{
	//clock_t start;
	//int n = 1024 * 1024 * 100;
	///*int a[12] = { 156, 15848, 1533, 148, 74, 5, 748, 156, 158747, 1258, 17878, 7811 };*/
	//int *a = (int*)malloc(sizeof(int)* n);
	//generateArray(a, n);
	//int *b = (int*)malloc(sizeof(int)* n);
	//memcpy(b, a, sizeof(int)* n);
	//start = clock();
	//qsort(a, n, sizeof(a[0]), cmp);
	//printf("快排时间：%dms\n", clock() - start);


	int n = 1024*1024*16;
	float *a = (float*)malloc(sizeof(float)*n*2);
	generateArray(a, 2*n);
	
	float *b = (float*)_mm_malloc(sizeof(float)*n*2,16);
	float *r = (float*)_mm_malloc(sizeof(float)*n * 2,16);
	memcpy(b, a, sizeof(float)*n * 2);
	clock_t start = clock();
	qsort(a, n*2, sizeof(a[0]), cmp);
	printf("快排时间：%dms\n",clock()-start);
	start = clock();
	blockMergeSort(b, n*2);
	printf("MergeSort时间：%dms\n", clock() - start);
	printf("误差%f\n\n", calError(a, b, n * 2));
	/*qsort(b, n, sizeof(a[0]), cmp);
	qsort(&b[n], n, sizeof(a[0]), cmp);*/
	//blockInnerParallelMerge(b, &b[n], r, n, 3);
	//merge2sq(b, &b[n], r, n, n);
	//printf("误差%f\n\n", calError(a, r, n * 2));
	/*for (int i = 0; i < n * 2; i++)
		printf("%f ", a[i]);
	printf("\n\n");
	for (int i = 0; i < n * 2; i++)
		printf("%f ", r[i]);
	printf("\n");*/


	/*float a[4] = { 1, 5, 4, 6 };
	float b[4] = { 2, 3, 7, 8 };
	float c[4] = { 1.2, 3.2, 4.5, 7.6 };
	float d[4] = { 2.2, 4.5, 6.6, 7.8 };
	__m128 m[4];
	m[0] = _mm_load_ps(a);
	m[1] = _mm_load_ps(b);
	m[2] = _mm_load_ps(c);
	m[3] = _mm_load_ps(d);
	registerSort(m,4);
	_mm_store_ps(a, m[0]);
	_mm_store_ps(b, m[1]);
	_mm_store_ps(c, m[2]);
	_mm_store_ps(d, m[3]);
	for (int i = 0; i < 4; i++)
		printf("%f ", a[i]);
	printf("\n");
	for (int i = 0; i < 4; i++)
		printf("%f ", b[i]);
	printf("\n");
	for (int i = 0; i < 4; i++)
		printf("%f ", c[i]);
	printf("\n");
	for (int i = 0; i < 4; i++)
		printf("%f ", d[i]);
	printf("\n");
	printf("\n");
	m[0] = _mm_load_ps(a);
	m[1] = _mm_load_ps(b);
	bitonicSort(&m[0], &m[1]);
	_mm_store_ps(a, m[0]);
	_mm_store_ps(b, m[1]);
	for (int i = 0; i < 4; i++)
		printf("%f ", a[i]);
	printf("\n");
	for (int i = 0; i < 4; i++)
		printf("%f ", b[i]);
	printf("\n");
	m[0] = _mm_load_ps(c);
	m[1] = _mm_load_ps(d);
	bitonicSort(&m[0], &m[1]);
	_mm_store_ps(c, m[0]);
	_mm_store_ps(d, m[1]);
	for (int i = 0; i < 4; i++)
		printf("%f ", c[i]);
	printf("\n");
	for (int i = 0; i < 4; i++)
		printf("%f ", d[i]);
	printf("\n");
	printf("\n");

	float *x = (float*)_mm_malloc(sizeof(float)* 8,16);
	float *y = (float*)_mm_malloc(sizeof(float)* 8,16);
	float *r = (float*)_mm_malloc(sizeof(float)* 16, 16);
	memcpy(x,a,sizeof(float)*4);
	memcpy(&x[4], b, sizeof(float)* 4);
	memcpy(y, c, sizeof(float)* 4);
	memcpy(&y[4], d, sizeof(float)* 4);
	blockInnerParallelMerge(x, y, r, 8, 3);

	for (int i = 0; i < 16; i++)
		printf("%f ",r[i]);
	printf("\n");*/


	int in;
	cin >> in;



}