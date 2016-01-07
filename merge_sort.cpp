#include<iostream>
#include<math.h>


using namespace std;

#define CACHESIZE 1024*1024*4
#define SIMDWIDTH 4
#define WORDLENGTH 16
#define PARALLEL_SIZE_LIMIT 1024*1024
#define THREAD_NUM 4
#define FIFO_SIZE 4*1024



void registerSort(__m128 *m,int n)
{
	for (int i = 0; i < n-1; i++)
	{
		for (int j = 0; j < n - i - 1; j++)
		{
			//bubble sort
			__m128 l = _mm_min_ps(m[j], m[j+1]);
			__m128 h = _mm_max_ps(m[j], m[j+1]);
			m[j] = l; m[j+1] = h;
		}
	}
	_MM_TRANSPOSE4_PS(m[0], m[1], m[2], m[3]);
}


void printM128(__m128 m)
{
	float a[4];
	_mm_store_ps(a, m);
	for (int i = 0; i < 4; i++)
		printf("%f ",a[i]);
}

void bitonicSort(__m128 *a,__m128 *b)
{
	*b = _mm_shuffle_ps(*b, *b, _MM_SHUFFLE(0,1,2,3));//reverse b

	__m128 l = _mm_min_ps(*a, *b);
	__m128 h = _mm_max_ps(*a, *b);
	__m128 l1p=_mm_shuffle_ps(l,h,_MM_SHUFFLE(1,0,1,0));
	__m128 h1p=_mm_shuffle_ps(l, h, _MM_SHUFFLE(3,2,3,2));

	l = _mm_min_ps(l1p,h1p);
	h = _mm_max_ps(l1p, h1p);
	l1p = _mm_shuffle_ps(l, h, _MM_SHUFFLE(0, 2, 0, 2));
	l1p = _mm_shuffle_ps(l1p, l1p, _MM_SHUFFLE(2,0,3,1));
	h1p = _mm_shuffle_ps(l, h, _MM_SHUFFLE(1,3,1,3));
	h1p = _mm_shuffle_ps(h1p, h1p, _MM_SHUFFLE(2,0, 3,1));

	l = _mm_min_ps(l1p, h1p);
	h = _mm_max_ps(l1p, h1p);

	*a = _mm_shuffle_ps(l,h,_MM_SHUFFLE(1,0,1,0));
	*a = _mm_shuffle_ps(*a, *a, _MM_SHUFFLE(3,1,2,0));
	*b = _mm_shuffle_ps(l, h, _MM_SHUFFLE(3,2,3,2));
	*b = _mm_shuffle_ps(*b, *b, _MM_SHUFFLE(3, 1, 2, 0));

}


void commonMerge(float *a, float *b, float *r, int n1, int n2)
{
	int num = 0, ai = 0, bi = 0;
	while (ai<n1&&bi<n2)
	{
		if (a[ai] <= b[bi]){
			r[num++] = a[ai++];
		}
		else{
			r[num++] = b[bi++];
		}
	}
	if (ai < n1){
		memcpy(&r[num], &a[ai], sizeof(float)*(n1 - ai));
	}
	else{
		memcpy(&r[num], &b[bi], sizeof(float)*(n2 - bi));
	}
}

//same length
void merge2sq(float *x,float *y,float *r,int n)
{
	float *tempF = (float*)_mm_malloc(sizeof(float)*SIMDWIDTH,WORDLENGTH);
	int xIndex = 0, yIndex = 0,rIndex=0;
	__m128 v1 = _mm_load_ps(x);
	__m128 v2 = _mm_load_ps(y);
	bitonicSort(&v1, &v2);
	_mm_store_ps(tempF, v1);
	memcpy(r, tempF, sizeof(float)*SIMDWIDTH);
	rIndex += 4;
	xIndex += 4;
	yIndex += 4;
	v1 = v2;
	while (xIndex<n&&yIndex<n)//if x or y is empty
	{
		if (x[xIndex] <= y[yIndex])
		{
			v2 = _mm_load_ps(&x[xIndex]);
			xIndex += 4;
		}
		else
		{
			v2 = _mm_load_ps(&y[yIndex]);
			yIndex += 4;
		}

		bitonicSort(&v1, &v2);
		_mm_store_ps(&r[rIndex], v1);
		rIndex += 4;
		v1 = v2;
	}
	//printf("xIndex=%d,yIndex=%d\n",xIndex,yIndex);
	_mm_store_ps(tempF, v1);
	if (n == 4){
		memcpy(&r[rIndex], tempF, sizeof(float)*SIMDWIDTH);
		return;
	}
	if (xIndex < n)
	{
		commonMerge(tempF, &x[xIndex], &r[rIndex], SIMDWIDTH, n - xIndex);
	}
	if (yIndex < n)
	{
		commonMerge(tempF, &y[yIndex], &r[rIndex], SIMDWIDTH, n -yIndex);
	}
	_mm_free(tempF);
}



//custom length
void merge2sq(float *x, float *y, float *r, int n1,int n2)
{
	float *temp = (float*)_mm_malloc(sizeof(float)* 4, 16);
	if (n1 < 4||n2<4){
		commonMerge(x, y, r, n1, n2);
		return;
	}
	int xi = 0, yi = 0, rIndex = 0;
	__m128 v1 = _mm_load_ps(x);
	__m128 v2 = _mm_load_ps(y);
	
	bitonicSort(&v1, &v2);
	_mm_store_ps(temp, v1);
	memcpy(r, temp, sizeof(float)* 4);
	rIndex += 4;
	xi += 4;
	yi += 4;
	v1 = v2;
	while (xi<n1&&yi<n2)//if x or y isn't empty
	{
		if (n1 - xi < 4 || n2 - yi < 4)
			break;
		if (x[xi] <= y[yi])
		{
			v2 = _mm_load_ps(&x[xi]);
			xi += 4;
		}
		else
		{
			v2 = _mm_load_ps(&y[yi]);
			yi += 4;
		}

		bitonicSort(&v1, &v2);
		_mm_store_ps(temp, v1);
		memcpy(&r[rIndex], temp, sizeof(float)* 4);
		//_mm_store_ps(&r[rIndex], v1);
		rIndex += 4;
		v1 = v2;
	}
		//v1 restore to temp
		_mm_store_ps(temp, v1);
		if (xi >= n1){
			commonMerge(temp, &y[yi], &r[rIndex], 4, n2 - yi);
			_mm_free(temp);
		}
		else{
			//merge x1 and temp to r1;
			float *r1 = (float*)_mm_malloc(sizeof(float)* (4 + n1 - xi), 16);
			commonMerge(temp, &x[xi], r1, 4, n1-xi);
			if (yi >= n2){
				memcpy(&r[rIndex],r1,sizeof(float)*(4+n1-xi));
			}
			else{
				commonMerge(r1, &y[yi], &r[rIndex], 4 + n1 - xi, n2 - yi);
			}
			_mm_free(temp); _mm_free(r1);
		}

}

void divide2sq(int *aDiv,int *bDiv,int *rDiv,float *a,float *b,int n,int threadNum)
{
	int bSize = ceil((double)n * 2 / threadNum);
	int divNum = 1, ai = 0, bi = 0, num = 0;
	aDiv[0] = 0; bDiv[0] = 0;
	while (divNum < threadNum)
	{
		if (ai < n&&bi < n)
		{
			if (a[ai] <= b[bi])
				ai++;
			else
				bi++;
			num++;
			if (num >= bSize)
			{
				aDiv[divNum] = ai;
				bDiv[divNum] = bi;
				divNum++;
				num = 0;
			}
		}
		else
		{
			if (ai + bSize-num >= n&&bi+bSize-num>=n)
			{
				divNum++;
			}
			if (ai < n)
			{
				aDiv[divNum] = ai + bSize-num;
				bDiv[divNum] = n;
				ai += bSize;
				divNum++;
			}
			else{
				aDiv[divNum] = n;
				bDiv[divNum] = bi + bSize - num;
				bi += bSize;
				divNum++;
			}
			num = 0;
		}
	}
	rDiv[0] = 0;
	for (int i = 1; i < threadNum; i++){
		rDiv[i] = i*bSize;
	}
}


void blockInnerParallelMerge(float *a,float *b,float *r,int n,int threadNum)
{
	if (threadNum <= 1 || n<PARALLEL_SIZE_LIMIT)
	{
		merge2sq(a, b, r, n);
		return;
	}
	int *aDiv = (int*)malloc(sizeof(int)*(threadNum));
	int *bDiv = (int*)malloc(sizeof(int)*(threadNum));//bDiv[0]!=0
	int *rDiv = (int *)malloc(sizeof(int)*(threadNum));//rDiv[0]=0;i:i*bSize
	divide2sq(aDiv,bDiv,rDiv,a,b,n,threadNum);

	//parallel
#pragma omp parallel for num_threads(threadNum) 
	for (int i = 0; i < threadNum; i++)
	{
		//get sq to handle
		int xLen = i == threadNum - 1 ? n - aDiv[i] : aDiv[i + 1] - aDiv[i];
		int yLen = i == threadNum - 1 ? n - bDiv[i] : bDiv[i + 1] - bDiv[i];
		//printf("xLen=%d,yLen=%d\n",xLen,yLen);
		if (xLen == 0)
		{
			memcpy(&r[rDiv[i]],&b[bDiv[i]],sizeof(float)*yLen);
		}
		else if (yLen == 0){
			memcpy(&r[rDiv[i]], &a[aDiv[i]], sizeof(float)*xLen);
		}
		else{
			float *x = (float*)_mm_malloc(sizeof(float)*xLen,16);
			float *y = (float*)_mm_malloc(sizeof(float)*yLen, 16);
			memcpy(x, &a[aDiv[i]], sizeof(float)*xLen);
			memcpy(y, &b[bDiv[i]], sizeof(float)*yLen);
			merge2sq(x, y, &r[rDiv[i]], xLen, yLen);
			_mm_free(x); _mm_free(y);
		}
	}
	free(aDiv); free(bDiv); free(rDiv);

}


void blockMergeSort(float *a, int n)
{
	/*
	1.RegisterSort
	2.iterate merging sort
	*/
	int bSize = SIMDWIDTH*SIMDWIDTH;
	//registerSort
	int bNum = n / bSize;
	float* tempF = (float*)_mm_malloc(sizeof(float)*SIMDWIDTH, WORDLENGTH);
	for (int i = 0; i < bNum; i++){
		__m128 m[SIMDWIDTH];
		for (int j = 0; j < SIMDWIDTH; j++){
			memcpy(tempF, &a[i*bSize+j*SIMDWIDTH], sizeof(float)*SIMDWIDTH);
			m[j] = _mm_load_ps(tempF);
		}
		registerSort(m, SIMDWIDTH);
		for (int j = 0; j < SIMDWIDTH; j++){
			_mm_store_ps(tempF, m[j]);
			memcpy(&a[i*bSize + j*SIMDWIDTH], tempF, sizeof(float)*SIMDWIDTH);
		}	
	}
	_mm_free(tempF);

	//merge
	int itr = log(n / SIMDWIDTH) / log(2);
	for (int i = 0; i < itr; i++){
		int bbSize = pow(2, i)*SIMDWIDTH;
		int bbNum = n / bbSize;
		float *tempR = (float*)_mm_malloc(sizeof(float)*bbSize * 2, WORDLENGTH);
		//printf("bbSize=%d,bbNum=%d\n", bbSize, bbNum);
		for (int j = 0; j < bbNum / 2; j++){
			blockInnerParallelMerge(&a[2 * j*bbSize], &a[(2 * j + 1)*bbSize], tempR, bbSize, THREAD_NUM);
			memcpy(&a[2 * j*bbSize], tempR, sizeof(float)*bbSize * 2);
		}
		_mm_free(tempR);
	}


}




