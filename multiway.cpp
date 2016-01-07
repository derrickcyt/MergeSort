#include<iostream>
#include"merge_sort.h"
#include<math.h>

#define CACHESIZE 1024*1024*4
#define SIMDWIDTH 4
#define WORDLENGTH 16
#define PARALLEL_SIZE_LIMIT 1024*1024
#define THREAD_NUM 4
#define FIFO_SIZE 4*1024


struct Root{
	float *data;
	int len;
};

struct Leaf{
	float *sortedData;
	int pos;
	int len;
};

struct Nodes{
	struct Node *nodes;
	int realNum;
	int nodeDataLen;
};

struct Node{
	float *data;
	int head;
	int tail;
	int isEnd;
	int nodeDataLen;
};


void initNode(struct Node *node, int nodeDataLen)
{
	node = (struct Node*)malloc(sizeof(struct Node));
	node->head = node->tail = node->isEnd = 0;
	node->data = (float*)malloc(sizeof(float)*nodeDataLen);
	node->nodeDataLen = nodeDataLen;
}

void initNodes(struct Nodes *nodes, int realNum, int nodeDataLen)
{
	nodes = (struct Nodes*)malloc(sizeof(struct Nodes));
	nodes->realNum = realNum;
	nodes->nodeDataLen = nodeDataLen;
	nodes->nodes = (struct Node*)malloc(sizeof(struct Node)*realNum);
	for (int i = 0; i < realNum; i++) initNode(&nodes->nodes[i], nodeDataLen);
}

int getNodeLen(struct Node *node,int nodeDataLen)
{	
	return (node->tail - node->head) % nodeDataLen;
}

void enQueue(float *data, int len, struct Node *node)
{
	int tailRest = node->nodeDataLen - node->tail;
	if (len <= tailRest){
		memcpy(&node->data[node->tail], data, sizeof(float)*len);
		node->tail = (node->tail+len)%node->nodeDataLen;
	}
	else{
		memcpy(&node->data[node->tail], data, sizeof(float)*tailRest);
		len -= tailRest;
		memcpy(node->data, &data[tailRest], sizeof(float)*len);
		node->tail = len;
	}
}

float outQueue(struct Node *node)
{
	int f = node->data[node->head];
	node->head = (node->head + 1) % node->nodeDataLen;
	return f;
}

void outQueue(float *data, int len, struct Node *node)
{
	for (int i = 0; i<len; i++)
		data[i] = outQueue(node);
}

int checkIfReady(struct Node *node,int nodeDataLen)
{
	return getNodeLen(node, nodeDataLen) > (int)nodeDataLen/2 ? 1 : 0;
}

//merge sons' data to node
void mergeNode(int nodeId, struct Nodes *nodes)
{
	//merge opt
	Node *son1 = &nodes->nodes[nodeId * 2];
	Node *son2 = &nodes->nodes[nodeId * 2+1];
	int nowlen = 0, limit = nodes->nodeDataLen / 2;
	float *tempF = (float*)_mm_malloc(sizeof(float)*SIMDWIDTH, WORDLENGTH);
	float *v1 = (float*)_mm_malloc(sizeof(float)*SIMDWIDTH, WORDLENGTH);
	float *v2 = (float*)_mm_malloc(sizeof(float)*SIMDWIDTH, WORDLENGTH);
	while (nowlen < limit){
		outQueue(v1, SIMDWIDTH, son1);
		outQueue(v2, SIMDWIDTH, son2);
		__m128 m1 = _mm_load_ps(v1);
		__m128 m2 = _mm_load_ps(v2);
		bitonicSort(&m1, &m2);
		_mm_store_ps(tempF, m1);
		enQueue()




	}

}

void fetch2Ready(int nodeId,struct Nodes *nodes)
{
	int sonId = nodeId*2;
	//check if this son is not leaf
	if (sonId < nodes->realNum + 1){
		if (!checkIfReady(&nodes->nodes[sonId], nodes->nodeDataLen))
			fetch2Ready(sonId, nodes);
		if (!checkIfReady(&nodes->nodes[sonId+1], nodes->nodeDataLen))
			fetch2Ready(sonId+1, nodes);
		mergeNode(nodeId, nodes);
	}else{
		//fetch data from L

	}




		
	
}

void mergeSort(float *a, int n)
{
	//divide into block.size of each block is C/2e.
	int bSize = CACHESIZE*0.5 / sizeof(float);
	int bNum = ceil((double)n / bSize);

	for (int i = 0; i < bNum; i++)
	{
		int start = i*bSize;
		blockMergeSort(&a[start], bSize);
	}






}