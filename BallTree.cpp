#include <cmath>
#include <stdlib.h>
#include <time.h>

#include "BallTree.h"
#include "Utility.h"

bool BallTree::buildTree(int n, int d, float** data) {
	int *dataIndex = new int[n];
	for (int i = 0; i < n; i++) {
		dataIndex[i] = i + 1;
	}
	//代表集合的索引

	root=NULL;
	//生成根节点。
	root = MakeBallTree(n, d, data, dataIndex);
	//在MakeBallTree里面delete dataIndex

}

bool BallTree::storeTree(const char index_path) {

}

bool BallTree::restoreTree(const char* index_path) {

}

int BallTree::mipSearch(int d, float* query) {

}

bool BallTree::insertData(int d, float* data) {

}

bool BallTree::deleteData(int d, float* data) {

}

bool BallTree::buidlQuadTree(int n, int d, float** data) {

}

treeNode* MakeBallTree(int n, int d, float** data, int *dataIndex) {
	treeNode *T;
	float* mean = getMean(data, dataIndex, n, d);
	float R = getRadius(mean, data, dataIndex, n, d);
	T=new treeNode(mean, R, n);

	if (n <= N0) {
		//把数据放到子叶节点里
		//给treeNode加数据
		return T;
	}
	else {
		AB ab = MakeBallTreeSplit(data, dataIndex, n, d);
		//得到a点 b点
		delete [] dataIndex;
		//ab.AIndex ab.BIndex空间在下一层释放
		T->left = MakeBallTree(ab.asize, d, data, ab.AIndex);
		T->right = MakeBallTree(ab.bsize, d, data, ab.BIndex);
		return T;
	}
}

AB MakeBallTreeSplit(float ** data, int* dataIndex, int size, int d) {
	AB ab;
	ab.AIndex=new int[size];
	ab.BIndex=new int[size];
	//在外面delete掉
	ab.asize=0;
	ab.bsize=0;

	srand((unsigned)time(NULL));
	int randomPointIndex = (rand() % (size-0))+ 0;
	//随机找到一个点的坐标
	float *tempA;
	float *tempB;
	int realIndex = dataIndex[randomPointIndex];
	tempA = getMaxPoint(data[realIndex], data, dataIndex, size, d);
	//找到A点
	tempB = getMaxPoint(a, data, dataIndex, size, d);
	//找到B点

	//分出A集合，B集合
	for (int i=0; i<size; i++){
		float r1=0; 
		float r2=0;
		for (int j=0; j<d; j++){
			r1+=pow(tempA[j]-data[i][j], 2);
			r2+=pow(tempB[j]-data[i][j], 2);
		}
		r1=sqrt(r1);
		r2=sqrt(r2);
		if (r1<=r2) {
			ab.AIndex[ab.asize]=i;
			ab.asize++;
		}
		else {
			ab.BIndex[ab.bsize]=i;
			ab.bsize++;
		}
	}


	// ab.d = d;
	return ab;
}


float* getMean(float **data, int *dataIndex, int size, int dim) {
	float mean=new float[dim];

}

float getRadius(float* mean, float **data, int *dataIndex, int size, int dim) {
	return getMax(mean, data, dataIndex, size, dim);
}

float getMax(float* mean, float **data, int *dataIndex, int size, int dim) {
	float max = 0;
	for (int i = 0; i < size; i++) {
		float radius = 0;
		int index = dataIndex[i];
		//得到相对于data里面的index
		for (int j = 0; j < d; j++) {
			radius += pow(fabs(mean[j] - data[index][j]));
		}

		radius = sqrt(radius);
		//得到根号值，即半径
		if (max < radius) {
			max = radius;
		}
	}
	return max;
}

float* getMaxPoint(float* mean, float **data, int *dataIndex, int size, int d) {
	float max = 0;
	int maxIndex = 0;
	for (int i = 0; i < size; i++) {
		float radius = 0;
		int index = dataIndex[i];
		//得到相对于data里面的index
		for (int j = 0; j < d; j++) {
			radius += pow(fabs(mean[j] - data[index][j]), 2);
		}

		radius = sqrt(radius);
		//得到根号值，即半径
		if (max < radius) {
			maxIndex = index;
			max = radius;
		}
	}
	return data[maxIndex];
}
