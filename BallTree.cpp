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

	treeNode root
	//生成根节点。
	root = MakeBallTree(n, d, data, dataIndex);

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
	T = new treeNode();
	float* mean = getMean(data, dataIndex, n);
	float R = getRadius(mean, data, dataIndex, n, d);
	T->treeNode(mean, R, n);

	if (n <= N0) {
		return T;
	}
	else {
		AB ab = MakeBallTreeSplit(data, dataIndex, n, d);
		int *leftSet;
		int *rightSet;
		int leftSetSize = createLeftSet(n, d, data, dataIndex, leftSet, ab);
		int rightSetSize = createRightSet(d, data, leftSet);
		delete [] dataIndex;
		T->left = MakeBallTree(leftSetSize, d, data, a, rightSet);
		T->right = MakeBallTree(rightSetSize, d, data, rightSet);
		return T;
	}
}

AB MakeBallTreeSplit(float ** data, int* dataIndex, int size, int d) {
	srand((unsigned)time(NULL));
	int randomPointIndex = (rand() % (size-0))+ 0;
	//随机找到一个点的坐标
	float *temp;
	float *a = new float [d];
	int realIndex = dataIndex[randomPointIndex];
	temp = getMaxPoint(data[realIndex], data, dataIndex, size, d);
	for (int i = 0; i < d; i++) {
		a[i] = temp[i];
	}
	//找到A点
	temp = getMaxPoint(a, data, dataIndex, size, d);
	float *b = new float [d];
	for (int i = 0; i < d; i++) {
		b[i] = temp[i];
	}
	//找到B点
	AB ab;
	ab.A = a;
	ab.B = b;
	ab.d = d;
	return ab;
}

int createLeftSet(int n, int d, float **data, int *dataIndex, int *leftSet, AB ab) {

}
//生成左子树的点集的索引
int createRightSet(int d, float **data, int *leftSet, int *rightSet) {

}

float* getMean(float **data, int *dataIndex, int size) {

}

float getRadius(float* mean, float **data, int *dataIndex, int size, int d) {
	return getMax(mean, data, dataIndex, size, d);
}

float getMax(float* mean, float **data, int *dataIndex, int size, int d) {
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
			radius += pow(fabs(mean[j] - data[index][j]));
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