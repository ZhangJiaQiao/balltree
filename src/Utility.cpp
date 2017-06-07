	#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include<cstring>

#include "Utility.h"

using namespace std;

ballTreeNode::ballTreeNode() {
	radius = 0;
	mean = NULL;
	leftRid = rightRid = upRid = downRid = myRid = Rid();
	isLeftLeaf = isRightLeaf = isUpLeaf = isDownLeaf = false;
	isLeft = isRight = isUp = isDown = false;
	table = NULL;
	tableSize = 0;
	left = right = up = down = father = NULL;
}
ballTreeNode::ballTreeNode(float r, float* m, int d) {
	radius = r;
	mean = new float[d - 1];
	memcpy(mean, m, (d - 1) * sizeof(float));
	leftRid = rightRid = upRid = downRid = myRid = Rid();
	isLeftLeaf = isRightLeaf = isUpLeaf = isDownLeaf = false;
	isLeft = isRight = isUp = isDown = false;
	table = NULL;
	tableSize = 0;
	left = right = up = down = father = NULL;
}

void ballTreeNode::delete_memory() {
	if (table != NULL) {
		for (int i = 0; i < tableSize; i++) {
			delete[] table[i];
		}
		delete[] table;
		table = NULL;
	}
}
void ballTreeNode::release_right() {
	if (isRightLeaf&&right != NULL) {
		right->delete_memory();
		delete[]right->mean;
		delete right;
		right = NULL;
	}
}
void ballTreeNode::release_left() {
	if (isLeftLeaf&&left != NULL) {
		left->delete_memory();
		delete[]left->mean;
		delete left;
		left = NULL;
	}
}
void ballTreeNode::release_up() {
	if (isUpLeaf&&up != NULL) {
		up->delete_memory();
		delete[]up->mean;
		delete up;
		up = NULL;
	}
}
void ballTreeNode::release_down() {
	if (isDownLeaf&&down != NULL) {
		down->delete_memory();
		delete[]down->mean;
		delete down;
		down = NULL;
	}
}
bool read_data(int n, int d, float** &data, const char* file_name) {
	FILE* fin = fopen(file_name, "r");
	if (!fin) {
		printf("%s doesn't exist!\n", file_name);
		return false;
	}

	data = new float*[n];
	for (int i = 0; i < n; i++) {
		data[i] = new float[d];
		for (int j = 0; j < d; j++) {
			fscanf(fin, "%f", &data[i][j]);
		}
        printf("Data %d loaded.\n", (int)data[i][0]);
	}

	printf("Finish reading %s\n", file_name);
	fclose(fin);

	return true;
}

float **parseFloatArr(std::vector<float*> v) {
    float **data = new float*[v.size()];
    for (int i = 0; i < v.size(); i++) {
        data[i] = v[i];
    }
    return data;
}

int parsePageId(std::streampos filePtr) {
    return (int)(filePtr / PAGE_SIZE);
}

int parseSlotId(std::streampos filePtr, int slotsize, int numSlot) {
    if (filePtr % PAGE_SIZE == 0)
        return 0;
    return (int)(((filePtr%PAGE_SIZE)-numSlot)/slotsize);
}
