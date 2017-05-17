#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>

#include "Utility.h"

using namespace std;

ballTreeNode::ballTreeNode() {
    radius = 0;
    mean = NULL;
    rid = Rid();
    left = right = NULL;
}
ballTreeNode::ballTreeNode(float r, float* m, int d) {
    rid = Rid();
    radius = r;
    mean = m;
    left = right = NULL;
}

bool read_data(int n, int d, float** &data, const char* file_name) {
	FILE* fin = fopen(file_name, "r");
	if (!fin) {
		printf("%s doesn't exist!\n", file_name);
		return false;
	}

	int id;
	data = new float*[n];
	for (int i = 0; i < n; i++) {
		data[i] = new float[d];
		fscanf(fin, "%d", &id);
		for (int j = 0; j < d; j++) {
			fscanf(fin, "%f", &data[i][j]);
		}
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

void storeNode(ballTreeNode *node, std::ofstream &output, int numSlot, int d) {
    if (node->rid.slotid == 0) {    // bitMap needs to be inserted.
        output.seekp(node->rid.pageid * 65536);
        bool *bitMap = new bool[numSlot];
        memset(bitMap, 0, numSlot);
        output.write((char*)bitMap, numSlot);
        delete[] bitMap;
    }
    output.seekp(node->rid.pageid * 65536 + node->rid.slotid);
    bool a = true;
    output.write((char*)&a, sizeof(bool));
    
    output.seekp(node->rid.pageid * 65536 + node->rid.slotid * (sizeof(int) * 7 + sizeof(float) * (d + 1)) + numSlot);
    float *floatArr = new float[d + 1];
    int *intArr = new int[7];
    intArr[0] = node->rid.pageid;
    intArr[1] = node->rid.slotid;
    intArr[2] = node->left_rid.pageid;
    intArr[3] = node->left_rid.slotid;
    intArr[4] = node->right_rid.pageid;
    intArr[5] = node->right_rid.slotid;
    intArr[6] = d;
    floatArr[0] = node->radius;
    for (int i = 0; i < d; i++)
        floatArr[i + 1] = node->mean[i];
    output.write((char*)intArr, 7 * sizeof(int));
    output.write((char*)floatArr, (d + 1) * sizeof(float));
    delete[] floatArr;
    delete[] intArr;
}