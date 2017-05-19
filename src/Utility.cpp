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
    left = right = NULL;
    id = fatherId = -1;
    isLeft = false;
    table = NULL;
}
ballTreeNode::ballTreeNode(float r, float* m, int d) {
    radius = r;
    mean = new float[d];
    memcpy(mean, m, d);
    left = right = NULL;
    id = fatherId = -1;
    isLeft = false;
    table = NULL;
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

void storeIndexNode(ballTreeNode *node, std::ofstream &output, int numSlot, int d, std::streampos &filePtr) {
    int pageid = parsePageId(filePtr);
    int slotid = parseSlotId(filePtr, INDEX_SLOTSIZE, numSlot);
    if (slotid == 0) {    // bitMap needs to be inserted.
        output.seekp(filePtr);
        bool *bitMap = new bool[numSlot];
        memset(bitMap, 0, numSlot);
        output.write((char*)bitMap, numSlot);
        filePtr = output.tellp();
        delete[] bitMap;
    }
    output.seekp(pageid * 65536 + slotid);
    bool a = true;
    output.write((char*)&a, sizeof(bool));
    
    output.seekp(pageid * 65536 + slotid * INDEX_SLOTSIZE + numSlot);
    float *floatArr = new float[d + 1];
    int *intArr = new int[3];
    intArr[0] = node->id;
    intArr[1] = node->fatherId;
    intArr[2] = d;
    floatArr[0] = node->radius;
    memcpy(floatArr + 1, node->mean, d);
    output.write((char*)intArr, 3 * sizeof(int));
    output.write((char*)floatArr, (d + 1) * sizeof(float));
    delete[] floatArr;
    delete[] intArr;
}

void storeDataNode(ballTreeNode *node, std::ofstream &output, int numSlot, int d, std::streampos &filePtr) {
    ;
}

int parsePageId(std::streampos filePtr) {
    return filePtr/65536;
}

int parseSlotId(std::streampos filePtr, int slotsize, int numSlot) {
    return (int)(((filePtr%65536)-numSlot)/slotsize);
}