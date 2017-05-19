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
    tableSize = 0;
}
ballTreeNode::ballTreeNode(float r, float* m, int d) {
    radius = r;
    mean = new float[d];
    memcpy(mean, m, d * sizeof(float));
    left = right = NULL;
    id = fatherId = -1;
    isLeft = false;
    table = NULL;
    tableSize = 0;
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
        printf("Data %d loaded.\n", id);
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
    return filePtr/PAGE_SIZE;
}

int parseSlotId(std::streampos filePtr, int slotsize, int numSlot) {
    if (filePtr % PAGE_SIZE == 0)
        return 0;
    return (int)(((filePtr%PAGE_SIZE)-numSlot)/slotsize);
}