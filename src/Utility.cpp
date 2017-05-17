#include <cstdio>
#include <cmath>
#include <algorithm>

#include "Utility.h"

using namespace std;

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

void storeNode(ballTreeNode *node, std::ofstream &output, int pageid, int slotid, int numSlot, int d) {
    if (slotid == 0) {    // bitMap needs to be inserted.
        output.seekp(pageid * 65536);
        bool *bitMap = new bool[numSlot];
        bitMap[0] = 1;
        memset(bitMap + 1, 0, numSlot - 1);
        output.write((char*)bitMap, numSlot);
        delete[] bitMap;
    }
    output.seekp(pageid * 65536 + slotid + numSlot);
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
    memcpy(floatArr + 1, node->mean, d);
    output.write((char*)intArr, 7 * sizeof(int));
    output.write((char*)floatArr, (d + 1) * sizeof(float));
    delete[] floatArr;
    delete[] intArr;
}