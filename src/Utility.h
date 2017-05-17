#ifndef __UTILITY_H
#define __UTILITY_H

#define L 256

#include <vector>
#include <fstream>
#include <cmath>

struct Rid {
    int pageid;
    int slotid;

    Rid() {
        pageid = slotid = -1;
    }
    Rid(int p, int s) {
        pageid = p;
        slotid = s;
    }
};

struct ballTreeNode {
    float radius;
    float *mean;
    Rid rid;

    ballTreeNode *left;
    ballTreeNode *right;
    Rid left_rid;
    Rid right_rid;

    ballTreeNode() {
        radius = 0;
        mean = nullptr;
        rid = Rid();
        left = right = nullptr;
    }
    ballTreeNode(float r, float* m, int d) {
        rid = Rid();
        radius = r;
        mean = new float[d];
        memcpy(mean, m, d);
        left = right = nullptr;
    }
};

bool read_data(int n, int d, float** &data, const char* file_name);

float **parseFloatArr(std::vector<float*> v);

void storeNode(ballTreeNode *node, std::ofstream &output, int pageid, int slotid, int numSlot, int d);

#endif
