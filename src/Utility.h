#ifndef __UTILITY_H
#define __UTILITY_H

#define L 256

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

    ballTreeNode();
    ballTreeNode(float r, float* m, int d);
};

bool read_data(int n, int d, float** &data, const char* file_name);

float **parseFloatArr(std::vector<float*> v);

void storeNode(ballTreeNode *node, std::ofstream &output, int numSlot, int d);

#endif
