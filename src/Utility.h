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
    int id;
    int fatherId;
    bool isLeft;
    float **table;

    ballTreeNode *left;
    ballTreeNode *right;

    ballTreeNode();
    ballTreeNode(float r, float* m, int d);
};

const int INDEX_SLOTSIZE = sizeof(float) * 2 + sizeof(int) * 3 + sizeof(bool);

bool read_data(int n, int d, float** &data, const char* file_name);

float **parseFloatArr(std::vector<float*> v);

void storeIndexNode(ballTreeNode *node, std::ofstream &output, int numSlot, int d, std::streampos &filePtr);
void storeDataNode(ballTreeNode *node, std::ofstream &output, int numSlot, int d, std::streampos &filePtr);

int parsePageId(std::streampos filePtr);
int parseSlotId(std::streampos filePtr, int slotsize, int numSlot);

#endif
