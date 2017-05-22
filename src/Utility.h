#ifndef __UTILITY_H
#define __UTILITY_H

#define L 256
#define PAGE_SIZE 65536

const int INDEX_INT_SIZE = 4;
const int INDEX_BOOL_SIZE = 2;
const int DATA_INT_SIZE = 1;

const int METADATA_INDEX_OFFSET = sizeof(int) * 2;
const int METADATA_DATA_OFFSET = sizeof(int);

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
    float radius;  //球的半径
    float *mean;  //球的圆心
    
    Rid leftRid;  //  The rid of the left child.
    Rid rightRid;  //  The rid of the right child.

    bool isLeftLeaf;  //左节点是否为叶子
    bool isRightLeaf;  //右节点是否为叶子
    
    float **table;  //如果为叶子，table存数据
    int tableSize;  //table的大小

    ballTreeNode *left;  //左节点的索引
    ballTreeNode *right;  //右节点的索引

    ballTreeNode();
    ballTreeNode(float r, float* m, int d);
};

bool read_data(int n, int d, float** &data, const char* file_name);

float **parseFloatArr(std::vector<float*> v);

int parsePageId(std::streampos filePtr);
int parseSlotId(std::streampos filePtr, int slotsize, int numSlot);

//--------------------ZJQ:20170521任务3与任务4实现--------------------------//
struct Mip {
    float product;
    int index;
}
#endif
