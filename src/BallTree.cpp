#include <cmath>
#include <fstream>
#include <vector>

#include "BallTree.h"

const int DATA_INT_SIZE = 5;
const int INDEX_INT_SIZE = 2;

BallTree::BallTree() {
    root = NULL;
    dimension = 0;
    numIndexSlot = 0;
    numDataSlot = 0;
    id = 0;
}

BallTree::~BallTree() {
    ;
}

bool BallTree::buildTree(int n, int d, float **data) {
    dimension = d;
    DATA_SLOTSIZE = sizeof(int) * DATA_INT_SIZE + sizeof(float) * dimension * N0 + sizeof(bool) * 2;
    INDEX_SLOTSIZE = sizeof(int) * INDEX_INT_SIZE + sizeof(float) * (dimension + 1) + sizeof(bool) * 2;
    numIndexSlot = floor(PAGE_SIZE / INDEX_SLOTSIZE);
    numDataSlot = floor(PAGE_SIZE / DATA_SLOTSIZE);

    curPageID = 0; //根节点页号槽号默认为(0,0)，剩下节点从(0,1)开始
    curSlotID = 1;
    //int fatherId = -1;
    buildBall(root, n, d, data, false, false);

    for (int i = 0; i < n; i++)
        delete[] data[i];
    delete[] data;
    if (root == NULL)
        return false;
    return true;
}

void BallTree::buildBall(ballTreeNode *&node, int n, int d, float **data, bool isLeftLeaf, bool isRightLeaf) {
    if (curSlotID >= numIndexSlot) {
        curPageID++;
    }

    float *mean = computeMean(n, d, data);
    node = new ballTreeNode(computeRadius(n, d, data, mean), mean, d);
    //node->id = id++;
    //node->fatherId = fatherId;
    //node->isLeft = isLeft;

    node->leftPageID
    delete[] mean;

    if (n <= N0) {
        node->table = new float*[n];
        node->tableSize = n;
        for (int i = 0; i < n; i++) {
            node->table[i] = new float[d];
            memcpy(node->table[i], data[i], d * sizeof(float));
        }
        printf("Tree node Data Entry %d built.\n", node->id);
        return;
    }

    printf("Tree node index Entry %d built.\n", node->id);

    float *A, *B;
    MakeBallTreeSplit(n, d, data, A, B);

    std::vector<float*> leftData, rightData;
    for (int i = 0; i < n; i++) {
        printf("In record %d of %d datas.\n", i, n);
        if (computeDistance(data[i], A) <= computeDistance(data[i], B)) {
            leftData.push_back(data[i]);
        } else {
            rightData.push_back(data[i]);
        }
    }

    float **leftData_ = parseFloatArr(leftData);
    float **rightData_ = parseFloatArr(rightData);
    buildBall(node->left, leftData.size(), d, leftData_, node->id, true);
    buildBall(node->right, rightData.size(), d, rightData_, node->id, false);
    delete[] leftData_;
    delete[] rightData_;
}

bool BallTree::storeTree(const char *index_path) {
    char indexEntryPath[L], dataEntryPath[L];
    sprintf(indexEntryPath, "%s/indexEntries.dat", index_path);
    sprintf(dataEntryPath, "%s/dataEntries.dat", index_path);
    std::ofstream indexOutput(indexEntryPath, std::ofstream::out | std::ofstream::binary);
    std::ofstream dataOutput(dataEntryPath, std::ofstream::out | std::ofstream::binary);

    if (!indexOutput || !dataOutput)
        return false;

    std::streampos indexPtr = 0, dataPtr = 0;
    preorderStore(root, indexOutput, dataOutput, indexPtr, dataPtr);
    indexOutput.close();
    dataOutput.close();
    return true;
}

void BallTree::preorderStore(ballTreeNode *node, std::ofstream &indexOutput, std::ofstream &dataOutput, std::streampos &indexPtr, std::streampos &dataPtr) {
    if (node == NULL)
        return;
    if (node->table == NULL)
        storeIndexNode(node, indexOutput, indexPtr);
    else
        storeDataNode(node, dataOutput, dataPtr);

    printf("Tree node %d stored.\n", node->id);

    preorderStore(node->left, indexOutput, dataOutput, indexPtr, dataPtr);
    preorderStore(node->right, indexOutput, dataOutput, indexPtr, dataPtr);
}

void BallTree::storeIndexNode(ballTreeNode *node, std::ofstream &output, std::streampos &filePtr) {
    int pageid = parsePageId(filePtr);
    int slotid = parseSlotId(filePtr, INDEX_SLOTSIZE, numIndexSlot);

    if (slotid == 0) {    // bitMap needs to be inserted.
        output.seekp(filePtr);
        bool *bitMap = new bool[numIndexSlot];
        memset(bitMap, 0, numIndexSlot);
        output.write((char*)bitMap, numIndexSlot);
        filePtr = output.tellp();
        delete[] bitMap;
    }
    output.seekp(pageid * PAGE_SIZE + slotid);
    bool a = true;
    output.write((char*)&a, sizeof(bool));

    output.seekp(pageid * PAGE_SIZE + slotid * INDEX_SLOTSIZE + numIndexSlot);
    float *floatArr = new float[dimension + 1];
    int *intArr = new int[3];
    intArr[0] = node->id;
    intArr[1] = node->fatherId;
    intArr[2] = dimension;
    floatArr[0] = node->radius;
    memcpy(floatArr + 1, node->mean, dimension * sizeof(float));
    bool isLeft = node->isLeft;
    output.write((char*)intArr, 3 * sizeof(int));
    output.write((char*)floatArr, (dimension + 1) * sizeof(float));
    output.write((char*)&isLeft, sizeof(bool));
    filePtr = output.tellp();

    if (slotid == numIndexSlot - 1) {
        filePtr = (pageid + 1) * PAGE_SIZE;
    }

    delete[] floatArr;
    delete[] intArr;
}

void BallTree::storeDataNode(ballTreeNode *node, std::ofstream &output, std::streampos &filePtr) {
    int pageid = parsePageId(filePtr);
    int slotid = parseSlotId(filePtr, DATA_SLOTSIZE, numDataSlot);

    if (slotid == 0) {    // bitMap needs to be inserted.
        output.seekp(filePtr);
        bool *bitMap = new bool[numDataSlot];
        memset(bitMap, 0, numDataSlot);
        output.write((char*)bitMap, numDataSlot);
        filePtr = output.tellp();
        delete[] bitMap;
    }
    output.seekp(pageid * PAGE_SIZE + slotid);
    bool a = true;
    output.write((char*)&a, sizeof(bool));

    output.seekp(pageid * PAGE_SIZE + slotid * DATA_SLOTSIZE + numDataSlot);
    int *intArr = new int[3];
    float *floatArr = new float[dimension * node->tableSize];
    intArr[0] = node->id;
    intArr[1] = node->fatherId;
    intArr[2] = node->tableSize;
    for (int i = 0; i < node->tableSize; i++)
        memcpy(floatArr + i * dimension, node->table[i], dimension * sizeof(float));
    bool isLeft = node->isLeft;
    output.write((char*)intArr, 3 * sizeof(int));
    output.write((char*)floatArr, dimension * node->tableSize * sizeof(float));
    output.write((char*)&isLeft, sizeof(bool));
    filePtr = pageid * PAGE_SIZE + numDataSlot + DATA_SLOTSIZE * (slotid + 1);

    if (slotid == numDataSlot - 1) {
        filePtr = (pageid + 1) * PAGE_SIZE;
    }

    delete[] intArr;
    delete[] floatArr;
}

float* BallTree::computeMean(int n, int d, float **data) {
    float *mean = new float[d - 1];

    for (int i = 0; i < d; i++) {
        float tempSum = 0;
        for (int j = 0; j < n; j++) {
            tempSum += data[j][i + 1];
        }
        mean[i] = tempSum / n;
    }
    return mean;
}

float BallTree::computeRadius(int n, int d, float **data, float *mean) {
    float max = 0;
    for (int i = 0; i < n; i++) {
        float radius = 0;
        for (int j = 0; j < d - 1; j++) {
            radius += pow(fabs(mean[j] - data[i][j + 1]), 2);
        }

        radius = sqrt(radius);
        if (max < radius) {
            max = radius;
        }
    }
    return max;
}

float BallTree::computeDistance(float *x, float *y) {
    float squareSum = 0;
    for (int i = 1; i < dimension; i++)
        squareSum += (x[i] - y[i]) * (x[i] - y[i]);
    return sqrt(squareSum);
}

bool BallTree::MakeBallTreeSplit(int n, int d, float **data, float *&A, float *&B) {
    float *pick = data[0];
    A = data[0];
    float maxDistance = 0;
    for (int i = 0; i < n; i++) {
        float curDistance = computeDistance(pick, data[i]);
        if (curDistance > maxDistance) {
            maxDistance = curDistance;
            A = data[i];
        }
    }
    maxDistance = 0;
    B = pick;
    for (int i = 0; i < n; i++) {
        float curDistance = computeDistance(A, data[i]);
        if (curDistance > maxDistance) {
            maxDistance = curDistance;
            B = data[i];
        }
    }
    return true;
}

//-----------------------ZJQ:20170521任务3与4实现-----------------------//
bool BallTree::restoreTree(const char* index_path) {
    std::ifstream infile(index_path, ios::binary);
    if (!infile) {
        cout << "open " << index_path << " failed!" << endl;
        return false;
    }
    //--------------下面实现的代码为读取根节点的所有信息，树的维度，槽长以及槽数---------------//

    //-------------------------------------------------------------------------------//
    return true;
}
//restore()函数将树的根节点设置即可，一次找一个节点-------uncompleted
int BallTree::mipSearch(int d, float* query) {
    dimension = d;

    Mip mip;
    mip.product = -1;
    mip.index = -1;
    //记录最大内积及目标序号
    TreeSearch(query, root, mip);
    return mip.index；
}

void BallTree::TreeSearch(float* query, ballTreeNode* node, Mip &mip) {
    if (mip.product < MIP(query, node)) {
        if (node->table == NULL) {
            LinearSearch(query, node, mip)
        }
        else {
            node->left = getNode(node->leftPageID, node->leftSlotID);
            node->right = getNode(node->rightPageID, node->rightSlotID);
            float leftProduct = MIP(query, node->left);
            float rightProduct = MIP(query, node->right);
            if (leftProduct < rightProduct) {
                TreeSearch(query, node->right, mip);
                delete node->right;
                TreeSearch(query, node->left, mip);
                delete node->left;
            }
            else {
                TreeSearch(query, node->left, mip);
                delete node->left;
                TreeSearch(query, node->right, mip);
                delete node->right;
            }
        }
    }
}
//mipSearch()与TreeSearch()分别对应论文的算法5,4的伪代码，实现检索功能--------completed
void BallTree::LinearSearch(float* query, ballTreeNode* node, Mip& mip) {
    for (int i = 0; i < node->tableSize; i++) {
        float newProduct = computeInnerProduct(query, node->table[i]);
        if (newProduct > mip.product) {
            mip.index = node->table[i][0];
            //这里的index位每一个数据项的第一个数
            mip.product = newProduct;
        }
    }
}
//对存有数据的叶子节点进行线性查找-------completed
float BallTree::computeInnerProduct(float* query, float* data) {
    float product = 0;
    for (int i = 0; i < dimension - 1; i++) {
        product += query[i] * data[i + 1];
    }
    return product;
}
//计算内积--------completed
float BallTree::MIP(float *query, ballTreeNode* node) {
    float product;
    for (int i = 0; i < dimension - 1; i++) {
        product += query[i] * node->mean[i];
    }
    product += node->radius;
    return product;
}
//估算球的上界与query的内积--------completed

ballTreeNode* BallTree::getNode(int pageID, int slotID) {
    ballTreeNode* node = new ballTreeNode();
    //---------------下面填写代码----------------//



    //-------------------------------------------//
    return node;
}
//根据页号和槽号获得节点数据--------------------umcompeleted