#include <cmath>
#include <fstream>
#include <vector>

#include "BallTree.h"

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
    DATA_SLOTSIZE = sizeof(int) * 3 + sizeof(float) * dimension * N0 + sizeof(bool);
    INDEX_SLOTSIZE = sizeof(int) * 3 + sizeof(float) * (dimension + 1) + sizeof(bool);
    numIndexSlot = floor(PAGE_SIZE / INDEX_SLOTSIZE);
    numDataSlot = floor(PAGE_SIZE / DATA_SLOTSIZE);

    id = 0;
    int fatherId = -1;
    buildBall(root, n, d, data, fatherId, true);

    for (int i = 0; i < n; i++)
        delete[] data[i];
    delete[] data;
    if (root == NULL)
        return false;
    return true;
}

void BallTree::buildBall(ballTreeNode *&node, int n, int d, float **data, int fatherId, bool isLeft) {
    float *mean = computeMean(n, d, data);
    node = new ballTreeNode(computeRadius(n, d, data, mean), mean, d);
    node->id = id++;
    node->fatherId = fatherId;
    node->isLeft = isLeft;

    delete[] mean;

    if (n <= N0) {
        node->table = new float*[n];
        node->tableSize = n;
        for (int i = 0; i < n; i++) {
            node->table[i] = new float[d];
            memcpy(node->table[i], data[i], d * sizeof(float));
        }
        return;
    }

    float *A, *B;
    MakeBallTreeSplit(n, d, data, A, B);

    std::vector<float*> leftData, rightData;
    for (int i = 0; i < n; i++) {
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
    float *mean = new float[d];

    for (int i = 0; i < d; i++) {
        float tempSum = 0;
        for (int j = 0; j < n; j++) {
            tempSum += data[j][i];
        }
        mean[i] = tempSum / n;
    }
    return mean;
}

float BallTree::computeRadius(int n, int d, float **data, float *mean) {
    float max = 0;
    for (int i = 0; i < n; i++) {
        float radius = 0;
        for (int j = 0; j < d; j++) {
            radius += pow(fabs(mean[j] - data[i][j]), 2);
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
    for (int i = 0; i < dimension; i++)
        squareSum += (x[i] - y[i]) * (x[i] - y[i]);
    return sqrt(squareSum);
}

bool BallTree::MakeBallTreeSplit(int n, int d, float **data, float *&A, float *&B) {
    float *pick = data[0];
    float maxDistance = 0;
    for (int i = 0; i < n; i++) {
        float curDistance = computeDistance(pick, data[i]);
        if (curDistance > maxDistance) {
            maxDistance = curDistance;
            A = data[i];
        }
    }
    maxDistance = 0;
    for (int i = 0; i < n; i++) {
        float curDistance = computeDistance(A, data[i]);
        if (curDistance > maxDistance) {
            maxDistance = curDistance;
            B = data[i];
        }
    }
    return true;
}