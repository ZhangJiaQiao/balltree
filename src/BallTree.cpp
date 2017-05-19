#include <cmath>
#include <fstream>
#include <vector>

#include "BallTree.h"

BallTree::BallTree() {
    root = NULL;
    dimension = 0;
    numSlot = 0;
    id = 0;
}

BallTree::~BallTree() {
    ;
}

bool BallTree::buildTree(int n, int d, float **data) {
    dimension = d;
    DATA_SLOTSIZE = sizeof(float) * dimension * N0 + sizeof(int);
    numSlot = floor(65536 / INDEX_SLOTSIZE);
    id = 0;
    int fatherId = -1;
    buildBall(root, n, d, data, fatherId);

    if (root == NULL)
        return false;
    return true;
}

void BallTree::buildBall(ballTreeNode *&node, int n, int d, float **data, int fatherId) {
    float *mean = computeMean(n, d, data);
    node = new ballTreeNode(computeRadius(n, d, data, mean), mean, d);
    node->id = id++;
    node->fatherId = fatherId;

    if (n <= N0) {
        node->table = new float*[n];
        for (int i = 0; i < n; i++) {
            node->table[i] = new float[d];
            memcpy(node->table[i], data[i], d);
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
    buildBall(node->left, leftData.size(), d, leftData_, node->id);
    buildBall(node->right, rightData.size(), d, rightData_, node->id);
    delete[] leftData_;
    delete[] rightData_;
}

bool BallTree::storeTree(const char *index_path) {
    char indexEntryPath[L], dataEntryPath[L];
    sprintf(indexEntryPath, "%s/indexEntries.dat", index_path);
    sprintf(dataEntryPath, "%s/dataEntries.dat", index_path);
    std::ofstream indexOutput(indexEntryPath, std::ofstream::out | std::ofstream::binary);
    std::ofstream dataOutput(indexEntryPath, std::ofstream::out | std::ofstream::binary);

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
        storeIndexNode(node, indexOutput, numSlot, dimension, indexPtr);
    else
        storeDataNode(node, dataOutput, numSlot, dimension, dataPtr);
    preorderStore(node->left, indexOutput, dataOutput, indexPtr, dataPtr);
    preorderStore(node->right, indexOutput, dataOutput, indexPtr, dataPtr);
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