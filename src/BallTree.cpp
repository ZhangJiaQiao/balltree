#include <cmath>
#include <fstream>
#include <vector>

#include "BallTree.h"

BallTree::BallTree() {
    root = NULL;
    dimension = 0;
    numSlot = 0;
    pageid = 0;
    slotid = 0;
}

BallTree::~BallTree() {
    ;
}

bool BallTree::buildTree(int n, int d, float **data) {
    dimension = d;
    numSlot = floor(65536 / (sizeof(int) * 7 + sizeof(float) * (d + 1)));
    buildBall(root, n, d, data);

    if (root == NULL)
        return false;
    return true;
}

bool BallTree::storeTree(const char *index_path) {
    std::ofstream output(index_path, std::ofstream::out | std::ofstream::binary);

    if (!output)
        return false;

    pageid = slotid = 0;
    preorderSetRid(root, NULL, true);
    preorderStore(root, output);
    output.close();
    return true;
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

void BallTree::preorderSetRid(ballTreeNode *node, ballTreeNode *father, bool isLeft) {
    if (node == NULL)
        return;
    node->rid.pageid = slotid == numSlot - 1 ? pageid++ : pageid;
    node->rid.slotid = slotid++;
    if (slotid == numSlot - 1) {
        node->rid.slotid = 0;
    }
    if (father != NULL) {
        if (isLeft)
            father->left_rid = node->rid;
        else
            father->right_rid = node->rid;
    }
    preorderSetRid(node->left, node, true);
    preorderSetRid(node->right, node, false);
}

void BallTree::preorderStore(ballTreeNode *node, std::ofstream &output) {
    if (node == NULL)
        return;
    storeNode(node, output, numSlot, dimension);
    preorderStore(node->left, output);
    preorderStore(node->right, output);
}

void BallTree::buildBall(ballTreeNode *&node, int n, int d, float **data) {
    float *mean = computeMean(n, d, data);
    node = new ballTreeNode(computeRadius(n, d, data, mean), mean, d);

    if (n <= N0)
        return;

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
    buildBall(node->left, leftData.size(), d, leftData_);
    buildBall(node->right, rightData.size(), d, rightData_);
    delete[] leftData_;
    delete[] rightData_;
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