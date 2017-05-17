#include "BallTree.h"
#include "Utility.h"

BallTree::BallTree() {
    root = nullptr;
    dimension = 0;
}

BallTree::~BallTree() {
    ;
}

bool BallTree::buildTree(int n, int d, float **data) {
    buildBall(root, n, d, data);
    dimension = d;

    if (root == nullptr)
        return false;
    return true;
}

bool BallTree::storeTree(const char *index_path) {
    std::ofstream output(index_path, std::ofstream::out | std::ofstream::binary);

    if (!output)
        return false;

    int curPage = 0, curSlot = 0;
    //int numSlot = sizeof(float) * (dimension + 1) + sizeof(Rid) * 3;
    preorderStore(root, output, curPage, curSlot);
    return true;
}

void BallTree::preorderStore(ballTreeNode *node, std::ofstream &output, int pageid, int slotid) {
    if (node == nullptr)
        return;
    storeNode(node, output, pageid, slotid);
    preorderStore(node->left, output, slotid == numSlot - 1 ? pageid + 1 : pageid, slotid == numSlot - 1 ? 0 : slotid + 1);
    preorderStore(node->right, output, slotid == numSlot - 1 ? pageid + 1 : pageid, slotid == numSlot - 1 ? 0 : slotid + 1);
}

void BallTree::buildBall(ballTreeNode *&node, int n, int d, float **data) {
    node = new ballTreeNode(computeRadius(n, d, data), computeMean(n, d, data), d);

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
}

float* BallTree::computeMean(int n, int d, float **data) {
    float *mean = new float[d];

    for (int i = 0; i < d; i++) {
        float tempSum = 0;
        for (int j = 0; j < n; j++) {
            tempSum += data[j][i];
        }
        tempSum /= n;
        mean[i] = tempSum;
    }
}

float BallTree::computeRadius(int n, int d, float **data, float *mean) {
    float max = 0;
    for (int i = 0; i < n; i++) {
        float radius = 0;
        for (int j = 0; j < d; j++) {
            radius += pow(fabs(mean[j] - data[index][j]), 2);
        }

        radius = sqrt(radius);
        //得到根号值，即半径
        if (max < radius) {
            max = radius;
        }
    }
    return max;
}