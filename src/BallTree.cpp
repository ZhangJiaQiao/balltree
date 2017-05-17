#include "BallTree.h"
#include "Utility.h"

BallTree::BallTree() {
    root = nullptr;
    dimension = 0;
    numSlot = 0;
}

BallTree::~BallTree() {
    ;
}

bool BallTree::buildTree(int n, int d, float **data) {
    buildBall(root, n, d, data);
    dimension = d;
    numSlot = floor(65536 / (3 * sizeof(Rid) + (d + 1) * sizeof(float) + 1));

    if (root == nullptr)
        return false;
    return true;
}

bool BallTree::storeTree(const char *index_path) {
    std::ofstream output(index_path, std::ofstream::out | std::ofstream::binary);

    if (!output)
        return false;

    int curPage = 0, curSlot = 0;
    preorderStore(root, output, curPage, curSlot);
    return true;
}

float BallTree::computeDistance(float *x, float *y) {
    float squareSum = 0;
    for (int i = 0; i < dimension; i++)
        squareSum += x[i] * x[i] - y[i] * y[i];
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
}

void BallTree::preorderStore(ballTreeNode *node, std::ofstream &output, int pageid, int slotid) {
    if (node == nullptr)
        return;
    storeNode(node, output, pageid, slotid, numSlot, dimension);
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
    delete[] leftData_;
    delete[] rightData_;
}
