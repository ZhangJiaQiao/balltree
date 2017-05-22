#include <cmath>
#include <fstream>
#include <vector>

#include "BallTree.h"

BallTree::BallTree() {
    INDEX_SLOTSIZE = DATA_SLOTSIZE = 0;
    curIndexRid = curDataRid = Rid();
    root = NULL;
    dimension = 0;
    numIndexSlot = 0;
    numDataSlot = 0;
}

BallTree::~BallTree() {
    ;
}

bool BallTree::buildTree(int n, int d, float **data) {
    /* Initialize some values at buiding time. */
    dimension = d;
    INDEX_SLOTSIZE = sizeof(int) * INDEX_INT_SIZE + sizeof(float) * (dimension + 1) + sizeof(bool) * 2;
    DATA_SLOTSIZE = sizeof(int) * DATA_INT_SIZE + sizeof(float) * dimension * N0;
    numIndexSlot = floor(PAGE_SIZE / INDEX_SLOTSIZE);
    numDataSlot = floor(PAGE_SIZE / DATA_SLOTSIZE);

    buildBall(root, n, d, data);

    for (int i = 0; i < n; i++)
        delete[] data[i];
    delete[] data;
    if (root == NULL)
        return false;
    return true;
}

void BallTree::buildBall(ballTreeNode *&node, int n, int d, float **data) {
    float *mean = computeMean(n, d, data);
    node = new ballTreeNode(computeRadius(n, d, data, mean), mean, d);

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
        printf("In record %d of %d datas.\n", i, n);
        if (computeDistance(data[i], A) <= computeDistance(data[i], B)) {
            leftData.push_back(data[i]);
        } else {
            rightData.push_back(data[i]);
        }
    }

    float **leftData_ = parseFloatArr(leftData);
    float **rightData_ = parseFloatArr(rightData);
    node->isLeftLeaf = leftData.size() <= N0;
    node->isRightLeaf = rightData.size() <= N0;
    buildBall(node->left, leftData.size(), d, leftData_);
    buildBall(node->right, rightData.size(), d, rightData_);
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

    curIndexRid = curDataRid = Rid(0, 0);
    preorderStore(root, NULL, indexOutput, dataOutput, false);
    indexOutput.close();
    dataOutput.close();
    return true;
}

void BallTree::preorderStore(ballTreeNode *node, ballTreeNode *father, std::ofstream &indexOutput, 
    std::ofstream &dataOutput, bool isLeft) {
    if (node == NULL)
        return;
    if (node->table == NULL) {
        if (father != NULL) {
            if (isLeft) father->leftRid = curIndexRid;
            else father->rightRid = curIndexRid;
        }
        storeIndexNode(node, indexOutput, curIndexRid);
        curIndexRid.pageid = curIndexRid.slotid == numIndexSlot - 1 ? curIndexRid.pageid + 1 : curIndexRid.pageid;
        curIndexRid.slotid = curIndexRid.slotid == numIndexSlot - 1 ? 0 : curIndexRid.slotid + 1;
    } else {
        if (father != NULL) {
            if (isLeft) father->leftRid = curDataRid;
            else father->rightRid = curDataRid;
        }
        curDataRid.pageid = curDataRid.slotid == numDataSlot - 1 ? curDataRid.pageid + 1 : curDataRid.pageid;
        curDataRid.slotid = curDataRid.slotid == numDataSlot - 1 ? 0 : curDataRid.slotid + 1;
        storeDataNode(node, dataOutput, curDataRid);
    }

    preorderStore(node->left, node, indexOutput, dataOutput, true);
    preorderStore(node->right, node, indexOutput, dataOutput, false);
}

void BallTree::storeIndexNode(ballTreeNode *node, std::ofstream &output, Rid &rid) {
    int pageid = rid.pageid;
    int slotid = rid.slotid;

    if (slotid == 0) {    // bitMap needs to be inserted.
        output.seekp(PAGE_SIZE * pageid);
        bool *bitMap = new bool[numIndexSlot];
        memset(bitMap, 0, numIndexSlot);
        output.write((char*)bitMap, numIndexSlot);
        delete[] bitMap;
    }
    output.seekp(pageid * PAGE_SIZE + slotid);
    bool a = true;
    output.write((char*)&a, sizeof(bool));

    output.seekp(pageid * PAGE_SIZE + slotid * INDEX_SLOTSIZE + numIndexSlot);
    float *floatArr = new float[dimension + 1];
    int *intArr = new int[INDEX_INT_SIZE];
    intArr[0] = node->leftRid.pageid;
    intArr[1] = node->leftRid.slotid;
    intArr[2] = node->rightRid.pageid;
    intArr[3] = node->rightRid.slotid;
    floatArr[0] = node->radius;
    memcpy(floatArr + 1, node->mean, dimension * sizeof(float));
    bool *boolArr = new bool[2];
    boolArr[0] = node->isLeftLeaf;
    boolArr[1] = node->isRightLeaf;
    output.write((char*)intArr, INDEX_INT_SIZE * sizeof(int));
    output.write((char*)floatArr, (dimension + 1) * sizeof(float));
    output.write((char*)boolArr, INDEX_BOOL_SIZE * sizeof(bool));

    delete[] floatArr;
    delete[] intArr;
}

void BallTree::storeDataNode(ballTreeNode *node, std::ofstream &output, Rid &rid) {
    int pageid = rid.pageid;
    int slotid = rid.slotid;

    if (slotid == 0) {    // bitMap needs to be inserted.
        output.seekp(PAGE_SIZE * pageid);
        bool *bitMap = new bool[numDataSlot];
        memset(bitMap, 0, numDataSlot);
        output.write((char*)bitMap, numDataSlot);
        delete[] bitMap;
    }
    output.seekp(pageid * PAGE_SIZE + slotid);
    bool a = true;
    output.write((char*)&a, sizeof(bool));

    output.seekp(pageid * PAGE_SIZE + slotid * DATA_SLOTSIZE + numDataSlot);
    int numTuples = node->tableSize;
    float *floatArr = new float[dimension * node->tableSize];
    for (int i = 0; i < node->tableSize; i++)
        memcpy(floatArr + i * dimension, node->table[i], dimension * sizeof(float));
    output.write((char*)&numTuples, sizeof(int));
    output.write((char*)floatArr, dimension * node->tableSize * sizeof(float));

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