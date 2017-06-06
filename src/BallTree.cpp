#include <cmath>
#include <fstream>
#include <vector>
#include <queue>
#include <cstdio>
#include <cstring>

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
    INDEX_SLOTSIZE = sizeof(int) * INDEX_INT_SIZE + sizeof(float) * dimension + sizeof(bool) * 2;
    DATA_SLOTSIZE = sizeof(int) * DATA_INT_SIZE + sizeof(float) * dimension * N0 + sizeof(float) * dimension;
    numIndexSlot = floor(PAGE_SIZE / (INDEX_SLOTSIZE + 1));
    numDataSlot = floor(PAGE_SIZE / (DATA_SLOTSIZE + 1));

    buildBall(root, NULL, n, d, data, false);

    if (root == NULL)
        return false;

    return true;
}
int leaf_num = 0;
void BallTree::buildBall(ballTreeNode *&node, ballTreeNode *father, int &n, int &d, float **&data, bool isLeft) {
    float *mean = new float[d - 1];
    computeMean(mean, n, d, data);
    node = new ballTreeNode(computeRadius(n, d, data, mean), mean, d);
    node->father = father;
    node->isLeft = isLeft;

    delete[] mean;

    if (n <= N0) {
        node->table = new float*[n];
        node->tableSize = n;
        for (int i = 0; i < n; i++) {
			node->table[i] = data[i];
           // node->table[i] = new float[d];
           // memcpy(node->table[i], data[i], d * sizeof(float));
        }
		if(father != NULL) delete[]data;
		leaf_num++;
		printf("Make %d leaf nodes.\n", leaf_num);
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
	if (father!=NULL) delete[]data;

	int leftn = leftData.size();
	int rightn = rightData.size();
	//释放vector里的内存
    float **leftData_ = parseFloatArr(leftData);
	leftData.swap(std::vector<float*>());
    float **rightData_ = parseFloatArr(rightData);
	rightData.swap(std::vector<float*>());

    node->isLeftLeaf = leftn <= N0;
    node->isRightLeaf = rightn<= N0;
    buildBall(node->left, node, leftn, d, leftData_, true);
	//delete[] leftData_;
    buildBall(node->right, node, rightn, d, rightData_, false);
    //delete[] rightData_;
}

bool BallTree::storeTree(const char *index_path) {
    char indexEntryPath[L], dataEntryPath[L];
    sprintf(indexEntryPath, "%s/indexEntries.dat", index_path);
    sprintf(dataEntryPath, "%s/dataEntries.dat", index_path);
    std::ofstream indexOutput(indexEntryPath, std::ios::out | std::ios::binary);
    std::ofstream dataOutput(dataEntryPath, std::ios::out | std::ios::binary);

    if (!indexOutput || !dataOutput)
        return false;

    /* Store some metadata */
    indexOutput.write((char*)&INDEX_SLOTSIZE, sizeof(int));
    indexOutput.write((char*)&dimension, sizeof(int));
    dataOutput.write((char*)&DATA_SLOTSIZE, sizeof(int));

    /* Preorder store the nodes. */
    curIndexRid = curDataRid = Rid(0, 0);
    bfsStore(indexOutput, dataOutput);
	//删除中间节点
	deleteIndexNode();
    //preorderStore(root, NULL, indexOutput, dataOutput, false);
    indexOutput.close();
    dataOutput.close();

    return true;
}

//用于删除节点空间占用的内存
bool BallTree::deleteIndexNode() {
	preDeleteIndexNode(root);
	return true;
}

void BallTree::preDeleteIndexNode(ballTreeNode* cur) {
	if (cur->left != NULL) preDeleteIndexNode(cur->left);
	if (cur->right != NULL) preDeleteIndexNode(cur->right);
	if (cur->mean != NULL) delete[]cur->mean;
	delete cur;
	cur = NULL;
}

void BallTree::bfsStore(std::ofstream &indexOutput, std::ofstream &dataOutput) {
    std::queue<ballTreeNode*> bfs;
    bfs.push(root);
    while (bfs.size()) {
        ballTreeNode *cur = bfs.front();
        if (cur->table == NULL) {
            if (cur->father != NULL) {
                if (cur->isLeft) cur->father->leftRid = curIndexRid;
                else cur->father->rightRid = curIndexRid;
                updateRid(cur->father, indexOutput);//更新父节点左右节点的rid
            }
            cur->myRid = curIndexRid;
            storeIndexNode(cur, indexOutput, curIndexRid);
            curIndexRid.pageid = curIndexRid.slotid == numIndexSlot - 1 ? curIndexRid.pageid + 1 : curIndexRid.pageid;
            curIndexRid.slotid = curIndexRid.slotid == numIndexSlot - 1 ? 0 : curIndexRid.slotid + 1;
        } else {
            if (cur->father != NULL) {
                if (cur->isLeft) cur->father->leftRid = curDataRid;
                cur->father->rightRid = curDataRid;
                updateRid(cur->father, indexOutput);//更新父节点左右节点的rid
            }
            cur->myRid = curDataRid;
            storeDataNode(cur, dataOutput, curDataRid);
			//cur->delete_memory();	//删除叶子节点上的table数据
            curDataRid.pageid = curDataRid.slotid == numDataSlot - 1 ? curDataRid.pageid + 1 : curDataRid.pageid;
            curDataRid.slotid = curDataRid.slotid == numDataSlot - 1 ? 0 : curDataRid.slotid + 1;
        }
        bfs.pop();
        if (cur->left != NULL)
            bfs.push(cur->left);
        if (cur->right != NULL)
            bfs.push(cur->right);
    }
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
        output.seekp(METADATA_INDEX_OFFSET + PAGE_SIZE * pageid);
        bool *bitMap = new bool[numIndexSlot];
        memset(bitMap, 0, numIndexSlot);
        output.write((char*)bitMap, numIndexSlot);
        delete[] bitMap;
    }
    output.seekp(METADATA_INDEX_OFFSET + pageid * PAGE_SIZE + slotid);
    bool a = true;
    output.write((char*)&a, sizeof(bool));

    output.seekp(METADATA_INDEX_OFFSET + pageid * PAGE_SIZE + slotid * INDEX_SLOTSIZE + numIndexSlot);
    float *floatArr = new float[dimension];
    int *intArr = new int[INDEX_INT_SIZE];
    intArr[0] = node->leftRid.pageid;
    intArr[1] = node->leftRid.slotid;
    intArr[2] = node->rightRid.pageid;
    intArr[3] = node->rightRid.slotid;
    floatArr[0] = node->radius;
    memcpy(floatArr + 1, node->mean, (dimension - 1) * sizeof(float));
    bool *boolArr = new bool[2];
    boolArr[0] = node->isLeftLeaf;
    boolArr[1] = node->isRightLeaf;
    output.write((char*)intArr, INDEX_INT_SIZE * sizeof(int));
    output.write((char*)floatArr, dimension * sizeof(float));
    output.write((char*)boolArr, INDEX_BOOL_SIZE * sizeof(bool));

    delete[] floatArr;
    delete[] intArr;
}
int dataNode_num = 0;
void BallTree::storeDataNode(ballTreeNode *node, std::ofstream &output, Rid &rid) {
    int pageid = rid.pageid;
    int slotid = rid.slotid;

    if (slotid == 0) {    // bitMap needs to be inserted.
        output.seekp(METADATA_DATA_OFFSET + PAGE_SIZE * pageid);
        bool *bitMap = new bool[numDataSlot];
        memset(bitMap, 0, numDataSlot);
        output.write((char*)bitMap, numDataSlot);
        delete[] bitMap;
    }
    output.seekp(METADATA_DATA_OFFSET + pageid * PAGE_SIZE + slotid);
    bool a = true;
    output.write((char*)&a, sizeof(bool));

    output.seekp(METADATA_DATA_OFFSET + pageid * PAGE_SIZE + slotid * DATA_SLOTSIZE + numDataSlot);
    int numTuples = node->tableSize;
    float *floatArr = new float[dimension * node->tableSize];
	float *mean = new float[dimension - 1];
	float radius;
	for (int i = 0; i < node->tableSize; i++) {
		memcpy(floatArr + i * dimension, node->table[i], dimension * sizeof(float));
	}
	memcpy(mean, node->mean, (dimension - 1) * sizeof(float));
	radius = node->radius;
    output.write((char*)&numTuples, sizeof(int));
	output.write((char*)mean, (dimension - 1) * sizeof(float));
	output.write((char*)&radius, sizeof(float));
    output.write((char*)floatArr, dimension * node->tableSize * sizeof(float));
	delete[] mean;
    delete[] floatArr;
	dataNode_num++;
	printf("Store %d leaf nodes\n", dataNode_num);
}

void BallTree::updateRid(ballTreeNode *node, std::ofstream &output) {
    output.seekp(METADATA_INDEX_OFFSET + PAGE_SIZE * node->myRid.pageid + numIndexSlot + INDEX_SLOTSIZE * node->myRid.slotid);
    int buf[INDEX_INT_SIZE];
    buf[0] = node->leftRid.pageid;
    buf[1] = node->leftRid.slotid;
    buf[2] = node->rightRid.pageid;
    buf[3] = node->rightRid.slotid;
    output.write((char*)buf, sizeof(int) * 4);
}

void BallTree::computeMean(float *&mean, int n, int d, float **data) {
    for (int i = 1; i < d; i++) {
        float tempSum = 0;
        for (int j = 0; j < n; j++) {
            tempSum += data[j][i];
        }
        mean[i - 1] = tempSum / n;
    }
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
	float *pick = data[n/2];
    A = pick;
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
    //restore()函数将树的根节点设置即可，一次找一个节点
    indexEntry_path = new char[L];
    dataEntry_path = new char[L];
    sprintf(indexEntry_path, "%s/indexEntries.dat", index_path);
    sprintf(dataEntry_path, "%s/dataEntries.dat", index_path);
    
    std::ifstream indexFile(indexEntry_path, std::ios::binary);
    std::ifstream dataFile(dataEntry_path, std::ios::binary);
    if (!indexFile || !dataFile) {
        printf("open failed!\n");
        return false;
    }

    int buffer;
    indexFile.read((char*)&buffer, sizeof(int));
    INDEX_SLOTSIZE = buffer;
    numIndexSlot = floor(PAGE_SIZE / (INDEX_SLOTSIZE + 1));
    indexFile.read((char*)&buffer, sizeof(int));
    dimension = buffer;
    dataFile.read((char*)&buffer, sizeof(int));
    DATA_SLOTSIZE = buffer;
    numDataSlot = floor(PAGE_SIZE / (DATA_SLOTSIZE + 1));

    indexFile.close();
    dataFile.close();

    root = getNode(0, 0, true);

	//add to debug
	buildNode(root);
    return true;
}
void BallTree::buildNode(ballTreeNode* node) {
	//if (node->table != NULL) return;
	if (!node->isLeftLeaf) {
		node->left = getNode(node->leftRid.pageid, node->leftRid.slotid, !(node->isLeftLeaf));
		buildNode(node->left);
	}
	else node->left = NULL;
	if (!node->isRightLeaf) {
		node->right = getNode(node->rightRid.pageid, node->rightRid.slotid, !(node->isRightLeaf));
		buildNode(node->right);
	}
	else node->right = NULL;
}

int BallTree::mipSearch(int d, float* query, Mip &mip) {
	//论文的算法5
	dimension = d;

	//记录最大内积及目标序号
	TreeSearch(query, root, mip);
	return mip.index;
}
void BallTree::TreeSearch(float* query, ballTreeNode* node, Mip &mip) {
	if (node->table != NULL) {
		LinearSearch(query, node, mip);
	}
	else {
		if (node->isLeftLeaf)
			node->left = getNode(node->leftRid.pageid, node->leftRid.slotid, !(node->isLeftLeaf));
		if (node->isRightLeaf)
			node->right = getNode(node->rightRid.pageid, node->rightRid.slotid, !(node->isRightLeaf));
		float leftProduct = MIP(query, node->left);
		float rightProduct = MIP(query, node->right);
		if (leftProduct < rightProduct) {
			if (mip.product < rightProduct) {
				TreeSearch(query, node->right, mip);
			}
			if (node->isRightLeaf&&node->right!=NULL) {
				node->right->delete_memory();
				delete []node->right->mean;
				delete node->right;
				node->right = NULL;
			}

			if (mip.product < leftProduct) {
				TreeSearch(query, node->left, mip);
			}
			if (node->isLeftLeaf&&node->left != NULL) {
				node->left->delete_memory();
				delete[]node->left->mean;
				delete node->left;
				node->left = NULL;
			}
		}
		else {
			if (mip.product < leftProduct) {
				TreeSearch(query, node->left, mip);
			}
			if (node->isLeftLeaf&&node->left!=NULL) {
				node->left->delete_memory();
				delete[]node->left->mean;
				delete node->left;
				node->left = NULL;
			}

			if (mip.product < rightProduct&&node->right!=NULL) {
				TreeSearch(query, node->right, mip);
			}
			if (node->isRightLeaf) {
				node->right->delete_memory();
				delete[]node->right->mean;
				delete node->right;
				node->right = NULL;
			}
		}

	}
}


void BallTree::LinearSearch(float* query, ballTreeNode* node, Mip& mip) {
//对存有数据的叶子节点进行线性查找
    for (int i = 0; i < node->tableSize; i++) {
        float newProduct = computeInnerProduct(query, node->table[i]);
        if (newProduct > mip.product) {
            mip.index = node->table[i][0];
            //这里的index位每一个数据项的第一个数
            mip.product = newProduct;
        }
    }
}

float BallTree::computeInnerProduct(float* query, float* data) {
//计算内积
    float product = 0;
    for (int i = 0; i < dimension - 1; i++) {
        product += query[i+1] * data[i + 1];
    }
    return product;
}

float BallTree::MIP(float *query, ballTreeNode* node) {
//估算球的上界与query的内积
    float product = 0;
	float radiusQuery = 0;
    for (int i = 0; i < dimension - 1; i++) {
        product += query[i+1] * node->mean[i];
		radiusQuery += query[i+1] * query[i+1];
    }
    product += node->radius*sqrt(radiusQuery);
    return product;
}


ballTreeNode* BallTree::getNode(int pageID, int slotID, bool isIndex) {
//根据页号和槽号获得节点数据
    ballTreeNode* node = new ballTreeNode();
    char *file_path = isIndex ? indexEntry_path : dataEntry_path;
    std::ifstream input(file_path, std::ios::binary);
    int offset = isIndex ? METADATA_INDEX_OFFSET : METADATA_DATA_OFFSET;
    int numSlot = isIndex ? numIndexSlot : numDataSlot;
    int slotsize = isIndex ? INDEX_SLOTSIZE : DATA_SLOTSIZE;
    input.seekg(offset + pageID * PAGE_SIZE + numSlot + slotsize * slotID);
    if (isIndex) {
        int intBuffer[INDEX_INT_SIZE];
        float *floatBuffer = new float[dimension];
        bool boolBuffer[INDEX_BOOL_SIZE];
        input.read((char*)intBuffer, sizeof(int) * INDEX_INT_SIZE);
        input.read((char*)floatBuffer, sizeof(float) * (dimension));
        input.read((char*)boolBuffer, sizeof(bool) * INDEX_BOOL_SIZE);
        float *mean = new float[dimension - 1];
        memcpy(mean, floatBuffer + 1, sizeof(float) * (dimension - 1));
        node->myRid = Rid(pageID, slotID);
        node->radius = floatBuffer[0];
        node->mean = mean;
        node->leftRid = Rid(intBuffer[0], intBuffer[1]);
        node->rightRid = Rid(intBuffer[2], intBuffer[3]);
        node->isLeftLeaf = boolBuffer[0];
        node->isRightLeaf = boolBuffer[1];
        
        delete[] floatBuffer;
    } else {
        int intBuffer;
        input.read((char*)&intBuffer, sizeof(int));
        node->table = new float*[intBuffer];
		node->mean = new float[dimension - 1];
        node->tableSize = intBuffer;
        float *floatBuffer = new float[dimension];
		float *mean = new float[dimension - 1];
		float radius;
		input.read((char*)mean, sizeof(float) * (dimension - 1));
		input.read((char*)&radius, sizeof(float));

		memcpy(node->mean, mean, sizeof(float) * (dimension - 1));
		node->radius = radius;

        for (int i = 0; i < intBuffer; i++) {
            node->table[i] = new float[dimension];
            input.read((char*)floatBuffer, sizeof(float) * dimension);
            memcpy(node->table[i], floatBuffer, sizeof(float) * dimension);
        }
        node->myRid = Rid(pageID, slotID);
        delete[] floatBuffer;
        delete[] mean;
    }
    return node;
}
