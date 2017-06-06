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
    INDEX_SLOTSIZE = sizeof(int) * INDEX_INT_SIZE + sizeof(float) * dimension + sizeof(bool) * INDEX_BOOL_SIZE;
	//子节点的页号&槽号 + 半径&圆心（半径+元组平均值）+ 对子节点是否为子叶节点的bool判断
	//半径占平均值数组的第一个位置，后面的位置都是各维度的平均值，所以用到dimension的大小，其实一个数据是dimension的大小，第一个位置是index
    DATA_SLOTSIZE = sizeof(int) * DATA_INT_SIZE + sizeof(float) * dimension * N0 + sizeof(float) * dimension;
	//叶子节点的tableSize数据个数 + 对子节点是否为子叶节点的bool判断 + 半径&圆心（半径+元组平均值）
    numIndexSlot = floor(PAGE_SIZE / (INDEX_SLOTSIZE + 1));
    numDataSlot = floor(PAGE_SIZE / (DATA_SLOTSIZE + 1));

    buildBall(root, NULL, n, d, data, false, false);

    if (root == NULL)
        return false;

    return true;
}

int leaf_num = 0;
void BallTree::buildBall(ballTreeNode *&node, ballTreeNode *father, int n, int d, float **data, bool dir1, bool dir2) {
    float *mean = new float[d - 1];
    computeMean(mean, n, d, data);
    node = new ballTreeNode(computeRadius(n, d, data, mean), mean, d);
	delete[] mean;
    node->father = father;
	if (dir1) {
		if (dir2)
			node->isLeft = true;
		else 
			node->isRight = true;
	}
	else {
		if (dir2)
			node->isUp = true;
		else
			node->isDown = true;
	}


    if (n <= N0) {
        node->table = new float*[n];
        node->tableSize = n;
        for (int i = 0; i < n; i++) {
			node->table[i] = data[i];
			// 通过修改注释部分减少了内存消耗
            // node->table[i] = new float[d];
            // memcpy(node->table[i], data[i], d * sizeof(float));
        }
		if (node != root) delete[]data;
		leaf_num++;
		printf("Make %d leaf nodes.\n", leaf_num);
        return;
    }

    float *A, *B, *C, *D;
    MakeBallTreeSplit(n, d, data, A, B);

    std::vector<float*> leftData, rightData, upData, downData;
    for (int i = 0; i < n; i++) {
        if (computeDistance(data[i], A) <= computeDistance(data[i], B)) {
            leftData.push_back(data[i]);
        } else {
            rightData.push_back(data[i]);
        }
    }

	int lsize = leftData.size();
	int rsize = rightData.size();
    float **leftData_ = parseFloatArr(leftData);
	leftData.swap(std::vector<float*>());
    float **rightData_ = parseFloatArr(rightData);
	rightData.swap(std::vector<float*>());
	MakeBallTreeSplit(lsize, d, leftData_, A, B);
	MakeBallTreeSplit(rsize, d, rightData_, C, D);
	delete[] leftData_;
	delete[] rightData_;
	
	for (int i = 0; i < n; i++) {
		float ADist, BDist, CDist, DDist;
		ADist = computeDistance(data[i], A);
		BDist = computeDistance(data[i], B);
		CDist = computeDistance(data[i], C);
		DDist = computeDistance(data[i], D);
		if (ADist<BDist) {
			if (ADist < CDist) {
				if (ADist < DDist) {
					//A最小
					leftData.push_back(data[i]);
				}
				else {
					//D最小
					downData.push_back(data[i]);
				}
			}
			else {
				if (CDist < DDist) {
					//C最小
					upData.push_back(data[i]);
				}
				else {
					//D最小
					downData.push_back(data[i]);
				}
			}
		}
		else {
			if (BDist < CDist) {
				if (BDist < DDist) {
					//B最小
					rightData.push_back(data[i]);
				}
				else {
					//D最小
					downData.push_back(data[i]);
				}
			}
			else {
				if (CDist < DDist) {
					//C最小
					upData.push_back(data[i]);
				}
				else {
					//D最小
					downData.push_back(data[i]);
				}
			}
		}
	}
	if (node != root) delete[]data;
	//将vector变为float数组，同时释放vector内存
	int leftn = leftData.size();
	int rightn = rightData.size();
	int upn = upData.size();
	int downn = downData.size();
	leftData_ = parseFloatArr(leftData);
	leftData.swap(std::vector<float*>());
	rightData_ = parseFloatArr(rightData);
	rightData.swap(std::vector<float*>());
	float **upData_ = parseFloatArr(upData);
	upData.swap(std::vector<float*>());
	float **downData_ = parseFloatArr(downData);
	downData.swap(std::vector<float*>());
	
	

    node->isLeftLeaf = leftn <= N0;
    node->isRightLeaf = rightn <= N0;
	node->isUpLeaf = upn <= N0;
	node->isDownLeaf = downn <= N0;
    buildBall(node->left, node, leftn, d, leftData_, true, true);
    buildBall(node->right, node, rightn, d, rightData_, true, false);
	buildBall(node->up, node, upn, d, upData_, false, true);
	buildBall(node->down, node, downn, d, downData_, false, false);
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
	//四叉树扩展部分
	if (cur->up != NULL) preDeleteIndexNode(cur->up);
	if (cur->down != NULL) preDeleteIndexNode(cur->down);
	//
	if (cur->mean != NULL) delete[]cur->mean;
	delete cur;
	cur = NULL;
}

void BallTree::bfsStore(std::ofstream &indexOutput, std::ofstream &dataOutput) {
    std::queue<ballTreeNode*> bfs;
    bfs.push(root);
    while (bfs.size()) {
        ballTreeNode *cur = bfs.front();
        if (cur->table == NULL) {	//当前节点为叶子节点
            if (cur->father != NULL) {
                if (cur->isLeft) cur->father->leftRid = curIndexRid;
                else if (cur->isRight) cur->father->rightRid = curIndexRid;
				else if (cur->isUp) cur->father->upRid = curIndexRid;
				else cur->father->downRid = curIndexRid;
                updateRid(cur->father, indexOutput);//更新父节点子节点的rid
            }
            cur->myRid = curIndexRid;
            storeIndexNode(cur, indexOutput, curIndexRid);
			//获得下个节点的Rid
            curIndexRid.pageid = curIndexRid.slotid == numIndexSlot - 1 ? curIndexRid.pageid + 1 : curIndexRid.pageid;
            curIndexRid.slotid = curIndexRid.slotid == numIndexSlot - 1 ? 0 : curIndexRid.slotid + 1;
        } else {
			if (cur->father != NULL) {
				if (cur->isLeft) cur->father->leftRid = curDataRid;
				else if (cur->isRight) cur->father->rightRid = curDataRid;
				else if (cur->isUp) cur->father->upRid = curDataRid;
				else cur->father->downRid = curDataRid;
                updateRid(cur->father, indexOutput);//更新父节点左右节点的rid
            }
            cur->myRid = curDataRid;
            storeDataNode(cur, dataOutput, curDataRid);
			//cur->delete_memory();	//删除叶子节点上的table数据
			//获得下个节点的Rid
            curDataRid.pageid = curDataRid.slotid == numDataSlot - 1 ? curDataRid.pageid + 1 : curDataRid.pageid;
            curDataRid.slotid = curDataRid.slotid == numDataSlot - 1 ? 0 : curDataRid.slotid + 1;
        }
        bfs.pop();
        if (cur->left != NULL)
            bfs.push(cur->left);
        if (cur->right != NULL)
            bfs.push(cur->right);
		//四叉树扩展部分
		if (cur->up != NULL)
			bfs.push(cur->up);
		if (cur->down != NULL)
			bfs.push(cur->down);
		//
    }
}
//
//void BallTree::preorderStore(ballTreeNode *node, ballTreeNode *father, std::ofstream &indexOutput, 
//    std::ofstream &dataOutput, bool isLeft) {
//    if (node == NULL)
//        return;
//    if (node->table == NULL) {
//        if (father != NULL) {
//            if (isLeft) father->leftRid = curIndexRid;
//            else father->rightRid = curIndexRid;
//        }
//        storeIndexNode(node, indexOutput, curIndexRid);
//        curIndexRid.pageid = curIndexRid.slotid == numIndexSlot - 1 ? curIndexRid.pageid + 1 : curIndexRid.pageid;
//        curIndexRid.slotid = curIndexRid.slotid == numIndexSlot - 1 ? 0 : curIndexRid.slotid + 1;
//    } else {
//        if (father != NULL) {
//            if (isLeft) father->leftRid = curDataRid;
//            else father->rightRid = curDataRid;
//        }
//        curDataRid.pageid = curDataRid.slotid == numDataSlot - 1 ? curDataRid.pageid + 1 : curDataRid.pageid;
//        curDataRid.slotid = curDataRid.slotid == numDataSlot - 1 ? 0 : curDataRid.slotid + 1;
//        storeDataNode(node, dataOutput, curDataRid);
//    }
//
//    preorderStore(node->left, node, indexOutput, dataOutput, true);
//    preorderStore(node->right, node, indexOutput, dataOutput, false);
//}

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
	//四叉树扩展内容
	intArr[4] = node->upRid.pageid;
	intArr[5] = node->upRid.slotid;
	intArr[6] = node->downRid.pageid;
	intArr[7] = node->downRid.slotid;
	//
    floatArr[0] = node->radius;
    memcpy(floatArr + 1, node->mean, (dimension - 1) * sizeof(float));
    bool *boolArr = new bool[INDEX_BOOL_SIZE];
    boolArr[0] = node->isLeftLeaf;
    boolArr[1] = node->isRightLeaf;
	//四叉树扩展内容
	boolArr[2] = node->isUpLeaf;
	boolArr[3] = node->isDownLeaf;
	//
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
	dataNode_num++;
	printf("Store %d leaf nodes\n", dataNode_num);
	delete[] mean;
    delete[] floatArr;
}

void BallTree::updateRid(ballTreeNode *node, std::ofstream &output) {
    output.seekp(METADATA_INDEX_OFFSET + PAGE_SIZE * node->myRid.pageid + numIndexSlot + INDEX_SLOTSIZE * node->myRid.slotid);
    int buf[INDEX_INT_SIZE];
    buf[0] = node->leftRid.pageid;
    buf[1] = node->leftRid.slotid;
    buf[2] = node->rightRid.pageid;
    buf[3] = node->rightRid.slotid;
	//四叉树扩展
	buf[4] = node->upRid.pageid;
	buf[5] = node->upRid.slotid;
	buf[6] = node->downRid.pageid;
	buf[7] = node->downRid.slotid;
	//
    output.write((char*)buf, sizeof(int) * INDEX_INT_SIZE);
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
	else node->up = NULL;
	if (!node->isUpLeaf) {
		node->up = getNode(node->upRid.pageid, node->upRid.slotid, !(node->isUpLeaf));
		buildNode(node->up);
	}
	else node->down = NULL;
	if (!node->isDownLeaf) {
		node->down = getNode(node->downRid.pageid, node->downRid.slotid, !(node->isDownLeaf));
		buildNode(node->down);
	}
	else node->down = NULL;

}

int BallTree::mipSearch(int d, float* query) {
	//论文的算法5
	dimension = d;

	Mip mip;
	mip.product = -1;
	mip.index = -1;
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
		if (node->isUpLeaf)
			node->up = getNode(node->upRid.pageid, node->upRid.slotid, !(node->isUpLeaf));
		if (node->isDownLeaf)
			node->down = getNode(node->downRid.pageid, node->downRid.slotid, !(node->isDownLeaf));
		float leftProduct = MIP(query, node->left);
		float rightProduct = MIP(query, node->right);
		float upProduct = MIP(query, node->up);
		float downProduct = MIP(query, node->down);
		if (leftProduct > rightProduct) {
			if (mip.product < leftProduct&&node->left != NULL) {
				TreeSearch(query, node->left, mip);
			}
			node->release_left();
			if (mip.product < rightProduct&&node->right != NULL) {
				TreeSearch(query, node->right, mip);
			}
			node->release_right();
		}
		else {
			if (mip.product < rightProduct&&node->right != NULL) {
				TreeSearch(query, node->right, mip);
			}
			node->release_right();
			if (mip.product < leftProduct&&node->left != NULL) {
				TreeSearch(query, node->left, mip);
			}
			node->release_left();
		}

		if (upProduct > downProduct) {
			if (mip.product < upProduct&&node->up != NULL) {
				TreeSearch(query, node->up, mip);
			}
			node->release_up();
			if (mip.product < downProduct&&node->down != NULL) {
				TreeSearch(query, node->down, mip);
			}
			node->release_down();
		}
		else {
			if (mip.product < downProduct&&node->down != NULL) {
				TreeSearch(query, node->down, mip);
			}
			node->release_down();
			if (mip.product < upProduct&&node->up != NULL) {
				TreeSearch(query, node->up, mip);
			}
			node->release_up();
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
		node->upRid = Rid(intBuffer[4], intBuffer[5]);
		node->downRid = Rid(intBuffer[6], intBuffer[7]);
        node->isLeftLeaf = boolBuffer[0];
        node->isRightLeaf = boolBuffer[1];
		node->isUpLeaf = boolBuffer[2];
		node->isDownLeaf = boolBuffer[3];
        
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
