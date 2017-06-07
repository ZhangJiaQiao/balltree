#include <cmath>
#include <fstream>
#include <vector>
#include <queue>
#include <cstdio>
#include <cstring>
#include<string.h>

#include "BallTree.h"
BallTree::BallTree() {
	ModelQuad = false;
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
/*
 *buildTree
 *@param: n, d, data
 *@return: bool
 *用于初始化建树参数等， 如dimension
 *其中的REAL_DATA_INT_SIZE,REAL_INDEX_INT_SIZE,REAL_INDEX_BOOL_SIZE等
 *都是用于转换而二叉树和四叉树不同文件存取的参数接口
 *buildBall 是实际实现递归建树的函数
 */
bool BallTree::buildTree(int n, int d, float **data) {
	/* Initialize some values at buiding time. */
	ModelQuad = false;
	//ModelQuad=false代表非四叉树模式
	updateSize();
	dimension = d;
	INDEX_SLOTSIZE = sizeof(int) * REAL_INDEX_INT_SIZE + sizeof(float) * dimension + sizeof(bool) * REAL_INDEX_BOOL_SIZE;
	DATA_SLOTSIZE = sizeof(int) * REAL_DATA_INT_SIZE + sizeof(float) * dimension * N0 + sizeof(float) * dimension;
	numIndexSlot = floor(PAGE_SIZE / (INDEX_SLOTSIZE + 1));
	numDataSlot = floor(PAGE_SIZE / (DATA_SLOTSIZE + 1));

	buildBall(root, NULL, n, d, data, false);

	if (root == NULL)
		return false;

	return true;
}
/*
 *buildBall
 *@param: node, father, n, d, data, isLeft
 *@return: void
 *建树算法核心部分
 *传入buildBall的data是上一个父节点的数据集
 *其实在本次实验中并没有按照论文里的吧每个节点包含的数据存在该节点
 *
 *具体实现：
 *传入节点的data，计算得到半径和圆心
 *以半径圆心构建节点
 *将data通过MakeBallTreeSplit计算获得两个点
 *通过对两个点距离的比较将data分离得到两个集合
 *再递归传入buildBall进行子节点的构建
 */
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
		if (father != NULL) delete[]data;
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
		}
		else {
			rightData.push_back(data[i]);
		}
	}
	if (father != NULL) delete[]data;

	int leftn = leftData.size();
	int rightn = rightData.size();
	//释放vector里的内存
	std::vector<float*> empty;
	float **leftData_ = parseFloatArr(leftData);
	leftData.swap(empty);
	float **rightData_ = parseFloatArr(rightData);
	rightData.swap(empty);

	node->isLeftLeaf = leftn <= N0;
	node->isRightLeaf = rightn <= N0;
	buildBall(node->left, node, leftn, d, leftData_, true);
	//delete[] leftData_;
	buildBall(node->right, node, rightn, d, rightData_, false);
	//delete[] rightData_;
}
/*
 *storeTree
 *@param: index_path
 *
 *存储节点数据
 *分为两个部分存入两个二进制文件
 *
 *文件结构：元数据+多个页
 *页结构：位图+多个槽
 *
 *打开两个输入流，存入文件的元数据，描述是否为四叉树模式以及文件槽大小
 *通过bfsStore进行广度优先的顺序存储整个树的节点到文件中
 *最后删除树留在内存中的节点
 *
 */
bool BallTree::storeTree(const char *index_path) {
	char indexEntryPath[L], dataEntryPath[L];
	sprintf(indexEntryPath, "%s/indexEntries.dat", index_path);
	sprintf(dataEntryPath, "%s/dataEntries.dat", index_path);
	std::ofstream indexOutput(indexEntryPath, std::ios::out | std::ios::binary);
	std::ofstream dataOutput(dataEntryPath, std::ios::out | std::ios::binary);

	if (!indexOutput || !dataOutput)
		return false;

	/* Store some metadata */
	indexOutput.write((char*)&ModelQuad, sizeof(bool));
	indexOutput.write((char*)&INDEX_SLOTSIZE, sizeof(int));
	indexOutput.write((char*)&dimension, sizeof(int));
	dataOutput.write((char*)&ModelQuad, sizeof(bool));
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
/*
 *preDeleteIndexNode
 *@param: cur
 *
 *前序遍历删除节点
 */
void BallTree::preDeleteIndexNode(ballTreeNode* cur) {
	if (cur->left != NULL) preDeleteIndexNode(cur->left);
	if (cur->right != NULL) preDeleteIndexNode(cur->right);
	if (ModelQuad) {
		//四叉树扩展部分
		if (cur->up != NULL) preDeleteIndexNode(cur->up);
		if (cur->down != NULL) preDeleteIndexNode(cur->down);
		//
	}
	if (cur->mean != NULL) delete[]cur->mean;
	delete cur;
	cur = NULL;
}

/*
 *bfsStore
 *@param: indexOutput, dataOutput
 *
 *获得两个文件指针用于写入中间结点和叶子节点
 *
 *广度优先遍历各个节点
 *通过curRid, curIndexRid, curDataRid实时赋予当前节点存入位置的Rid
 *
 *第一部分：中间节点index
 *调用storeIndexNode()存储数据
 *
 *第二部分：叶子节点data
 *调用storeDataNode()存储数据
 *
 *在对节点进行存储是需要对其父节点对应的存储位置的子节点的Rid进行更新：updateRid()
 *
 */
void BallTree::bfsStore(std::ofstream &indexOutput, std::ofstream &dataOutput) {
	std::queue<ballTreeNode*> bfs;
	bfs.push(root);
	while (bfs.size()) {
		ballTreeNode *cur = bfs.front();
		if (cur->table == NULL) {
			if (cur->father != NULL) {
				if (ModelQuad) {
					if (cur->isLeft) cur->father->leftRid = curIndexRid;
					else if (cur->isRight) cur->father->rightRid = curIndexRid;
					else if (cur->isUp) cur->father->upRid = curIndexRid;
					else cur->father->downRid = curIndexRid;
				}
				else {
					if (cur->isLeft) cur->father->leftRid = curIndexRid;
					else cur->father->rightRid = curIndexRid;
				}
				updateRid(cur->father, indexOutput);//更新父节点左右节点的rid
			}
			cur->myRid = curIndexRid;
			storeIndexNode(cur, indexOutput, curIndexRid);
			curIndexRid.pageid = curIndexRid.slotid == numIndexSlot - 1 ? curIndexRid.pageid + 1 : curIndexRid.pageid;
			curIndexRid.slotid = curIndexRid.slotid == numIndexSlot - 1 ? 0 : curIndexRid.slotid + 1;
		}
		else {
			if (cur->father != NULL) {
				if (ModelQuad) {
					if (cur->isLeft) cur->father->leftRid = curDataRid;
					else if (cur->isRight) cur->father->rightRid = curDataRid;
					else if (cur->isUp) cur->father->upRid = curDataRid;
					else cur->father->downRid = curDataRid;
				}
				else {
					if (cur->isLeft) cur->father->leftRid = curDataRid;
					else cur->father->rightRid = curDataRid;
				}
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
		if (ModelQuad) {
			//四叉树扩展部分
			if (cur->up != NULL)
				bfs.push(cur->up);
			if (cur->down != NULL)
				bfs.push(cur->down);
			//
		}
	}
}

/*
 *storeIndexNode
 *@param： node， output， rid
 *
 *对于第一个槽的情况，先在开头存入该页的bitmap位图
 *槽里存储节点的半径圆心+左右节点的页号和槽号+左右节点是否为叶子节点的标记
 *当为四叉树时加入up， down节点的页号和槽号以及标记
 *
 */
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
	int *intArr = new int[REAL_INDEX_INT_SIZE];
	intArr[0] = node->leftRid.pageid;
	intArr[1] = node->leftRid.slotid;
	intArr[2] = node->rightRid.pageid;
	intArr[3] = node->rightRid.slotid;
	if (ModelQuad) {
		//四叉树扩展内容
		intArr[4] = node->upRid.pageid;
		intArr[5] = node->upRid.slotid;
		intArr[6] = node->downRid.pageid;
		intArr[7] = node->downRid.slotid;
		//
	}
	floatArr[0] = node->radius;
	memcpy(floatArr + 1, node->mean, (dimension - 1) * sizeof(float));
	bool *boolArr = new bool[REAL_INDEX_BOOL_SIZE];
	boolArr[0] = node->isLeftLeaf;
	boolArr[1] = node->isRightLeaf;
	if (ModelQuad) {
		//四叉树扩展内容
		boolArr[2] = node->isUpLeaf;
		boolArr[3] = node->isDownLeaf;
		//
	}
	output.write((char*)intArr, REAL_INDEX_INT_SIZE * sizeof(int));
	output.write((char*)floatArr, dimension * sizeof(float));
	output.write((char*)boolArr, REAL_INDEX_BOOL_SIZE * sizeof(bool));

	delete[] floatArr;
	delete[] intArr;
}
/*
 *storeDataNode
 *@param: node, output, rid
 *
 *槽里存储叶子节点的半径圆心+元组大小+元组数据
 *
 */
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
/*
 *updateRid
 *@param: node, output
 *
 *写入当前节点对于从槽里的左右节点的rid
 *当为四叉时则加入up， down节点的更新
 *
 */
void BallTree::updateRid(ballTreeNode *node, std::ofstream &output) {
	output.seekp(METADATA_INDEX_OFFSET + PAGE_SIZE * node->myRid.pageid + numIndexSlot + INDEX_SLOTSIZE * node->myRid.slotid);
	int *buf=new int[REAL_INDEX_INT_SIZE];
	buf[0] = node->leftRid.pageid;
	buf[1] = node->leftRid.slotid;
	buf[2] = node->rightRid.pageid;
	buf[3] = node->rightRid.slotid;
	if (ModelQuad) {
		//四叉树扩展
		buf[4] = node->upRid.pageid;
		buf[5] = node->upRid.slotid;
		buf[6] = node->downRid.pageid;
		buf[7] = node->downRid.slotid;
		//
	}
	output.write((char*)buf, sizeof(int) * REAL_INDEX_INT_SIZE);
	delete[]buf;
}
/*
 *computeMean
 *@param: mean, n, d, data
 *
 *计算data数据里的平均值
 *同时存储到mean里返回
 *
 */
void BallTree::computeMean(float *&mean, int n, int d, float **data) {
	for (int i = 1; i < d; i++) {
		float tempSum = 0;
		for (int j = 0; j < n; j++) {
			tempSum += data[j][i];
		}
		mean[i - 1] = tempSum / n;
	}
}
/*
 *computeRadius
 *@param： n, d, data, mean
 *@return: max
 *
 *计算data中与mean最远的点的距离，作为半径
 *返回半径值
 *
 */
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
/*
 *computeDistance
 *@param: x, y
 *@return: distance
 *
 *计算d维里x与y的距离
 *返回距离值
 *
 *注：从1开始计算的原因是元组数据里第一个位置是index，即数据的下标
 *因此需要忽略0号位置
 *
 */
float BallTree::computeDistance(float *x, float *y) {
	float squareSum = 0;
	for (int i = 1; i < dimension; i++)
		squareSum += (x[i] - y[i]) * (x[i] - y[i]);
	return sqrt(squareSum);
}
/*
 *MakeBallTreeSplit
 *@param: n, d, data, A, B
 *
 *取data中间位置的点
 *寻找其对应的最远的点记作A
 *寻找A点对应最远的点记作B
 *
 */
bool BallTree::MakeBallTreeSplit(int n, int d, float **data, float *&A, float *&B) {
	float *pick = data[n / 2];
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
/*
 *restoreTree
 *@param: index_path
 *
 *重新从文件载入树
 *这里先读取了元数据，获得槽大小以及相关信息
 *再getNode函数获得根节点root
 *root节点里有左右节点的rid以及是否为叶子节点的标记信息
 *用于后面buildNode的递归构建树
 *
 */
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
	indexFile.read((char*)&ModelQuad, sizeof(bool));
	indexFile.read((char*)&buffer, sizeof(int));
	INDEX_SLOTSIZE = buffer;
	numIndexSlot = floor(PAGE_SIZE / (INDEX_SLOTSIZE + 1));
	indexFile.read((char*)&buffer, sizeof(int));
	dataFile.read((char*)&ModelQuad, sizeof(bool));
	dimension = buffer;
	dataFile.read((char*)&buffer, sizeof(int));
	DATA_SLOTSIZE = buffer;
	numDataSlot = floor(PAGE_SIZE / (DATA_SLOTSIZE + 1));

	indexFile.close();
	dataFile.close();
	updateSize();

	root = getNode(0, 0, true);

	//add to debug
	buildNode(root);
	return true;
}
/*
 *buildNode
 *@param： node
 *
 *前序遍历，结合getNode从文件里获取节点值
 *对于叶子节点不获取，通过文件对树的重建只把树的中间节点载入内存
 *而不是整棵树的所有数据，同时叶子节点具有的数据如果载入内存必定会对内存消耗巨大
 *
 */
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
	if (ModelQuad) {
		if (!node->isUpLeaf) {
			node->up = getNode(node->upRid.pageid, node->upRid.slotid, !(node->isUpLeaf));
			buildNode(node->up);
		}
		else node->up = NULL;
		if (!node->isDownLeaf) {
			node->down = getNode(node->downRid.pageid, node->downRid.slotid, !(node->isDownLeaf));
			buildNode(node->down);
		}
		else node->down = NULL;
	}
}
/*
 *mipSearch
 *@param：d， query
 *
 *调用QuadTreeSearch或者TreeSearch
 *进行分支限界查找最大内积
 */
int BallTree::mipSearch(int d, float* query) {
	//论文的算法5
	dimension = d;

	Mip mip;
	mip.product = -1;
	mip.index = -1;
	//记录最大内积及目标序号
	if (ModelQuad) QuadTreeSearch(query, root, mip);
	else TreeSearch(query, root, mip);
	return mip.index;
}
/*
 *TreeSearch
 *@param: query, node, mip
 *
 *当node为叶子节点时，调用LinearSearch进行叶子节点的线性遍历查找最大内积
 *当为中间节点时，判断其左右节点是否为叶子节点，若是叶子节点则getNode将对应叶子节点的数据
 *从内存读入
 *
 *再调用MIP计算左右节点的内积上限
 *对MIP小的分支进行优先查找
 *对于MIP小于当前最大内积的分支不予查找
 *
 *最后做好内存的控制，将node中的堆内存delete清空
 *
 */
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
			if (node->isRightLeaf&&node->right != NULL) {
				node->right->delete_memory();
				delete[]node->right->mean;
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
			if (node->isLeftLeaf&&node->left != NULL) {
				node->left->delete_memory();
				delete[]node->left->mean;
				delete node->left;
				node->left = NULL;
			}

			if (mip.product < rightProduct&&node->right != NULL) {
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

/*
 *LinearSearch
 *@param: query, node, mip
 *
 *计算query与node节点中table里的元组进行内存计算
 *同时比较当前最大内积mip，以更新mip最大内积值，以及对应元组index
 */
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
/*
 *computeInnerProduct
 *@param: query, data
 *@return: product
 *
 *计算query与data之间的内积并返回
 */
float BallTree::computeInnerProduct(float* query, float* data) {
	//计算内积
	float product = 0;
	for (int i = 0; i < dimension - 1; i++) {
		product += query[i + 1] * data[i + 1];
	}
	return product;
}
/*
 *MIP
 *@param： query， node
 *@return： product
 *
 *这里是算法的核心部分
 *由于论文中证明了query与集合中元组内积最大上限的计算公式：
 *MAX=query·mean+|query|*radius
 *函数返回该上限
 *
 */
float BallTree::MIP(float *query, ballTreeNode* node) {
	//估算球的上界与query的内积
	float product = 0;
	float radiusQuery = 0;
	for (int i = 0; i < dimension - 1; i++) {
		product += query[i + 1] * node->mean[i];
		radiusQuery += query[i + 1] * query[i + 1];
	}
	product += node->radius*sqrt(radiusQuery);
	return product;
}

/*
 *getNode
 *@param： pageID， slotID， isIndex
 *@return： node
 *
 *这里是从index和data文件里重新载入节点的关键模块
 *通过传入的页号和槽号计算偏移值，从而获得所需数据对应槽的位置
 *
 *通过isIndex判断节点是否为叶子节点，因为叶子节点和中间节点的存储格式以及文件不同
 *
 */
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
		int *intBuffer=new int[REAL_INDEX_INT_SIZE];
		float *floatBuffer = new float[dimension];
		bool *boolBuffer=new bool[REAL_INDEX_BOOL_SIZE];
		input.read((char*)intBuffer, sizeof(int) * REAL_INDEX_INT_SIZE);
		input.read((char*)floatBuffer, sizeof(float) * (dimension));
		input.read((char*)boolBuffer, sizeof(bool) * REAL_INDEX_BOOL_SIZE);
		float *mean = new float[dimension - 1];
		memcpy(mean, floatBuffer + 1, sizeof(float) * (dimension - 1));
		node->myRid = Rid(pageID, slotID);
		node->radius = floatBuffer[0];
		node->mean = mean;
		node->leftRid = Rid(intBuffer[0], intBuffer[1]);
		node->rightRid = Rid(intBuffer[2], intBuffer[3]);
		node->isLeftLeaf = boolBuffer[0];
		node->isRightLeaf = boolBuffer[1];
		if (ModelQuad) {
			node->upRid = Rid(intBuffer[4], intBuffer[5]);
			node->downRid = Rid(intBuffer[6], intBuffer[7]);
			node->isUpLeaf = boolBuffer[2];
			node->isDownLeaf = boolBuffer[3];
		}

		delete[]boolBuffer;
		delete[]intBuffer;
		delete[] floatBuffer;
	}
	else {
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
//--------------------------------quad ball tree by zzy------------------------------
/*
 *updateSize
 *
 *由于四叉树和二叉树数据存取格式不同
 *因此需要更新REAL_INDEX_INT_SIZE, REAL_INDEX_BOOL_SIZE, REAL_DATA_INT_SIZE
 *确保二叉树和四叉树的正常转换
 *
 */
void BallTree::updateSize() {
	REAL_INDEX_INT_SIZE = (ModelQuad) ? QUAD_INDEX_INT_SIZE : INDEX_INT_SIZE;
	REAL_INDEX_BOOL_SIZE = (ModelQuad) ? QUAD_INDEX_BOOL_SIZE : INDEX_BOOL_SIZE;
	REAL_DATA_INT_SIZE = (ModelQuad) ? QUAD_DATA_INT_SIZE : DATA_INT_SIZE;

}
/*
 *buildQuadTree
 *@param: n, d, data
 *
 *四叉树构建
 *ModelQuad=true代表树模式是四叉树
 *同时调用updateSize()更新存取格式
 */
bool BallTree::buildQuadTree(int n, int d, float **data) {
	ModelQuad = true;
	updateSize();
	dimension = d;
	INDEX_SLOTSIZE = sizeof(int) * REAL_INDEX_INT_SIZE + sizeof(float) * dimension + sizeof(bool) * REAL_INDEX_BOOL_SIZE;
	DATA_SLOTSIZE = sizeof(int) * REAL_DATA_INT_SIZE + sizeof(float) * dimension * N0 + sizeof(float) * dimension;
	numIndexSlot = floor(PAGE_SIZE / (INDEX_SLOTSIZE + 1));
	numDataSlot = floor(PAGE_SIZE / (DATA_SLOTSIZE + 1));

	buildQuadBall(root, NULL, n, d, data, false, false);

	if (root == NULL)
		return false;

	return true;
}
/*
 *buildQuadBall
 *@param node, father, n, d, data, dir1, dir2
 *
 *类似于buildBall
 *构建节点信息半径圆心等
 *对于data数量小于N0的进行叶子节点的构建
 *
 *不同的是传入dir1和dir2判断四个方向的节点位置
 *同时调用MakeQuadBallTreeSplit的函数可以获取data里发布较为均匀的四个点
 *再对data按四个点进行划分
 *获得四个部分进行子节点的迭代
 *
 */
void BallTree::buildQuadBall(ballTreeNode *&node, ballTreeNode *father, int n, int d, float **data, bool dir1, bool dir2) {

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
	std::vector<float*> empty;
	std::vector<float*> leftData, rightData, upData, downData;
	MakeQuadBallTreeSplit(n, d, data, A, B, C, D);

	leftData.push_back(A);
	rightData.push_back(B);
	upData.push_back(C);
	downData.push_back(D);
	for (int i = 0; i < n; i++) {
		float ADist, BDist, CDist, DDist;
		if (data[i]==A)	continue;
		if (data[i]==B)	continue;
		if (data[i]==C) continue;
		if (data[i]==D) continue;

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
	float **leftData_ = parseFloatArr(leftData);
	leftData.swap(empty);
	float **rightData_ = parseFloatArr(rightData);
	rightData.swap(empty);
	float **upData_ = parseFloatArr(upData);
	upData.swap(empty);
	float **downData_ = parseFloatArr(downData);
	downData.swap(empty);



	node->isLeftLeaf = leftn <= N0;
	node->isRightLeaf = rightn <= N0;
	node->isUpLeaf = upn <= N0;
	node->isDownLeaf = downn <= N0;
	buildQuadBall(node->left, node, leftn, d, leftData_, true, true);
	buildQuadBall(node->right, node, rightn, d, rightData_, true, false);
	buildQuadBall(node->up, node, upn, d, upData_, false, true);
	buildQuadBall(node->down, node, downn, d, downData_, false, false);
}
/*
 *QuadTreeSearch
 *@param: query, node, mip
 *
 *类似于TreeSearch进行MIP求值
 *
 *由于按四个节点大小顺序安排访问顺序的实现较为复杂
 *对于四个子节点所采取的优先决定搜索方式是：
 *先判断left，right节点的大小比较，先按大小顺序访问left，right
 *再按大小顺序访问up，down
 *同时对于节点上限小于当前最大mip的情况不访问
 *
 */
void BallTree::QuadTreeSearch(float* query, ballTreeNode* node, Mip &mip) {
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
				QuadTreeSearch(query, node->left, mip);
			}
			node->release_left();
			if (mip.product < rightProduct&&node->right != NULL) {
				QuadTreeSearch(query, node->right, mip);
			}
			node->release_right();
		}
		else {
			if (mip.product < rightProduct&&node->right != NULL) {
				QuadTreeSearch(query, node->right, mip);
			}
			node->release_right();
			if (mip.product < leftProduct&&node->left != NULL) {
				QuadTreeSearch(query, node->left, mip);
			}
			node->release_left();
		}

		if (upProduct > downProduct) {
			if (mip.product < upProduct&&node->up != NULL) {
				QuadTreeSearch(query, node->up, mip);
			}
			node->release_up();
			if (mip.product < downProduct&&node->down != NULL) {
				QuadTreeSearch(query, node->down, mip);
			}
			node->release_down();
		}
		else {
			if (mip.product < downProduct&&node->down != NULL) {
				QuadTreeSearch(query, node->down, mip);
			}
			node->release_down();
			if (mip.product < upProduct&&node->up != NULL) {
				QuadTreeSearch(query, node->up, mip);
			}
			node->release_up();
		}


	}
}
/*
 *MakeQuadBallTreeSplit
 *@param: n, d, data, A, B, C, D
 *
 *寻找离原点最近的点作为A
 *再找离A点最远的两个点为B，C
 *最后找B点最远的点D
 *
 */
bool BallTree::MakeQuadBallTreeSplit(int n, int d, float **data, float *&A, float *&B, float *&C, float *&D){
		A = data[0];
		float *empty=new float[d];
		memset(empty, 0, d*sizeof(float));
		float maxDistance = 0;
		for (int i = 0; i < n; i++) {
			float curDistance = computeDistance(empty, data[i]);
			if (curDistance > maxDistance) {
				maxDistance = curDistance;
			  A = data[i];
			}
		}
		maxDistance = 0;
		B = data[n/2];
		C = data[n/2];
		for (int i = 0; i < n; i++) {
			float curDistance = computeDistance(A, data[i]);
			if (curDistance > maxDistance) {
				maxDistance = curDistance;
				C = B;
				B = data[i];
			}
		}
		maxDistance = 0;
		for (int i = 0; i < n; i++) {
			float curDistance = computeDistance(B, data[i]);
			if (curDistance > maxDistance&&data[i]!=A&&data[i]!=C) {
				maxDistance = curDistance;
				D = data[i];
			}
		}
		delete []empty;
		return true;
}
