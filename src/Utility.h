#ifndef __UTILITY_H
#define __UTILITY_H

#define L 256
#define PAGE_SIZE 65536

const int INDEX_INT_SIZE = 8;//改成8变成4叉树， 对应子节点pageid， slotid
const int INDEX_BOOL_SIZE = 4; //改成4变成4叉树，对应isRightLeaf， isLeftLeaf， isUpLeaf， isDownLeaf
const int DATA_INT_SIZE = 1;

const int METADATA_INDEX_OFFSET = sizeof(int) * 2;
const int METADATA_DATA_OFFSET = sizeof(int);

struct Rid {
    int pageid;
    int slotid;

    Rid() {
        pageid = slotid = -1;
    }
    Rid(int p, int s) {
        pageid = p;
        slotid = s;
    }
};

struct ballTreeNode {
    float radius;  //球的半径
    float *mean;  //球的圆心

    Rid myRid;  //  The rid of this node.
    Rid leftRid;  //  The rid of the left child.
    Rid rightRid;  //  The rid of the right child.
	//四叉树扩展
	Rid upRid;
	Rid downRid;
	//

    bool isLeftLeaf;  //左节点是否为叶子
    bool isRightLeaf;  //右节点是否为叶子
	//四叉树拓展
	bool isUpLeaf;
	bool isDownLeaf;
	//

    float **table;  //如果为叶子，table存数据
    int tableSize;  //table的大小

    ballTreeNode *left;  //左节点的索引
    ballTreeNode *right;  //右节点的索引
	//四叉树扩展
	ballTreeNode *up;
	ballTreeNode *down;
	//
    ballTreeNode *father;  //  pointer to father.
    bool isLeft;  //  whether the node is a left node.
	//四叉树扩展
	bool isUp;
	bool isDown;
	bool isRight;
	//

    ballTreeNode();
    ballTreeNode(float r, float* m, int d);
	void delete_memory() {
		if (table != NULL) {
			for (int i = 0; i < tableSize; i++) {
				delete[] table[i];
			}
			delete[] table;
			table = NULL;
		}
	}
	void release_right() {
		if (isRightLeaf&&right != NULL) {
			right->delete_memory();
			delete[]right->mean;
			delete right;
			right = NULL;
		}
	}
	void release_left() {
		if (isLeftLeaf&&left != NULL) {
			left->delete_memory();
			delete[]left->mean;
			delete left;
			left = NULL;
		}
	}
	void release_up() {
		if (isUpLeaf&&up != NULL) {
			up->delete_memory();
			delete[]up->mean;
			delete up;
			up = NULL;
		}
	}
	void release_down() {
		if (isDownLeaf&&down != NULL) {
			down->delete_memory();
			delete[]down->mean;
			delete down;
			down = NULL;
		}
	}
};

bool read_data(int n, int d, float** &data, const char* file_name);

float **parseFloatArr(std::vector<float*> v);

int parsePageId(std::streampos filePtr);
int parseSlotId(std::streampos filePtr, int slotsize, int numSlot);

//--------------------ZJQ:20170521任务3与任务4实现--------------------------//
struct Mip {
    float product;
    int index;
};
#endif
