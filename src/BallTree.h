#ifndef __BALL_TREE_H
#define __BALL_TREE_H

#define N0 20

#include "Utility.h"

class BallTree {
private:
    int INDEX_SLOTSIZE;
    int DATA_SLOTSIZE;

    int curPageID;   //记录当前节点应该存在哪页
    int curSlotID;   //记录单签节点应该存在哪槽

    ballTreeNode *root;
    int dimension;
    int numIndexSlot;
    int numDataSlot;
    int id;

    float *computeMean(int n, int d, float **data);
    float computeRadius(int n, int d, float **data, float *mean);
    float computeDistance(float *x, float *y);
    bool MakeBallTreeSplit(int n, int d, float **data, float *&A, float *&B);
    void buildBall(ballTreeNode *&node, int n, int d, float **data, int fatherId, bool isLeft);

    void preorderStore(ballTreeNode *node, std::ofstream &indexOutput, 
        std::ofstream &dataOutput, std::streampos &indexPtr, std::streampos &dataPtr);
    void storeIndexNode(ballTreeNode *node, std::ofstream &output, std::streampos &filePtr);
    void storeDataNode(ballTreeNode *node, std::ofstream &output, std::streampos &filePtr);
public:
	BallTree();
	~BallTree();

	bool buildTree(int n, int d, float** data);

	bool storeTree(const char* index_path);

//-----------------------ZJQ:20170521任务3与4实现-----------------------//
	bool restoreTree(const char* index_path);
	
	int mipSearch(int d, float* query);

    void LinearSearch(float query, ballTreeNode *node, Mip& mip);

    float computeInnerProduct(int d, float* query, float* data);

    float MIP(int d, float* query, ballTreeNode* node);

    ballTreeNode* getNode(int pageID, int slotID);
//----------------------------------------------------------------------//
	// optional
	//bool insertData(int d, float* data);

	// optional
	//bool deleteData(int d, float* data);

	// optional
	//bool buildQuadTree(int n, int d, float** data);

    int getNumIndexSlot() { return numIndexSlot; }
    int getNumDataSlot() { return numDataSlot; }
    int getIndexSlotSize() { return INDEX_SLOTSIZE; }
    int getDataSlotSize() { return DATA_SLOTSIZE; }
};
#endif
