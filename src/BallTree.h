#ifndef __BALL_TREE_H
#define __BALL_TREE_H

#define N0 20

#include "Utility.h"

class BallTree {
private:
    int INDEX_SLOTSIZE;
    int DATA_SLOTSIZE;

    Rid curIndexRid;    //  Current rid of the non-leaf node.
    Rid curDataRid;    //  Current rid of the leaf node.

    ballTreeNode *root;
    int dimension;
    int numIndexSlot;
    int numDataSlot;

    char *indexEntry_path;
    char *dataEntry_path;

    void computeMean(float *&mean, int n, int d, float **data);
    float computeRadius(int n, int d, float **data, float *mean);
    float computeDistance(float *x, float *y);
    bool MakeBallTreeSplit(int n, int d, float **data, float *&A, float *&B);
    void buildBall(ballTreeNode *&node, ballTreeNode *father, int n, int d, float **data, bool isLeft);

    void preorderStore(ballTreeNode *node, ballTreeNode *father, std::ofstream &indexOutput, 
        std::ofstream &dataOutput, bool isLeft);
    void bfsStore(std::ofstream &indexOutput, std::ofstream &dataOutput);
    void storeIndexNode(ballTreeNode *node, std::ofstream &output, Rid &rid);
    void storeDataNode(ballTreeNode *node, std::ofstream &output, Rid &rid);
    void updateRid(ballTreeNode *node, std::ofstream &output);
    
    void TreeSearch(float *query, ballTreeNode *node, Mip &mip);
public:
	BallTree();
	~BallTree();

	bool buildTree(int n, int d, float** data);

	bool storeTree(const char* index_path);

//-----------------------ZJQ:20170521任务3与4实现-----------------------//
	bool restoreTree(const char* index_path);
	
	int mipSearch(int d, float* query);

    void LinearSearch(float *query, ballTreeNode *node, Mip& mip);

    float computeInnerProduct(float* query, float* data);

    float MIP(float* query, ballTreeNode* node);

    ballTreeNode* getNode(int pageID, int slotID, bool isIndex);
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
