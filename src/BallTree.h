#ifndef __BALL_TREE_H
#define __BALL_TREE_H

#define N0 20

#include "Utility.h"

class BallTree {
private:
    ballTreeNode *root;
    int dimension;
    int numSlot;

    float *computeMean(int n, int d, float **data);
    float computeRadius(int n, int d, float **data, float * mean);
    float computeDistance(float *x, float *y);
    bool MakeBallTreeSplit(int n, int d, float **data, float *&A, float *&B);
    void buildBall(ballTreeNode *&node, int n, int d, float **data);

    void preorderStore(ballTreeNode *node, std::ofstream &output, int pageid, int slotid);
public:
	BallTree();
	~BallTree();

	bool buildTree(int n, int d, float** data);

	bool storeTree(const char* index_path);

	bool restoreTree(const char* index_path);
	
	int mipSearch(int d, float* query);

	// optional
	bool insertData(int d, float* data);

	// optional
	bool deleteData(int d, float* data);

	// optional
	bool buildQuadTree(int n, int d, float** data);
};
#endif
