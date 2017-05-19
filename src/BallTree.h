#ifndef __BALL_TREE_H
#define __BALL_TREE_H

#define N0 20

#include "Utility.h"

class BallTree {
private:
    int DATA_SLOTSIZE;

    ballTreeNode *root;
    int dimension;
    int numSlot;
    int id;

    float *computeMean(int n, int d, float **data);
    float computeRadius(int n, int d, float **data, float *mean);
    float computeDistance(float *x, float *y);
    bool MakeBallTreeSplit(int n, int d, float **data, float *&A, float *&B);
    void buildBall(ballTreeNode *&node, int n, int d, float **data, int fatherId);

    void preorderStore(ballTreeNode *node, std::ofstream &indexOutput, 
        std::ofstream &dataOutput, std::streampos &filePtr, std::streampos &filePtr1);
public:
	BallTree();
	~BallTree();

	bool buildTree(int n, int d, float** data);

	bool storeTree(const char* index_path);

	//bool restoreTree(const char* index_path);
	
	//int mipSearch(int d, float* query);

	// optional
	//bool insertData(int d, float* data);

	// optional
	//bool deleteData(int d, float* data);

	// optional
	//bool buildQuadTree(int n, int d, float** data);

    int getNumSlot() { return numSlot; }
};
#endif
