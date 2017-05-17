#ifndef __BALL_TREE_H
#define __BALL_TREE_H

#define N0 20

class BallTree {
public:
<<<<<<< HEAD
=======
	treeNode* root;

>>>>>>> 030243d7d465e50ec1594e75aa710da193b75b05
	BallTree();
	~BallTree();

	bool buildTree(
		int n,
		int d,
		float** data);

	bool storeTree(
		const char* index_path);

	bool restoreTree(
		const char* index_path);
	
	int mipSearch(
		int d,
		float* query);

	// optional
	bool insertData(
		int d,
		float* data);

	// optional
	bool deleteData(
		int d,
		float* data);

	// optional
	bool buildQuadTree(
		int n,
		int d,
		float** data);
};

<<<<<<< HEAD
=======

class treeNode {
public:
	float* mean;
	float R;
	//float *dataIndex;
	int dataSize;

	treeNode *left;
	treeNode *right;

	treeNode();
	~treeNode(){
		delete []mean;
	}
	treeNode(float* m, float r, int d) {
		mean = new float [d];
		for (int i = 0; i < d; i++) {
			mean[i] = m[i];
		}
		delete [] m;
		R = r;
		size = s;
		/*dataIndex = new int[s];
		for (int i = 0; i < sizel i++) {
			dataIndex[i] = d[i];
		}*/
	}
}

struct AB {
	int* AIndex;
	int* BIndex;
	int asize;
	int bsize;
}
>>>>>>> 030243d7d465e50ec1594e75aa710da193b75b05
#endif
