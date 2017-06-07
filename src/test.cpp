#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <fstream>
#include <cmath>
#include<time.h>

#include "Utility.h"
#include "BallTree.h"

#define MNIST

#ifdef MNIST
char dataset[L] = "datasets/Mnist";
int n = 60000, d = 51;
int qn = 1000;
#endif

#ifdef NETFLIX
char dataset[L] = "datasets/Netflix";
int n = 17770, d = 51;
int qn = 1000;
#endif

int main() {
  time_t Begin_t, End_t;
	Begin_t=time(NULL);
	char data_path[L], query_path[L];
	char index_path[L], output_path[L];
	float** data = NULL;
	float** query = NULL;

	sprintf(data_path, "%s/src/dataset.txt", dataset);
	sprintf(query_path, "%s/src/query.txt", dataset);
	sprintf(index_path, "%s/index", dataset);
	sprintf(output_path, "%s/dst/answer.txt", dataset);

	if (!read_data(n, d, data, data_path)) {
		return 1;
	}

	printf("read data finished\n");
	//system("pause");
	BallTree ball_tree1;
	ball_tree1.buildQuadTree(n, d, data);
	printf("buildTree finished\n");
	//system("pause");
	ball_tree1.storeTree(index_path);
	printf("storeTree finished\n");
	//system("pause");

	if (!read_data(qn, d, query, query_path)) {
		printf("Read query faild.\n");
	}
	printf("read query finished\n");
	//system("pause");
	FILE* fout = fopen(output_path, "w");
	if (!fout) {
		printf("can't open %s!\n", output_path);
		return 1;
	}

	BallTree ball_tree2;
	ball_tree2.restoreTree(index_path);
	printf("restoreTree finished\n");
	//system("pause");
	for (int i = 0; i < qn; i++) {
		int index = ball_tree2.mipSearch(d, query[i]);
		printf("In query %d, %d found.\n", i+1, index);
		fprintf(fout, "%d\n", index);
	}
	fclose(fout);
	printf("MIP search finished\n");
	//system("pause");

	for (int i = 0; i < n; i++) {
		delete[] data[i];
	}

	for (int i = 0; i < qn; i++) {
		delete[] query[i];
	}
	End_t=time(NULL);
	printf("Program takes %d s.\n", End_t-Begin_t);

	//system("pause");
	return 0;
}
