# BallTree for Maximum Inner-Product Search
## The idea

The idea is raised by a paper **Maximum Inner-Product Search Using Cone Trees**, which is given by Parikshit and Alexander. The idea is that we use a ball tree to partition the dataset until a particular number of tuples(N0) is reached. Thus the leaf node must contain equal or less tuples than the number. For each ball, we can easily compute the upper bound of the inner-product between the ball and the query by the mean and the radius. Therefore, we can traverse from the root to a specific leaf and finally find the inner-product under N0 times computations between a tuple and a query.

## Implementation

This implementation is the outer-storage version. The ball trees are build in the memory and are written into binary files. When using the trees to find the maximum inner-product, just load the trees from the files page by page and then do the process work. The implementation is divided into several tasks.

## Start

Use `make` in the root directory to build the excution file.

Then execute `bin/main` to launch the application in the root directory. 

Note that you should always stay in the root directory, otherwise the relative paths won't work.

## Tasks

1. Tree Building. --- Finished.
2. Tree Storage. --- Finished.
3. Tree Loading. --- Finished.
4. Maximum Inner-product and Tree Trimming. --- Finished.
5. (optional) Tree insertion and deletion. --- 0%.
6. (optional) Quadratic Ball Tree Implementation. --- 90%.

### Tree Storage

1. Store the nodes through pages and slots with each page size of 64KB.
2. Every slot only stores one tree node.
3. Every node has an ID, and access as well as storage are implemented by ID.
4. (Reference) Sperate the storage of leaf nodes and non-leaf nodes.
5. (Reference) Tuples are straightly stored into non-leaf nodes.

### Tree Loading

Finish the queries with as little memory as possible and as fast as possible.

## Datasets

Mnist (For debugging) and Netflix

Datasets | \#Objects | \#Dimension | \#Query | B | Data Size
-------- | -------- | ---------- | ------ | - | ---------
Mnist | 60,000 | 50 | 1000 | 64KB | 9.4MB
Netflix | 17770 | 50 | 1000 | 64KB | 8.5MB

## Storage Format

Non-Leaf Node

{ int\[4\](leftRid and rightRid), float(radius), float\[dimension\](mean), bool\[2\](isLeftChild, isRightChild) }

Leaf Node

{ int(k), float\[dimension\](mean), float(radius), float\[N0\]\[dimension\](k tuples) }
