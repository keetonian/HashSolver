#include <vector>
#include "l2hashtable.hpp"
#include "hashtable.hpp"
#include <string>

using namespace std;

class HobbesSolver {
public:
	HobbesSolver();
	~HobbesSolver();
	//void loadTree(string treeFileName);
	//void generateTree(string treeFileName);
	void init(int readLength, int seedNum, int seedLength);
	void reset();
	unsigned int solveDNA(string DNA);
	void print();
	void loadHashTables(string name);

private:
	//HashTree tree;
	int *invertedList;
	int **dynamicMatrix;
	int readLength;
	int seedNum;
	int seedLength;
	int *defaultInvertedList;
	int *defaultDynamicMatrix;
	int dynamicMatrixWidth;
	Hashtable hashtable;
	L2Hashtable l2hashtable;
};
