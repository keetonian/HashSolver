#include <vector>
#include "l1_hashtable.hpp"
#include "seed_solver.hpp"
#include <string>
#include <stdint.h>

using namespace std;

class HobbesSolver : public SeedSolver {
public:
	HobbesSolver();
	~HobbesSolver();
	//void loadTree(string treeFileName);
	//void generateTree(string treeFileName);
	void init(uint32_t readLength, uint32_t seedNum, uint32_t seedLength);
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
};
