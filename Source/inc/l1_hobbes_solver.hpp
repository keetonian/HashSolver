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
	void init(uint32_t readLength, uint32_t seedNum, uint32_t seedLength, uint32_t limit);
	void reset();
	int solveDNA(string DNA, uint8_t * seeds);
	int solveDNA(std::string DNA, std::vector<uint32_t>&) {return 0;}
	void print();
	void loadTables(void * hashtable);
	uint32_t get_seed_size();

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
	uint32_t limit;
	Hashtable * hashtable;
};
