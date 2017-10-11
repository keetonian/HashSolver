#include <vector>
#include <string>
#include "l1_hashtable.hpp"
#include "seed_solver.hpp"
#include <stdint.h>

using namespace std;

class BasicSolver : public SeedSolver {
public:
	BasicSolver();
	~BasicSolver();

	//void loadTree(string treeFileName);
	void init(uint32_t readLength, uint32_t seedNum, uint32_t seedLength);
	uint32_t solveDNA(string DNA);
	void loadHashTables(string name);

private:
	//HashTree tree;
	uint32_t readLength;
	uint32_t seedNum;
	uint32_t seedLength;
	Hashtable hashtable;

};

