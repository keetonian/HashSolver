#include "l1_hashtable.hpp"
#include <vector>
#include <string>
#include "seed_solver.hpp"
#include <stdint.h>

using namespace std;

class FastHASHSolver : public SeedSolver {
public:
	FastHASHSolver();
	~FastHASHSolver();

	//void loadTree(string treeFileName);
	void init(unsigned int readLength, unsigned int seedNum, unsigned int seedLength, uint32_t lmit);
	int solveDNA(const string &DNA, uint8_t * seeds) {return 0x0;}
	int solveDNA(const string &DNA, std::vector<uint32_t>& locations);
	int solveDNA(const string &DNA, uint8_t *, std::vector<uint32_t>& locations) {return 0;};
	void loadTables(void * hashtable);
	uint32_t get_seed_size();

private:
	vector<unsigned int> seedFreq;
	vector<unsigned int> offset;
	vector<unsigned int> idx;
	//HashTree tree;
	unsigned int readLength;
	unsigned int seedNum;
	unsigned int seedLength;
	Hashtable * hashtable;
};

