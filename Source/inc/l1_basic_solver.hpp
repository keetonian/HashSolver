#ifndef BASIC_SOLVER_H_
#define BASIC_SOLVER_H_

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
	void init(uint32_t readLength, uint32_t seedNum, uint32_t seedLength, uint32_t limit);
	int solveDNA(const string &DNA, uint8_t * seeds);
	int solveDNA(const std::string &DNA, std::vector<uint32_t>&) {return 0;}
	int solveDNA(const std::string &DNA, uint8_t *, std::vector<uint32_t>&) {return 0;}
	void loadTables(void * hashtable);
	uint32_t get_seed_size();

private:
	//HashTree tree;
	uint32_t readLength;
	uint32_t seedNum;
	uint32_t seedLength;
	uint32_t limit;
	Hashtable * hashtable;

};


#endif // BASIC_SOLVER_H_
