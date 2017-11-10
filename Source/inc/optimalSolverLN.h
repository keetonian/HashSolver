#include <string>
#include <ostream>
#include <vector>
#include <algorithm>
#include "bwt.h"
#include "seed_solver.hpp"

using namespace std;

struct Cell {
	Cell() {
		start = -1;
		end = -1;
		isleaf = false;
		frequency = -1;
		sa_begin = 0;
		sa_end = 0;
		lstart = -1;
		lend = -1;
		lfreq = -1;
		rstart = -1;
		rend = -1;
		rfreq = -1;
	};

	unsigned int start;
	unsigned int end;
	bool isleaf;
	unsigned int frequency;
	bwtint_t sa_begin;
	bwtint_t sa_end;

	unsigned int lstart;
	unsigned int lend;
	unsigned int lfreq;
	unsigned int rstart;
	unsigned int rend;
	unsigned int rfreq;

};

class OptimalSolverLN : public SeedSolver {
public:
	OptimalSolverLN();
	~OptimalSolverLN();
	void loadBwt(bwt_t *bwt, int min_length);
	void init(uint32_t readLength, uint32_t seedNum, uint32_t, uint32_t);
	void reset();
	int solveDNA(string DNA, std::vector<uint32_t> &);
	int solveDNA(string DNA, uint8_t * seeds) {return 0;}
	uint32_t get_seed_size() {return 0;}
	void loadTables(void * table);
	void fillMatrix(string DNA);
	unsigned int calculateLastDiv();
	unsigned int calcualteFreq();
	void printSeeds(ostream& stream);
	void printFreqs(ostream& stream);
	void printLength(ostream& stream);
	void printStats(ostream& stream);
	void sortOfFreq();
	void sortOfLength();
	
	void backtrack();
	
	//For debugging
	void setMinLength(int minLength);
	void feedL0();

private:
	const int MIN_LENGTH = 10;
	//Internal functions
	//Load the first level (level[0])
	void loadL0(string DNA);
	//Returns the location of the first optimal div in read
	int solveFirstOptimal(int opt_div, int pos, int l);
	//Sort the seeds
	template<class T>
	void sortSeeds(T relation);
	static bool compFreq(Cell left, Cell right);
	static bool compLength(Cell left, Cell right);

	//This is for debugging
	bool L0Loaded;
	
	bwt_t *bwt;

	//Internal data structures
	int readLength;
	int seedNum;
	int minLength;

	Cell* matrix;
	Cell* base;
	vector<Cell*> level;
	vector<Cell*> level0;
	Cell* defaultMatrix;
	Cell* defaultBase;
	int matrixSize;
	int baseSize;

	vector<Cell> seeds;

	unsigned int finalDiv;

	//Statistics
	vector<unsigned long long> divTravel;
	unsigned int processedReads;
};
