#include <vector>
#include <string>
#include "hashtable.hpp"

using namespace std;

class BasicSolver {
public:
	BasicSolver();
	~BasicSolver();

	//void loadTree(string treeFileName);
	void init(unsigned int readLength, unsigned int seedNum, unsigned int seedLength);
	unsigned int solveDNA(string DNA);
	void loadHashTables(string name);

private:
	//HashTree tree;
	unsigned int readLength;
	unsigned int seedNum;
	unsigned int seedLength;
	Hashtable hashtable;

};

