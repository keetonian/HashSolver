#ifndef FASTHASHSOLVER_H_
#define FASTHASHSOLVER_H_

#include "l1_hashtable.hpp"
#include <vector>
#include <string>
#include "bwt.h"

using namespace std;

extern vector<uint32_t> SeedFrequencies;
extern vector<uint32_t> SeedOffsets;

class FastHashSolver {
  public:
    FastHashSolver();
    ~FastHashSolver();

    //void loadTree(string treeFileName);
    void init(unsigned int readLength, unsigned int seedNum, unsigned int seedLength, unsigned int totalSeeds, unsigned int algo);
    unsigned int solveDNA(uint8_t* seeds, vector<uint32_t> & locations);
    void load_hashtable(Hashtable * table) {hashtable = table;}
    void load_bwt(bwt_t * b) {bwt = b;}
    bool binaryRangeSearch(uint32_t, uint32_t, uint32_t, uint32_t);

  private:
    vector<unsigned int> idx;
    //HashTree tree;
    unsigned int readLength;
    unsigned int seedNum;
    unsigned int seedLength;
    unsigned int totalSeeds;
    Hashtable * hashtable;
    unsigned int algo;
    bwt_t * bwt;

};


#endif //FASTHASHSOLVER_H_
