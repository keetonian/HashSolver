#ifndef SEED_SOLVER_H_
#define SEED_SOLVER_H_

#include <stdint.h>
#include "l1_hashtable.hpp"

class SeedSolver {
  public:
    virtual int solveDNA(std::string DNA, uint8_t * seeds) = 0;
    virtual void init(uint32_t readLength, uint32_t seedNum, uint32_t seedLength, uint32_t limit) = 0;
    virtual void loadHashTables(Hashtable * hashtable) = 0;
    virtual uint32_t get_seed_size() = 0;

};
#endif //SEED_SOLVER_H_
