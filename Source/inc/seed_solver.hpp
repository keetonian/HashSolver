#ifndef SEED_SOLVER_H_
#define SEED_SOLVER_H_

#include <set>
#include <stdint.h>

class SeedSolver {
  public:
    virtual uint32_t solveDNA(std::string DNA) = 0;
    virtual void init(uint32_t readLength, uint32_t seedNum, uint32_t seedLength) = 0;
    virtual void loadHashTables(std::string name) = 0;

};
#endif //SEED_SOLVER_H_
