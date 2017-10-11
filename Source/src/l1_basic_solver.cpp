#include "l1_basic_solver.hpp"
#include <cassert>
#include <algorithm>

BasicSolver::BasicSolver() {
  hashtable = Hashtable();
}

BasicSolver::~BasicSolver() {
}

void BasicSolver::loadHashTables(string name) {
  hashtable.read_from_file(name.c_str());
}

void BasicSolver::init(uint32_t readLength, uint32_t seedNum, uint32_t seedLength) {
  this->seedNum = seedNum;
  this->seedLength = seedLength;
  this->readLength = readLength;
}

uint32_t BasicSolver::solveDNA(string DNA) {
  assert(DNA.length() == readLength);

  uint32_t interval = readLength / seedNum;

  uint32_t totalFreq = 0;

  for (int i = 0; i < seedNum; i++) {
    string seed = DNA.substr(i * interval, seedLength);
    totalFreq += hashtable.get_frequency(hashtable.get_hash(seed.c_str()));
  }

  return totalFreq;
}

