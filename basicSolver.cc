#include "basicSolver.h"
#include <cassert>
#include <algorithm>

BasicSolver::BasicSolver() {
}

BasicSolver::~BasicSolver() {
}

void BasicSolver::init(unsigned int readLength, unsigned int seedNum, unsigned int seedLength) {
  this->seedNum = seedNum;
  this->seedLength = seedLength;
  this->readLength = readLength;
}

unsigned int BasicSolver::solveDNA(string DNA) {
  assert(DNA.length() == readLength);

  unsigned int interval = readLength / seedNum;

  unsigned int totalFreq = 0;

  for (int i = 0; i < seedNum; i++) {
    string seed = DNA.substr(i * interval, seedLength);
    totalFreq += get_frequency(get_hash(seed.c_str()));
  }

  return totalFreq;
}

