#include "l1_basic_solver.hpp"
#include <cassert>
#include <algorithm>

BasicSolver::BasicSolver() {
  hashtable = NULL;
}

BasicSolver::~BasicSolver() {
}

void BasicSolver::loadTables(void * hashtable) {
  this->hashtable = (Hashtable*)hashtable;
}

uint32_t BasicSolver::get_seed_size() {
  return hashtable->get_seed_size();
}

void BasicSolver::init(uint32_t readLength, uint32_t seedNum, uint32_t seedLength, uint32_t limit) {
  this->seedNum = seedNum;
  this->limit = limit;
  this->seedLength = seedLength;
  this->readLength = readLength;
}

int BasicSolver::solveDNA(string DNA, uint8_t * seeds) {
  assert(DNA.length() == readLength);

  uint32_t interval = readLength / seedNum;

  int totalFreq = 0;

  for (int i = 0; i < seedNum; i++) {
    string seed = DNA.substr(i * interval, seedLength);
    totalFreq += hashtable->get_frequency(hashtable->get_hash(seed.c_str()));
    seeds[i] = i * interval;
  }

  return totalFreq;
}

