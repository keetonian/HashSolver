#include "l1_basic_solver.hpp"
#include <cassert>
#include <algorithm>

extern vector<uint32_t> SeedFrequencies;
extern vector<uint32_t> SeedOffsets;

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

int BasicSolver::solveDNA(const string &DNA, uint8_t * seeds) {
  assert(DNA.length() == readLength);

  uint32_t interval = readLength / seedNum;

  int totalFreq = 0;


  for (uint32_t i = 0; i < seedNum; i++) {
    uint32_t hash = hashtable->get_hash(&(DNA[i*interval]));
    uint32_t frequency = hashtable->get_frequency(hash);
    totalFreq += frequency;
    SeedFrequencies[i] = frequency; // FAST HASH
    SeedOffsets[i] = hashtable->get_offset(hash); // FAST HASH
    seeds[i] = i * interval;
  }

  return totalFreq;
}

