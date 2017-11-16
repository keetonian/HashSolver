#include "fast_hash_solver.hpp"
#include <cassert>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
 
vector<uint32_t> SeedFrequencies;
vector<uint32_t> SeedOffsets;

FastHashSolver::FastHashSolver() {
}

FastHashSolver::~FastHashSolver() {
}

void FastHashSolver::init(unsigned int readLength, unsigned int seedNum, unsigned int seedLength, unsigned int totalSeeds, unsigned int algo) {
  this->seedNum = seedNum;
  this->seedLength = seedLength;
  this->readLength = readLength;
  this->totalSeeds = totalSeeds;
  this->algo = algo;

  SeedFrequencies.resize(totalSeeds);
  SeedOffsets.resize(totalSeeds);
  idx.resize(totalSeeds);
  iota(idx.begin(), idx.end(), 0);
}

bool FastHashSolver::binaryRangeSearch(uint32_t min, uint32_t max, uint32_t offset, uint32_t frequency) {
  uint64_t low = offset, high = offset + frequency - 1;
  uint64_t mid, location;
  while(low <= high) {
    mid = (low+high)/2;
    if(algo == 0){
      location = hashtable->get_location(mid);
      if(location >= min && location <= max) 
	return true;
      else if(hashtable->get_location(low) > max) 
	return false;
      else if(hashtable->get_location(high) < min) 
	return false;
      else if(location > max)
	if(mid>0) 	
	  high = mid - 1;
	else
	  return false;
      else
	low = mid + 1;
    }
    else {
      location = bwt_sa(bwt, mid);
      if(location >= min && location <= max) 
	return true;
      else if(bwt_sa(bwt, low) > max) 
	return false;
      else if(bwt_sa(bwt, high) < min) 
	return false;
      else if(location > max)
	if(mid>0) 	
	  high = mid - 1;
	else
	  return false;
      else
	low = mid + 1;
    }
  }
  return false;
}

unsigned int FastHashSolver::solveDNA(uint8_t* seeds, vector<uint32_t> & locations) {

  uint32_t location;
  uint32_t recur_location;
  uint32_t min, max;
  unsigned int totalFreq = 0;
  int error = seedNum - 1;
  int matches = 0;

  sort(idx.begin(), idx.end(), [&](const size_t i1, const size_t i2) {return SeedFrequencies[i1] < SeedFrequencies[i2];});

  for(int i = 0; i < totalSeeds; i++) {
    if(idx[i] >= seedNum) continue;

    for(int j=0; j < SeedFrequencies[idx[i]]; j++) {
      location = hashtable->get_location(SeedOffsets[idx[i]] + j);
      matches = 1;	
      int count = 0;
      for(int k=0; k < readLength / seedLength; k++) {
	count++;
	if(k==i) continue;
	min = location + seeds[idx[k]] - seeds[idx[i]] - error;
	max = location + seeds[idx[k]] - seeds[idx[i]] + error;
	if(binaryRangeSearch(min, max, SeedOffsets[idx[k]], SeedFrequencies[idx[k]]))
	  matches++;
	if(matches >= totalSeeds - error) {
	  locations.push_back(location-seedLength * idx[i]);
	  // Location is found here
	  totalFreq++;
	  break;
	}
	if(totalSeeds - error - matches > totalSeeds - k + 1)
	  break;
      }	
    }

  }

  return totalFreq;
}

