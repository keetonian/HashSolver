#include "fastHASHSolver.h"
#include <numeric>
#include <cassert>
#include <algorithm>
#include <cstring>

FastHASHSolver::FastHASHSolver() {
}

FastHASHSolver::~FastHASHSolver() {
}

void FastHASHSolver::init(unsigned int readLength, unsigned int seedNum, unsigned int seedLength, uint32_t limit) {
  this->seedNum = seedNum;
  this->seedLength = seedLength;
  this->readLength = readLength;

  seedFreq.resize(readLength / seedLength);
  offset.resize(readLength / seedLength);
  idx.resize(readLength / seedLength);
  iota(idx.begin(), idx.end(), 0);
}

bool binaryRangeSearch(uint32_t min, uint32_t max, uint32_t offset, uint32_t frequency, Hashtable * hashtable) {
  uint64_t low = offset, high = offset + frequency - 1;
  uint64_t mid, location;
  while(low <= high) {
    mid = (low+high)/2;
    location = hashtable->get_location(mid);
    if(location >= min && location <= max) 
      return true;
    else if(hashtable->get_location(low) > max)   return false;
    else if(hashtable->get_location(high) < min)  return false;
    else if(location > max)
      if(mid>0) 	high = mid - 1;
      else		return false;
    else
      low = mid + 1;
  }
  return false;
}

void FastHASHSolver::loadTables(void * hashtable) {
  this->hashtable = (Hashtable*)hashtable;
}

uint32_t FastHASHSolver::get_seed_size() {
  return hashtable->get_seed_size();
}

int FastHASHSolver::solveDNA(const string &DNA, std::vector<uint32_t> &locations) {
  assert(DNA.length() == readLength);

  uint32_t location;
  uint32_t recur_location;
  uint32_t min, max;
  unsigned int totalFreq = 0;
  int error = seedNum - 1;
  int matches = 0;

  //vector<unsigned int> seedFreq;	
  //seedFreq.resize(readLength/seedLength);

  for (int i = 0; i < readLength / seedLength; i++) {
    string seed = DNA.substr(i * seedLength, seedLength);
    uint64_t hash = hashtable->get_hash(seed.c_str()); 

    seedFreq[i] = hashtable->get_frequency(hash);
    offset[i] = hashtable->get_offset(hash);
  }

  sort(idx.begin(), idx.end(), [&](const size_t i1, const size_t i2) {return seedFreq[i1] < seedFreq[i2];});


  for(int i = 0; i < readLength / seedLength; i++) {
    if(idx[i] >= seedNum) continue;

    for(int j=0; j < seedFreq[idx[i]]; j++) {
      location = hashtable->get_location(offset[idx[i]] + j);
      matches = 1;		

      for(int k=0; k < readLength / seedLength; k++) {
	if(k==i) continue;
	min = location + (idx[k] - idx[i])*seedLength - error;
	max = location + (idx[k] - idx[i])*seedLength + error;
	if(binaryRangeSearch(min, max, offset[idx[k]], seedFreq[idx[k]], hashtable))
	  matches++;
	/*
	   for(int l=0; l < seedFreq[idx[k]]; l++) {
	   recur_location = hashtable->get_location(offset[idx[k]] + l);
	   if(recur_location >= (location + (idx[k] - idx[i])*seedLength - error) && recur_location <= (location + (idx[k] - idx[i])*seedLength + error)) {
	   matches++;
	   break;
	   }
	   }
	   */
	if(matches >= readLength/seedLength - error) {
	  totalFreq++;
	  locations.push_back(location - seedLength * idx[i]);
	  break;
	}
	if(readLength/seedLength - error - matches > readLength/seedLength - k + 1)
	  break;
      }	

    }

  }

  return totalFreq;
}

