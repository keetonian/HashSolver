#include <cstring>
#include <cassert>
#include <iostream>
#include <climits>
#include <limits>
#include "l1_hobbes_solver.hpp"
#include "l1_hashtable.hpp"

HobbesSolver::HobbesSolver() {
  invertedList = NULL;
  defaultInvertedList = NULL;
  dynamicMatrix = NULL;
  defaultDynamicMatrix = NULL;
  readLength = 0;
  seedNum = 0;
  seedLength = 0;
  limit = 0;
  hashtable = NULL;
}

HobbesSolver::~HobbesSolver() {
  if  (invertedList != NULL) {
    delete [] invertedList;
    invertedList = NULL;
  }
  if  (defaultInvertedList != NULL) {
    delete [] defaultInvertedList;
    defaultInvertedList = NULL;
  }
  if  (dynamicMatrix != NULL) {
    delete [] dynamicMatrix;
    dynamicMatrix = NULL;
  }
  if  (defaultDynamicMatrix != NULL) {
    delete [] defaultDynamicMatrix;
    defaultDynamicMatrix = NULL;
  }
}

void HobbesSolver::init(uint32_t readLength, uint32_t seedNum, uint32_t seedLength, uint32_t limit) {
  //Initialize the matrix
  if  (invertedList != NULL) {
    delete [] invertedList;
    invertedList = NULL;
  }
  if  (defaultInvertedList != NULL) {
    delete [] defaultInvertedList;
    defaultInvertedList = NULL;
  }
  if  (dynamicMatrix != NULL) {
    delete [] dynamicMatrix;
    dynamicMatrix = NULL;
  }
  if  (defaultDynamicMatrix != NULL) {
    delete [] defaultDynamicMatrix;
    defaultDynamicMatrix = NULL;
  }

  this->limit = limit;
  this->readLength = readLength;
  this->seedNum = seedNum;
  this->seedLength = seedLength;
  dynamicMatrixWidth = readLength + 1 - seedNum * seedLength - limit;

  invertedList = new int [readLength - seedLength + 1];
  defaultInvertedList = new int [readLength - seedLength + 1];

  dynamicMatrix = new int* [seedNum + 1];
  dynamicMatrix[0] = new int [dynamicMatrixWidth * (seedNum + 1)];
  defaultDynamicMatrix = new int [dynamicMatrixWidth * (seedNum + 1)];

  for (unsigned int i = 1; i <= seedNum; i++)
    dynamicMatrix[i] = dynamicMatrix[i-1] + dynamicMatrixWidth;

  for (unsigned int i = 0; i < readLength - seedLength + 1; i++)
    defaultInvertedList[i] = 0;

  for (unsigned int i = 0; i < dynamicMatrixWidth * (seedNum + 1); i++)
    defaultDynamicMatrix[i] = 0;

}

void HobbesSolver::reset() {
  memcpy(dynamicMatrix[0], defaultDynamicMatrix, dynamicMatrixWidth * ( this->seedNum + 1));
  memcpy(invertedList, defaultInvertedList, readLength - this->seedLength + 1);
}

void HobbesSolver::loadTables(void * hashtable) {
  this->hashtable = (Hashtable *)hashtable;
}

uint32_t HobbesSolver::get_seed_size() {
  return hashtable->get_seed_size();
}

int HobbesSolver::solveDNA(string DNA, uint8_t * seeds) {
  assert(DNA.length() == readLength);
  reset();

  uint64_t hash;
  for (unsigned int i = 0; i < readLength - seedLength + 1; i++) {
    if (!i) {
      string seed = DNA.substr(i, seedLength);
      hash = hashtable->get_hash(seed.c_str());
    }
    else {
      hash = hashtable->get_hash(DNA[i+seedLength-1], hash);
    }
    //invertedList[i] = hashtable->get_frequency(hash);
  }

  for (int counter = 1; counter < seedNum + 1; counter++) {
    dynamicMatrix[counter][0] = invertedList[(counter-1) * this->seedLength] + dynamicMatrix[counter - 1][0];
  }
  for (int heightCount = 1; heightCount < seedNum + 1; heightCount++) {
    for (int widthCount = 1; widthCount < dynamicMatrixWidth; widthCount++) {
      int indexNum = (heightCount - 1) * seedLength + widthCount;
      dynamicMatrix[heightCount][widthCount] = dynamicMatrix[heightCount-1][widthCount] + invertedList[indexNum];
      if (dynamicMatrix[heightCount][widthCount-1] < dynamicMatrix[heightCount][widthCount]) 
	dynamicMatrix[heightCount][widthCount] = dynamicMatrix[heightCount][widthCount-1];

    }   
  }

  int width = dynamicMatrixWidth - 1;
  int invertedListIndex = readLength - seedLength - limit;
  int prevAcc, seedFreq;

  for(int heightCount = seedNum; heightCount >= 0; heightCount --) 
  {
    while(dynamicMatrix[heightCount][width] == dynamicMatrix[heightCount][width - 1] && width > 0)
    {
      --width;
    }

    if(heightCount < seedNum) {
      seedFreq = prevAcc - dynamicMatrix[heightCount][width];
      while(seedFreq != invertedList[invertedListIndex]) {
	--invertedListIndex;
      }
      seeds[heightCount] = invertedListIndex;
      invertedListIndex -= seedLength;
    }
    prevAcc = dynamicMatrix[heightCount][width];
  }

  return dynamicMatrix[seedNum][dynamicMatrixWidth-1];
}
