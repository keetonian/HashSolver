#include "l2_hashtable.hpp"
#include "l1_hashtable.hpp"
#include "stdio.h"
#include <iostream>

L2Hashtable::L2Hashtable() {
}

L2Hashtable::L2Hashtable(uint32_t threshold){
  this->l2threshold = threshold;
}

L2Hashtable::L2Hashtable(uint64_t num_tables, uint64_t locations, uint64_t overflow){
  
}

L2Hashtable::L2Hashtable(std::string filename){
  this->l2_read_from_file(filename.c_str());
}

L2Hashtable::~L2Hashtable() {
  if (l2hashtable)
    free(l2hashtable);
  if (l2locations)
    free(l2locations);
  if (l2overflow_values)
    free(l2overflow_values);
}

std::string L2Hashtable::get_name() {
  return "Baseline";
}

uint8_t L2Hashtable::pearsonHash(const char * seed) {
  static const uint8_t T[256] = {
    // 0-255 shuffled in any (random) order suffices
    98,  6, 85,150, 36, 23,112,164,135,207,169,  5, 26, 64,165,219, //  1
    61, 20, 68, 89,130, 63, 52,102, 24,229,132,245, 80,216,195,115, //  2
    90,168,156,203,177,120,  2,190,188,  7,100,185,174,243,162, 10, //  3
    237, 18,253,225,  8,208,172,244,255,126,101, 79,145,235,228,121, //  4
    123,251, 67,250,161,  0,107, 97,241,111,181, 82,249, 33, 69, 55, //  5
    59,153, 29,  9,213,167, 84, 93, 30, 46, 94, 75,151,114, 73,222, //  6
    197, 96,210, 45, 16,227,248,202, 51,152,252,125, 81,206,215,186, //  7
    39,158,178,187,131,136,  1, 49, 50, 17,141, 91, 47,129, 60, 99, //  8
    154, 35, 86,171,105, 34, 38,200,147, 58, 77,118,173,246, 76,254, //  9
    133,232,196,144,198,124, 53,  4,108, 74,223,234,134,230,157,139, // 10
    189,205,199,128,176, 19,211,236,127,192,231, 70,233, 88,146, 44, // 11
    183,201, 22, 83, 13,214,116,109,159, 32, 95,226,140,220, 57, 12, // 12
    221, 31,209,182,143, 92,149,184,148, 62,113, 65, 37, 27,106,166, // 13
    3, 14,204, 72, 21, 41, 56, 66, 28,193, 40,217, 25, 54,179,117, // 14
    238, 87,240,155,180,170,242,212,191,163, 78,218,137,194,175,110, // 15
    43,119,224, 71,122,142, 42,160,104, 48,247,103, 15, 11,138,239  // 16
  };

  size_t i;
  //size_t j;
  uint8_t h=0, index;

  for(i=0; i<l2seed_size; i++) {
    index = h^seed[i];
    h = T[index];
  }

  return h;
}

uint32_t L2Hashtable::l2_get_threshold() {
  return l2threshold;
}

uint32_t L2Hashtable::l2_get_l2_seed_size() {
  return l2seed_size;
}

void L2Hashtable::l2_set_num_tables(uint64_t number){
  l2num_tables = number;
}

uint64_t L2Hashtable::l2_get_num_tables(){
  return l2num_tables;
}

uint64_t L2Hashtable::l2_get_table_size() {
  return l2table_size;
}

size_t L2Hashtable::l2_get_locations_size(){
  return l2locations_size;
}

void L2Hashtable::l2_init_hashtable(uint64_t num_tables){
  this->l2num_tables = num_tables;
  this->l2hashtable = (uint32_t*)malloc((1 + num_tables * l2table_size) * sizeof(uint32_t));
}

void L2Hashtable::l2_init_locations(uint64_t size){
  this->l2locations_size = size;
  this->l2locations = (uint32_t*)malloc(l2locations_size * sizeof(uint32_t));
}

void L2Hashtable::l2_init_overflow(uint64_t size){
  l2overflow_size = size;
  if(!size)
    return;
  l2overflow_values = (uint32_t*)malloc(size * sizeof(uint32_t));
}

uint64_t L2Hashtable::l2_get_index(uint32_t offset, uint32_t I, uint32_t J, uint8_t hash){
  return ((I*6)<<8) + ((J)<<8) + ((offset * 36)<<8) + hash;
}

uint64_t L2Hashtable::l2_get_overflow(uint64_t index){
  uint64_t overflow = 0;
  for(uint32_t i = 0; i < l2overflow_size; i++){
    if(index >= l2overflow_values[i])
      overflow += L2_OVERFLOW;
  }
  return overflow;
}

void L2Hashtable::l2_set_overflow(uint32_t i, uint64_t index){
  l2overflow_values[i] = index;
}

uint32_t L2Hashtable::l2_get_frequency(uint32_t offset, uint32_t I, uint32_t J, uint8_t hash){
  return l2hashtable[l2_get_index(offset, I, J, hash+1)] - l2hashtable[l2_get_index(offset, I, J, hash)];
}

uint64_t L2Hashtable::l2_get_offset(uint32_t offset, uint32_t I, uint32_t J, uint8_t hash){
  uint64_t index = l2_get_index(offset, I, J, hash);
  return l2hashtable[index] + l2_get_overflow(index);
}

uint32_t L2Hashtable::l2_get_frequency2(uint64_t index){
  return l2hashtable[index+1] - l2hashtable[index];
}

uint64_t L2Hashtable::l2_get_offset2(uint64_t index){
  return l2hashtable[index] + l2_get_overflow(index);
}

void L2Hashtable::l2_set_frequency(uint64_t index, uint32_t frequency){
  l2hashtable[index] = frequency;
}

void L2Hashtable::l2_set_offset(uint64_t index, uint32_t offset){
  l2hashtable[index] = offset;
}

uint32_t L2Hashtable::l2_get_location(uint64_t offset){
  return l2locations[offset];
}

void L2Hashtable::l2_set_location(uint64_t offset, uint32_t location){
  l2locations[offset] = location;
}

void L2Hashtable::l2_write_to_file(const char* name) {
  this->l2_write_hashtable_to_file(name);
  this->l2_write_locations_to_file(name);
  this->l2_write_overflow_to_file(name);
}

/**
 * Writes the L2 tables to disk.
 * This file contains:
 *  l2 threshold
 *  l2 seed size
 *  l2 number of tables
 *  Table Data (all of the tables)
 */
void L2Hashtable::l2_write_hashtable_to_file(const char * name){
  std::string filename = name + l2hashtable_extension;
  FILE *f = fopen(filename.c_str(), "wb");

  // Write how many tables are in this file
  // Tables are of size 6*6*256 = 9216 buckets
  size_t data = fwrite(&l2threshold, sizeof(uint32_t), 1, f);
  if(data != 1) {
    printf("Not enough room on disk\n");
  }
  printf("Writing seed size: %zu\n", l2seed_size);
  data = fwrite(&l2seed_size, sizeof(uint64_t), 1, f);
  if(data != 1){
    printf("Not enough room on disk\n");
  }
  printf("Writing number of tables: %zu\n", l2num_tables);
  data = fwrite(&l2num_tables, sizeof(uint64_t), 1, f);
  if(data != 1){
    printf("Not enough room on disk\n");
  }

  printf("Writing hashtable: total size is %zu\n", l2num_tables*l2table_size);
  data = fwrite(l2hashtable, sizeof(uint32_t), l2num_tables*l2table_size + 1, f);
  printf("Hashtable written.\n");
  if(data != l2table_size * l2num_tables + 1)
    printf("Not enough room on disk. Elements written: %zu/%zu\n", data, l2table_size * l2num_tables);
  fclose(f);
}

void L2Hashtable::l2_write_locations_to_file(const char * name){
  std::string filename = name + l2locations_extension;
  FILE *f = fopen(filename.c_str(), "wb");
  fwrite(&l2locations_size, sizeof(size_t), 1, f);

  size_t data = fwrite(l2locations, sizeof(uint32_t), l2locations_size, f);
  if(data != l2locations_size)
    printf("Not enough room on disk. Elements written: %zu/%zu\n", data, l2locations_size);
  fclose(f);
}

void L2Hashtable::l2_write_overflow_to_file(const char * name){
  std::string filename = name + l2overflow_extension;
  FILE *f = fopen(filename.c_str(), "wb");
  fwrite(&l2overflow_size, sizeof(size_t), 1, f);

  size_t data = fwrite(l2overflow_values, sizeof(uint32_t), l2overflow_size, f);
  if(data != l2overflow_size)
    printf("Overflow table not written to disk\n");
  fclose(f);
}

void L2Hashtable::l2_read_from_file(const char* name) {
  this->l2_read_hashtable_from_file(name);
  this->l2_read_locations_from_file(name);
  this->l2_read_overflow_from_file(name);
}

void L2Hashtable::l2_read_hashtable_from_file(const char * name){
  std::string filename = name + l2hashtable_extension;
  FILE *f = fopen(filename.c_str(), "rb");
  uint32_t threshold[1];
  uint64_t l2seedsize[1];
  uint64_t size[1];

  size_t data = fread(threshold, sizeof(uint32_t), 1, f);
  if(data != 1)
    printf("L2 Threshold not read from file\n");
  this->l2threshold = threshold[0];
  // Get the l2 seed size 
  data = fread(l2seedsize, sizeof(uint64_t), 1, f);
  if(data != 1)
    printf("L2 seed size not read from file\n");
  this->l2seed_size = l2seedsize[0];
  // Get the number of tables
  data = fread(size, sizeof(uint64_t), 1, f);
  if(data != 1)
    printf("Number of tables not read from file\n");
  l2_init_hashtable(size[0]);

  data = fread(l2hashtable, sizeof(uint32_t), l2table_size * l2num_tables + 1, f);
  if(data != l2table_size * l2num_tables + 1)
    printf("Error: Unable to read l2 hashtable\n");
  //for(uint64_t i = 0; i < size[0]*l2table_size+1; i++) {
    //std::cout << l2hashtable[i] << std::endl;
  //}
  fprintf(stderr, "L2 Hashtable.offsets: %p-%p\n", l2hashtable, l2hashtable+l2table_size * l2num_tables + 1);
  fclose(f);
}

void L2Hashtable::l2_read_locations_from_file(const char * name){
  std::string filename = name + l2locations_extension;
  FILE *f = fopen(filename.c_str(), "rb");
  size_t size[1];

  // Get the size of the list
  size_t data = fread(size, sizeof(size_t), 1, f);
  if(data != 1)
    printf("Error: unable to read location list size\n");

  l2_init_locations(size[0]);

  data = fread(l2locations, sizeof(uint32_t), l2locations_size, f);
  if(data != l2locations_size)
    printf("Error: locations list not read\n");
  fprintf(stderr, "L2 Hashtable.locations: %p-%p\n", l2locations, l2locations + l2locations_size);
  fclose(f);
}

void L2Hashtable::l2_read_overflow_from_file(const char * name){
  std::string filename = name + l2overflow_extension;
  FILE *f = fopen(filename.c_str(), "rb");
  size_t size[1];

  // Get overflow size
  size_t data = fread(size, sizeof(size_t), 1, f);
  if(data != 1)
    printf("Error: unable to read overflow size\n");

  //printf("Overflow: %zu\n", size[0]);
  l2_init_overflow(size[0]);

  data = fread(l2overflow_values, sizeof(uint32_t), l2overflow_size, f);
  if(data != l2overflow_size)
    printf("Error: overflow list not read\n");
  fprintf(stderr, "L2 Hashtable.overflow: %p-%p\n", l2overflow_values, l2overflow_values + l2overflow_size);
  fclose(f);
}

void L2Hashtable::l2_free_memory(){
  if(l2locations != 0)
    free(l2locations);
  if(l2hashtable != 0)
    free(l2hashtable);
  if(l2overflow_values != 0)
    free(l2overflow_values);
  l2locations = 0;
  l2hashtable = 0;
  l2overflow_values = 0;
}


/*
 * Constructs the L2 tables
 * */
void L2Hashtable::construct_l2_tables(L2Hashtable & l2hashtable, Hashtable & hashtable, uint32_t big_buckets, std::vector<uint32_t> * buckets, std::vector<std::vector<char>* > * genome, std::vector<uint32_t> ** l1table){
  std::cout << "Starting L2 hashing" << std::endl;
  std::cout << "Elements in L2 hashing: " << big_buckets << std::endl;
  uint32_t l2seed_size = 10;
  uint32_t l1seed_size = 14;
  uint64_t i = 0;

  // Repurpose table for 2nd level hashing.
  // i: 6
  // j: 6
  // buckets: 256
  // big_buckets: number of these tables for L2
  // Each bucket contains a pointer to a vector where the locations will be stored.
  size_t tsize = 6 * 6 * 256;
  std::vector<uint32_t> ** table = (std::vector<uint32_t> **)malloc(tsize * sizeof(std::vector<uint32_t> * ));

  for(i = 0; i < tsize; i++){
    table[i] = 0;
  }

  std::vector<uint32_t> l2temploc;
  l2hashtable.l2_init_hashtable(big_buckets);


  std::cout << "Calculating hash tables for each seed" << std::endl;
  // Calculate hash tables for each location of each seed
  for(i = 0; i < buckets->size(); i++){
    if(i%(buckets->size()>>3) == 0){
      std::cout << (i) << "/" << buckets->size() << std::endl;
    }
    // Initialize all locations to 0
    for(uint64_t ii = 0; ii < tsize; ii++){
      if(table[ii] != 0)
	free(table[ii]);
      table[ii] = 0;
    }

    // Get the current seed, frequency, offset
    uint64_t seed = buckets->at(i);
    uint32_t frequency = l1table[seed]->size(); 

    // Calculations for position in reference genome
    uint32_t index = 0;
    uint32_t max = 0; 
    uint32_t min = 0;
    max += genome->at(index)->size();
    for(uint32_t j = 0; j < frequency; j++){

      // Select a location to use
      uint32_t location = l1table[seed]->at(j);

      // Calculate reference genome location
      for(; index < genome->size(); index++){
	if(min <= location && max > location)
	  break;
	if(index == genome->size() -1){
	  std::cerr << "ERROR: Index out of bounds in reference genome" << std::endl;
	  break;
	}

	// Update the min, max
	min += genome->at(index)->size();
	max += genome->at(index+1)->size();
      }

      // Sanity check: seeds match
      std::string l1string;
      for(uint32_t k = location-min; k < (location-min) + l1seed_size; k++){
	if(k >= genome->at(index)->size())
	  std::cout << "here" << std::endl;
	l1string += genome->at(index)->at(k);
      }
      if(seed != hashtable.get_hash(l1string.c_str()))
	std::cout << "Seeds don't match: " << Hashtable::reverse_hash(seed, hashtable.get_seed_size()) << '\t' << l1string << '\t' << location << std::endl;

      // For each location, populate all 6 I/J combinations
      uint32_t Ioffset = 0;
      for(uint32_t I = 0; I < 7; I++){
	if(I==3){
	  Ioffset = 1;
	  continue;
	}
	for(uint32_t J = 0; J < 6; J++){
	  std::string l2string;
	  uint32_t start; 

	  if(J < I){ // Seeds before I
	    // Check boundaries
	    if(location - min < (I-J) * l2seed_size)
	      continue;

	    start = (location - min) - (I-J) * l2seed_size;
	  }
	  else { // Seeds after I
	    // Check boundaries
	    if(location-min + l1seed_size + (J-I)*l2seed_size + l2seed_size >= genome->at(index)->size())
	      continue;

	    // Find the starting location
	    start = (location-min) + l1seed_size + (J-I)*l2seed_size;
	  }

	  // Get the seed
	  for(uint32_t k = 0; k < l2seed_size; k++){
	    l2string += genome->at(index)->at(start+k);
	  }

	  // Get the hash value and bucket
	  uint8_t hash = l2hashtable.pearsonHash(l2string.c_str());
	  uint32_t bucket = l2hashtable.l2_get_index(0, I-Ioffset, J, hash);
	  if(table[bucket] == 0)
	    table[bucket] = new std::vector<uint32_t>();

	  // Fill the bucket. 
	  table[bucket]->push_back(location);
	}
      }

    }

    uint32_t startindex = l2hashtable.l2_get_index(i, 0, 0, 0);
    for(uint64_t j = 0; j < (36<<8); j++){
      uint32_t frequency = (table[j] == 0) ? 0 : table[j]->size();
      l2hashtable.l2_set_offset(j+startindex, l2temploc.size());
      if(!frequency)
	continue;
      for(auto it = table[j]->begin(); it != table[j]->end(); it++){
	l2temploc.push_back(*it);
	//cout << *it << std::endl;
      }
    }
  }
  // Set very last index.
  l2hashtable.l2_set_offset(l2hashtable.l2_get_index(i, 0, 0, 0), l2temploc.size());

  // Calculate necessary location array size
  std::cout << "Location size: " << l2temploc.size() << std::endl;

  // Initialize l2 tables
  l2hashtable.l2_init_locations(l2temploc.size());
  for(i = 0; i < l2temploc.size(); i++){
    l2hashtable.l2_set_location(i, l2temploc.at(i));
    //cout << check_genome(l2temploc.at(i), genome) << "\n";
  }

  // Calculate the overflow values
  uint32_t previous_value = 0;
  std::vector<uint32_t> overflow_values;
  std::cout << l2hashtable.l2_get_num_tables() << " " << buckets->size() << std::endl;
  for (i = 0; i < l2hashtable.l2_get_num_tables() * l2hashtable.l2_get_table_size() + 1; i++) {
    uint32_t current_value = l2hashtable.l2_get_offset2(i);
    if (current_value < previous_value)
      overflow_values.push_back(i);
    previous_value = current_value;
  }

  l2hashtable.l2_init_overflow(overflow_values.size());
  for (uint32_t i = 0; i < overflow_values.size(); i++) {
    l2hashtable.l2_set_overflow(i, overflow_values.at(i));
  }

  std::cout << "Finished populating l2 tables" << std::endl;

  // Free l2 table used in this function
  free(table);
}
