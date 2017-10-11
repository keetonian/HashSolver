#include "l2_hashtable.hpp"
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

size_t L2Hashtable::l2_get_locations_size(){
  return l2locations_size;
}

void L2Hashtable::l2_init_hashtable(uint64_t num_tables){
  this->l2num_tables = num_tables;
  this->l2hashtable = (info*)malloc(num_tables*l2table_size*sizeof(info));
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
    if(index > l2overflow_values[i])
      overflow += L2_OVERFLOW;
  }
  return overflow;
}

void L2Hashtable::l2_set_overflow(uint32_t i, uint64_t index){
  l2overflow_values[i] = index;
}

uint32_t L2Hashtable::l2_get_frequency(uint32_t offset, uint32_t I, uint32_t J, uint8_t hash){
  return l2hashtable[l2_get_index(offset, I, J, hash)].frequency;
}

uint64_t L2Hashtable::l2_get_offset(uint32_t offset, uint32_t I, uint32_t J, uint8_t hash){
  uint64_t index = l2_get_index(offset, I, J, hash);
  return l2hashtable[index].offset + l2_get_overflow(index);
}

uint32_t L2Hashtable::l2_get_frequency2(uint64_t index){
  return l2hashtable[index].frequency;
}

uint64_t L2Hashtable::l2_get_offset2(uint64_t index){
  return l2hashtable[index].offset + l2_get_overflow(index);
}

void L2Hashtable::l2_set_frequency(uint64_t index, uint32_t frequency){
  l2hashtable[index].frequency = frequency;
}

void L2Hashtable::l2_set_offset(uint64_t index, uint32_t offset){
  l2hashtable[index].offset = offset;
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
  data = fwrite(l2hashtable, sizeof(info), l2num_tables*l2table_size, f);
  printf("Hashtable written.\n");
  if(data != l2table_size * l2num_tables)
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

  data = fread(l2hashtable, sizeof(info), l2table_size * l2num_tables, f);
  if(data != l2table_size * l2num_tables)
    printf("Error: Unable to read l2 hashtable\n");
  //for(uint64_t i = 0; i < size[0]*l2table_size+1; i++) {
    //std::cout << l2hashtable[i].frequency << std::endl;
  //}
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
