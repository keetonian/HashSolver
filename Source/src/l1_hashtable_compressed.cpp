#include "l1_hashtable.hpp"
#include <iostream>
#include "stdio.h"

Hashtable::Hashtable() {
  this->set_seed_size(14);
}

Hashtable::Hashtable(uint64_t seed) {
  if(seed > 17 || seed < 10){
    std::cerr << "Proper seed size is between 10 and 17 bp, inclusive" << std::endl;
    exit(1);
  }
  this->set_seed_size(seed);
}

Hashtable::Hashtable(std::string filename) {
  this->read_from_file(filename.c_str());
}

Hashtable::~Hashtable() {
  this->free_memory();
}

std::string Hashtable::get_name() {
  return "Compressed";
}

void Hashtable::set_l2_hashtable(L2Hashtable * l2table) {
  this->l2hashtable = l2table;
}

uint64_t Hashtable::get_hash(const char * seed) {
  uint64_t hash = 0;
  //printf("seed: %zu\n", seed_size);
  uint64_t i = 0;
  for(; i < seed_size; i++)
    hash |= (uint64_t)(char_values[(unsigned char)(seed[i])]) << (2ULL*((seed_size-i)-1));
  return hash;
}

uint64_t Hashtable::get_hash(char c, uint64_t prev_hash) {
  prev_hash &= ~((uint64_t)0x3 << (2ULL * (seed_size-1) ) );
  prev_hash = prev_hash << 2;
  prev_hash |= (uint64_t)(char_values[(unsigned char)c]);
  return prev_hash;
}


void Hashtable::set_seed_size(uint64_t seed) {
  seed_size = seed;
  table_size = 4ULL << (2*(seed-1));
  offsets_size = table_size;
  // 1/3 the size of the normal table, make sure to add the remainder to the size.
  frequencies_size = (offsets_size + (3-(offsets_size%3)))/3;
}

uint64_t Hashtable::get_seed_size() {
  return seed_size;
}

size_t Hashtable::get_table_size() {
  return table_size;
}

size_t Hashtable::get_frequencies_size() {
  return frequencies_size;
}

size_t Hashtable::get_locations_size() {
  return locations_size;
}

const uint64_t MASK = 0x1FFFFF;

uint32_t Hashtable::get_frequency(uint64_t hash){
  return (frequencies[hash/3]>>(21*(hash%3))) & MASK;
}

uint32_t Hashtable::get_offset(uint64_t hash){
  return offsets[hash];
}

void Hashtable::set_frequency(uint64_t hash, uint64_t frequency){
  uint64_t old_value = frequencies[hash/3] & ~(MASK<<(21*(hash%3)));
  uint64_t new_value = (frequency<<(21*(hash%3))) ^ old_value;
  frequencies[hash/3]= new_value;
}

void Hashtable::set_offset(uint64_t hash, uint32_t offset){
  offsets[hash]= offset;
}

// Need to check for improper sizes
void Hashtable::initialize_location(size_t num_elements){ 
  if(locations != 0){
    free(locations);
  }
  locations_size = num_elements;
  locations = (uint32_t *)malloc(locations_size * sizeof(uint32_t));
}

void Hashtable::initialize_hashtable(){
  if(frequencies != 0){
    free(frequencies);
  }
  if(offsets != 0) {
    free(offsets);
  }

  frequencies = (uint64_t*)malloc(frequencies_size * sizeof(uint64_t));
  offsets = (uint32_t*)malloc(offsets_size * sizeof(uint32_t));
}

void Hashtable::set_location(uint32_t index, uint32_t location){
  locations[index] = location;
}

uint32_t Hashtable::get_location(uint32_t index){
  return locations[index];
}

void Hashtable::write_to_file(const char * name) {
  this->write_frequencies_to_file(name);
  this->write_offsets_to_file(name);
  this->write_locations_to_file(name);
}

void Hashtable::read_from_file(const char * name) {
  this->read_frequencies_from_file(name);
  this->read_offsets_from_file(name);
  this->read_locations_from_file(name);
}

void Hashtable::write_frequencies_to_file(const char * name){
  std::string filename = name + frequencies_extension;
  FILE *f = fopen(filename.c_str(), "wb");
  // Write the table seed size first:
  size_t data = fwrite(&seed_size, sizeof(uint64_t), 1, f);
  if(data != 1){
    printf("Not enough room on disk\n");
  }

  data = fwrite(frequencies, sizeof(uint64_t), frequencies_size, f);
  if(data != frequencies_size)
    printf("Elements written: %zu/%zu\n", data, frequencies_size);
  fclose(f);
}

void Hashtable::read_frequencies_from_file(const char * name){
  std::string filename = name + frequencies_extension;
  FILE *f = fopen(filename.c_str(), "rb");
  uint64_t size[1];
  size_t data = fread(size, sizeof(uint64_t), 1, f);
  if(data != 1)
    printf("Seed size not read from the table\n");
  set_seed_size(size[0]);

  initialize_hashtable();
  data = fread(frequencies, sizeof(uint64_t), frequencies_size, f);
  if(data != frequencies_size)
    printf("Elements read: %zu/%zu\n", data, frequencies_size);
  fprintf(stderr, "L1 Hashtable.frequencies: %p-%p\n", frequencies, frequencies+frequencies_size);
  fclose(f);
}

void Hashtable::write_offsets_to_file(const char * name){
  std::string filename = name + offsets_extension;
  FILE *f = fopen(filename.c_str(), "wb");
  size_t data = fwrite(offsets, sizeof(uint32_t), offsets_size, f);
  if(data != offsets_size)
    printf("Elements written: %zu/%zu\n", data, offsets_size);
  fclose(f);
}

void Hashtable::read_offsets_from_file(const char * name){
  std::string filename = name + offsets_extension;
  FILE *f = fopen(filename.c_str(), "rb");
  size_t data = fread(offsets, sizeof(uint32_t), offsets_size, f);
  if(data != offsets_size)
    printf("Elements read: %zu/%zu\n", data, offsets_size);
  fprintf(stderr, "L1 Hashtable.offsets: %p-%p\n", offsets, offsets+offsets_size);
  fclose(f);
}

void Hashtable::free_memory(){
  if(frequencies != 0)
    free(frequencies);
  if(offsets != 0)
    free(offsets);
  if(locations != 0)
    free(locations);
  locations = 0;
  frequencies = 0;
  offsets = 0;
}

void Hashtable::write_locations_to_file(const char * name){
  std::string filename = name + locations_extension;
  FILE *f = fopen(filename.c_str(), "wb");
  fwrite(&locations_size, sizeof(size_t), 1, f);
  //printf("%zu\n", locations_size);

  // Write the size of the list as the first element
  size_t data = fwrite(locations, sizeof(uint32_t), locations_size, f);

  if(data != locations_size)
    printf("Elements written: %zu/%zu\n", data, locations_size);
  fclose(f);
}

void Hashtable::read_locations_from_file(const char * name){
  std::string filename = name + locations_extension;
  FILE *f = fopen(filename.c_str(), "rb");
  size_t size[1];

  // Get the size of the list (first element)
  size_t data = fread(size, sizeof(size_t), 1, f);
  if(data != 1)
    printf("Size of locations not read \n");

  initialize_location(size[0]);

  data = fread(locations, sizeof(uint32_t), locations_size, f);
  if(data != locations_size)
    printf("Not enough memory. Elements read: %zu/%zu\n", data, locations_size);
  fprintf(stderr, "L1 Hashtable.locations: %p-%p\n", locations, locations+locations_size);
  fclose(f);
}

void Hashtable::construct_table(Hashtable *hashtable, uint32_t total_size, std::vector<uint32_t> ** &table) {
  
  hashtable->initialize_hashtable();
  hashtable->initialize_location(total_size);

  total_size = 0;
  std::vector<uint64_t> buckets;
  uint64_t table_size = hashtable->get_table_size();
  for(uint64_t i = 0; i < table_size; i++){
    uint32_t frequency = table[i] == 0 ? 0 : table[i]->size();

    hashtable->set_frequency(i, frequency);
    hashtable->set_offset(i, total_size);
    for(uint64_t j = 0; j < frequency; j++){
      hashtable->set_location(j+total_size, table[i]->at(j));
    }
    total_size += frequency;
  }
}

std::vector<uint32_t> Hashtable::construct_table2(Hashtable *hashtable, uint32_t total_size, std::vector<uint32_t> ** &table, uint32_t l2_threshold) {
  
  hashtable->initialize_hashtable();
  hashtable->initialize_location(total_size);

  total_size = 0;
  std::vector<uint32_t> buckets;
  uint32_t l2_index = 0;
  uint64_t table_size = hashtable->get_table_size();
  for(uint64_t i = 0; i < table_size; i++){
    uint32_t frequency = table[i] == 0 ? 0 : table[i]->size();
    hashtable->set_frequency(i, frequency);

    if(l2_threshold && frequency >= l2_threshold) {
      buckets.push_back(i);
      hashtable->set_offset(i, l2_index);
      l2_index++;
      continue;
    }
    hashtable->set_offset(i, total_size);
    for(uint64_t j = 0; j < frequency; j++){
      hashtable->set_location(j+total_size, table[i]->at(j));
    }
    total_size += frequency;
  }

  return buckets;
}
