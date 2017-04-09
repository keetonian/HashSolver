#include "hashtable.h"
#include "stdio.h"

uint64_t get_hash(const char * seed){
  uint64_t hash = 0;
  //printf("seed: %zu\n", seed_size);
  uint64_t i = 0;
  for(; i < seed_size; i++)
    hash += (uint64_t)(char_values[seed[i]]) << (2ULL*((seed_size-i)-1));
  return hash;
}

void set_seed_size(uint64_t seed){
  seed_size = seed;
  table_size = 4ULL << (2*(seed-1));
}

uint64_t get_seed_size(){
  return seed_size;
}

uint64_t get_table_size(){
  return table_size;
}

uint32_t get_frequency(uint64_t hash){
  return hashtable[hash].frequency;
}

uint32_t get_offset(uint64_t hash){
  return hashtable[hash].offset;
}

info get_info(uint64_t hash){
  return hashtable[hash];
}

void set_frequency(uint64_t hash, uint32_t frequency){
  hashtable[hash].frequency = frequency;
}

void set_offset(uint64_t hash, uint32_t offset){
  hashtable[hash].offset = offset;
}

void set_info(uint64_t hash, info information){
  hashtable[hash] = information;
}
    
// Need to check for improper sizes
void initialize_location(size_t num_elements){ 
  if(locations != 0){
    free(locations);
  }
  locations_size = num_elements;
  locations = (uint32_t *)malloc(locations_size * sizeof(uint32_t));
}

void initialize_hashtable(){
  if(hashtable != 0){
    free(hashtable);
  }
  hashtable = (info*)malloc(table_size * sizeof(info));
}

void set_location(uint32_t index, uint32_t location){
  locations[index] = location;
}

uint32_t get_location(uint32_t index){
  return locations[index];
}

void write_table_to_file(const char * filename){
  FILE *f = fopen(filename, "wb");
  // Write the table seed size first:
  size_t data = fwrite(&seed_size, sizeof(uint64_t), 1, f);
  if(data != 1){
    printf("Not enough room on disk\n");
  }

  data = fwrite(hashtable, sizeof(info), table_size, f);
  if(data != table_size)
    printf("Elements written: %zu/%zu\n", data, table_size);
  fclose(f);
}

void read_table_from_file(const char * filename){
  FILE *f = fopen(filename, "rb");
  uint64_t size[1];
  size_t data = fread(size, sizeof(uint64_t), 1, f);
  if(data != 1)
    printf("Seed size not read from the table\n");
  set_seed_size(size[0]);

  initialize_hashtable();
  data = fread(hashtable, sizeof(info), table_size, f);
  if(data != table_size)
    printf("Elements read: %zu/%zu\n", data, table_size);
  fclose(f);
}

void free_memory(){
  if(hashtable != 0)
    free(hashtable);
  if(locations != 0)
    free(locations);
}

void write_locations_to_file(const char * filename){
  FILE *f = fopen(filename, "wb");
  fwrite(&locations_size, sizeof(size_t), 1, f);
  //printf("%zu\n", locations_size);

  // Write the size of the list as the first element
  size_t data = fwrite(locations, sizeof(uint32_t), locations_size, f);

  if(data != locations_size)
    printf("Elements written: %zu/%zu\n", data, locations_size);
  fclose(f);
}

void read_locations_from_file(const char * filename){
  FILE *f = fopen(filename, "rb");
  size_t size[1];

  // Get the size of the list (first element)
  size_t data = fread(size, sizeof(size_t), 1, f);
  if(data != 1)
    printf("Size of locations not read \n");

  initialize_location(size[0]);

  data = fread(locations, sizeof(uint32_t), locations_size, f);
  if(data != locations_size)
    printf("Not enough memory. Elements read: %zu/%zu\n", data, locations_size);
  fclose(f);
}
