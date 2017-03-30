#include "l2hashtable.h"
#include "stdio.h"

uint8_t pearsonHash(const char * seed) {
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
  size_t j;
  uint8_t h=0, index;

  for(i=0; i<l2seed_size; i++) {
    index = h^seed[i];
    h = T[index];
  }

  return h;
}

void l2_set_num_tables(uint64_t number){
  l2num_tables = number;
}

void l2_init_hashtable(uint64_t num_tables){
  l2num_tables = num_tables;
  l2hashtable = (info*)malloc(num_tables*l2table_size*sizeof(info));
}

void l2_init_locations(uint64_t size){
  l2locations_size = size;
  l2locations = (uint32_t*)malloc(l2locations_size * sizeof(uint32_t));
}

uint32_t l2_get_index(uint32_t offset, uint32_t I, uint32_t J, uint8_t hash){
  return ((I*6)<<8) + ((J)<<8) + ((offset * 36)<<8) + hash;
}

uint32_t l2_get_frequency(uint32_t offset, uint32_t I, uint32_t J, uint8_t hash){
  return l2hashtable[l2_get_index(offset, I, J, hash)].frequency;
}

uint32_t l2_get_offset(uint32_t offset, uint32_t I, uint32_t J, uint8_t hash){
  return l2hashtable[l2_get_index(offset, I, J, hash)].offset;
}

uint32_t l2_get_frequency2(uint32_t index){
  return l2hashtable[index].frequency;
}

uint32_t l2_get_offset2(uint32_t index){
  return l2hashtable[index].offset;
}

void l2_set_frequency(uint32_t index, uint32_t frequency){
  l2hashtable[index].frequency = frequency;
}

void l2_set_offset(uint32_t index, uint32_t offset){
  l2hashtable[index].offset = offset;
}

uint32_t l2_get_location(uint32_t offset){
  return l2locations[offset];
}

void l2_set_location(uint32_t offset, uint32_t location){
  l2locations[offset] = location;
}

void l2_write_hashtable_to_file(const char * name){
  FILE *f = fopen(name, "wb");

  // Write how many tables are in this file
  // Tables are of size 6*6*256 = 9216 buckets
  printf("Writing number of tables: %zu\n", l2num_tables);
  size_t data = fwrite(&l2num_tables, sizeof(uint64_t), 1, f);
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

void l2_write_locations_to_file(const char * name){
  FILE *f = fopen(name, "wb");
  fwrite(&l2locations_size, sizeof(size_t), 1, f);

  size_t data = fwrite(l2locations, sizeof(uint32_t), l2locations_size, f);
  if(data != l2locations_size)
    printf("Not enough room on disk. Elements written: %zu/%zu\n", data, l2locations_size);
  fclose(f);
}

void l2_read_hashtable_from_file(const char * name){
  FILE *f = fopen(name, "rb");
  uint64_t size[1];

  // Get the number of tables
  size_t data = fread(size, sizeof(uint64_t), 1, f);
  if(data != 1)
    printf("Number of tables not read from file\n");
  l2_init_hashtable(size[0]);

  data = fread(l2hashtable, sizeof(info), l2table_size * l2num_tables, f);
  if(data != l2table_size * l2num_tables)
    printf("Error: Unable to read l2 hashtable\n");
  fclose(f);
}

void l2_read_locations_from_file(const char * name){
  FILE *f = fopen(name, "rb");
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

void l2_free_memory(){
  if(l2locations != 0)
    free(l2locations);
  if(l2hashtable != 0)
    free(l2hashtable);
}
