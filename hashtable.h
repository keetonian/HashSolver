#ifndef HASHTABLE_H_
#define HASHTABLE_H_

#include "stdint.h"
#include "stdlib.h"

static size_t table_size = 268435456; // 4^14
static uint64_t seed_size = 14;

typedef struct {
  uint32_t frequency;
  uint32_t offset;
} info;

static info * hashtable = 0; // table that contains frequency/offset data
static uint32_t * locations;       // list that contains seed locations at offsets in hashtable
static size_t locations_size = 0;

const static unsigned char char_values[128] = 
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

#ifdef __cplusplus
extern "C" {
#endif 

  uint64_t get_hash(const char * seed);

  void set_seed_size(uint64_t seed);

  uint32_t get_frequency(uint64_t hash);
  uint32_t get_offset(uint64_t hash);
  info get_info(uint64_t hash);

  void set_frequency(uint64_t hash, uint32_t frequency);
  void set_offset(uint64_t hash, uint32_t offset);
  void set_info(uint64_t hash, info information);

  void initialize_location(uint64_t num_elements);
  void initialize_hashtable();

  void set_location(uint32_t index, uint32_t location);
  uint32_t get_location(uint32_t index);

  void write_table_to_file(const char * filename);
  void read_table_from_file(const char * filename);
  void write_locations_to_file(const char * filename);
  void read_locations_from_file(const char * filename);

  void free_memory();

#ifdef __cplusplus
}
#endif

#endif //HASHTABLE_H_
