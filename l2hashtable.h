#ifndef L2HASHTABLE_H_
#define L2HASHTABLE_H_

#include "stdint.h"
#include "stdlib.h"
#include "hashtable.h"

static uint64_t l2seed_size = 10;
static uint64_t l2num_tables = 0;
static uint64_t l2table_size = (36<<8);
static info * l2hashtable = 0;
static uint32_t * l2locations = 0;
static size_t l2locations_size = 0;

static size_t l2overflow_size = 0;
static uint32_t * l2overflow_values = 0;
const uint64_t L2_OVERFLOW = 0x100000000;

#ifdef __cplusplus
extern "C" {
#endif
  uint8_t pearsonHash(const char * seed);
  void l2_set_num_tables(uint64_t number);
  uint64_t l2_get_num_tables();
  size_t l2_get_locations_size();

  void l2_init_hashtable(uint64_t num_tables);
  void l2_init_locations(uint64_t size);
  void l2_init_overflow(uint64_t size);

  uint32_t l2_get_index(uint32_t offset, uint32_t I, uint32_t J, uint8_t hash);

  uint64_t l2_get_overflow(uint32_t index);
  void l2_set_overflow(uint32_t i, uint32_t index);

  uint32_t l2_get_frequency(uint32_t offset, uint32_t I, uint32_t J, uint8_t hash);
  uint64_t l2_get_offset(uint32_t offset, uint32_t I, uint32_t J, uint8_t hash);
  uint32_t l2_get_frequency2(uint32_t index);
  uint64_t l2_get_offset2(uint32_t index);

  void l2_set_frequency(uint32_t index, uint32_t frequency);
  void l2_set_offset(uint32_t index, uint32_t offset);

  uint32_t l2_get_location(uint64_t offset);
  void l2_set_location(uint64_t offset, uint32_t location);

  void l2_write_hashtable_to_file(const char * name);
  void l2_write_locations_to_file(const char * name);
  void l2_write_overflow_to_file(const char * name);

  void l2_read_hashtable_from_file(const char * name);
  void l2_read_locations_from_file(const char * name);
  void l2_read_overflow_from_file(const char * name);

  void l2_free_memory();

#ifdef __cplusplus
}
#endif

#endif //L2HASHTABLE_H_
