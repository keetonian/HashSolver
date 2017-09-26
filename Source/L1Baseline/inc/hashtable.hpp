#ifndef HASHTABLE_H_
#define HASHTABLE_H_

#include "stdint.h"
#include "stdlib.h"
#include <string>
#include "hash_common.hpp"

class Hashtable{
  protected:

    size_t table_size;
    uint64_t seed_size;

    info * hashtable = 0;	    // table that contains frequency/offset data
    uint32_t * locations = 0;       // list that contains seed locations at offsets in hashtable
    size_t locations_size = 0;

  public:
    Hashtable();
    Hashtable(uint64_t seed_size);
    Hashtable(std::string filename);
    ~Hashtable();

    uint64_t get_hash(const char * seed);

    void set_seed_size(uint64_t seed);
    uint64_t get_seed_size();
    uint64_t get_table_size();
    size_t get_locations_size();

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

    void write_to_file(const char * name);
    void read_from_file(const char * name);

    void write_table_to_file(const char * filename);
    void read_table_from_file(const char * filename);
    void write_locations_to_file(const char * filename);
    void read_locations_from_file(const char * filename);

    void free_memory();

};
#endif //HASHTABLE_H_
