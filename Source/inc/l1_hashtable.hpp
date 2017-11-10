#ifndef L1_HASHTABLE_H_
#define L1_HASHTABLE_H_

#include "stdint.h"
#include "stdlib.h"
#include <string>
#include <vector>
#include "hash_common.hpp"

class L2Hashtable;

class Hashtable{
  protected:

    size_t table_size;
    size_t frequencies_size;
    size_t offsets_size;
    uint64_t seed_size;

    info * hashtable = 0;	    // table that contains frequency/offset data
    uint64_t * frequencies = 0;	    // table that contains frequency data
    uint32_t * offsets = 0;	    // table that contains offset data
    uint32_t * locations = 0;       // list that contains seed locations at offsets in hashtable
    size_t locations_size = 0;

    L2Hashtable * l2hashtable;

  public:
    Hashtable();
    Hashtable(uint64_t seed_size);
    Hashtable(std::string filename);
    ~Hashtable();

    std::string get_name();
    void set_l2_hashtable(L2Hashtable * l2hashtable);
    uint64_t get_hash(const char * seed);
    uint64_t get_hash(char c, uint64_t prev_hash);

    void set_seed_size(uint64_t seed);
    uint64_t get_seed_size();
    uint64_t get_table_size();
    uint64_t get_frequencies_size();
    uint64_t get_offsets_size();
    size_t get_locations_size();

    uint32_t get_frequency(uint64_t hash);
    uint32_t get_offset(uint64_t hash);
    info get_info(uint64_t hash);

    void set_frequency(uint64_t hash, uint64_t frequency);
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
    void write_frequencies_to_file(const char * filename);
    void read_frequencies_from_file(const char * filename);
    void write_offsets_to_file(const char * filename);
    void read_offsets_from_file(const char * filename);
    void write_locations_to_file(const char * filename);
    void read_locations_from_file(const char * filename);

    void free_memory();

    static std::string reverse_hash(uint64_t hash, uint64_t seed_size) {
      std::string reverse;
      for (uint32_t i = 0; i < seed_size; i++){
	uint8_t c = hash & 0x3;
	switch(c){
	  case 0x0: reverse = 'A' + reverse;
		    break;
	  case 0x1: reverse = 'C' + reverse;
		    break;
	  case 0x2: reverse = 'G' + reverse;
		    break;
	  case 0x3: reverse = 'T' + reverse;
		    break;
	}   
	hash = hash >> 2;
      }
      return reverse;
    }

    static void construct_table(Hashtable *hashtable, uint32_t locations_size, std::vector<uint32_t> ** &table);
    static std::vector<uint32_t> construct_table2(Hashtable *hashtable, uint32_t locations_size, std::vector<uint32_t> ** &table, uint32_t l2_threshold);


};
#endif //L1_HASHTABLE_H_
