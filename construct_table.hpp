#ifndef CONSTRUCT_TABLE_H_
#define CONSTRUCT_TABLE_H_

#include <string>
#include <vector>
#include "stdint.h"
#include "genome.hpp"
#include "tbb/concurrent_vector.h"
#include "tbb/tbb.h"

#define PLACE1(src, dest) ((src) | dest)
#define PLACE2(src, dest) ((src<<21) | dest)
#define PLACE3(src, dest) ((src<<42) | dest)
#define GET_PLACE1(src) (src&(~0x1FFFFF))
#define GET_PLACE2(src) ((src&(~(0x1FFFFF<<21)))>>21)
#define GET_PLACE3(src) ((src&(~(0x1FFFFF<<42)))>>42)

char* file=0;
bool write_genome_to_file = true;
size_t seed=0;
uint32_t replace_n=10;
uint32_t l2_threshold=0;
Genome *genome;
int contigs = 0;
std::string name;
std::vector<std::vector<char> * > genome_vector;
std::vector<std::string> genome_names;
uint64_t genome_size;
std::mutex mtx;


void init_table(uint64_t start, uint64_t end, std::vector<uint32_t> ** table);
void fill_table(uint64_t seed, tbb::concurrent_vector<uint32_t> ** table, std::vector<char> *genome, uint64_t offset);
void construct_l2_table(uint32_t big_buckets, std::vector<uint64_t> * buckets, std::vector<std::vector<char> * > * genome, tbb::concurrent_vector<uint32_t> ** l1table);
uint32_t get_genome_index(uint32_t location, std::vector<std::vector<char> * > * genome);
std::string check_genome(uint32_t location, std::vector<std::vector<char> * > * genome);
std::string reverse_hash(uint64_t hash);
void initialize_genome();
void create_name();

class L1Constructor {
  uint64_t my_seed;
  tbb::concurrent_vector<uint32_t> ** my_table;
  std::vector<char> * const my_genome;
  uint64_t my_offset;

  public:
  void operator()( const tbb::blocked_range<size_t> & r) const {
    for( size_t i = r.begin(); i != r.end(); ++i){
      std::stringstream st;
      for(uint64_t j = i; j < my_seed + i; j++){
	st << my_genome->at(j);
      }
      std::string result = st.str();
      size_t num = result.find("N");
      if(num < seed)
	continue;
      uint64_t hash = get_hash(result.c_str());
      if(my_table[hash] == 0){
	mtx.lock();
	my_table[hash] = new tbb::concurrent_vector<uint32_t>();
	mtx.unlock();
      }
      my_table[hash]->push_back(i+my_offset);
    }
  }

  L1Constructor(uint64_t seed, tbb::concurrent_vector<uint32_t> ** table, std::vector<char> * const genome, uint64_t offset):
    my_seed(seed), my_table(table), my_genome(genome), my_offset(offset) {}

};

#endif //CONSTRUCT_TABLE_H_
