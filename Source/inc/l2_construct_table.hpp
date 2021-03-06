#ifndef L2_CONSTRUCT_TABLE_H_
#define L2_CONSTRUCT_TABLE_H_

#include <string>
#include <vector>
#include "stdint.h"
#include "genome.hpp"
#include "l1_hashtable.hpp"
#include "l2_hashtable.hpp"
#include "tbb/tbb.h"

#define PLACE1(src, dest) ((src) | dest)
#define PLACE2(src, dest) ((src<<21) | dest)
#define PLACE3(src, dest) ((src<<42) | dest)
#define GET_PLACE1(src) (src&(~0x1FFFFF))
#define GET_PLACE2(src) ((src&(~(0x1FFFFF<<21)))>>21)
#define GET_PLACE3(src) ((src&(~(0x1FFFFF<<42)))>>42)

char* genome_file=0;
char* directory = 0;
bool write_genome_to_file = true;
size_t seed=0;
uint32_t replace_n=10;
uint32_t l2_threshold=0;
Genome *genome;
Hashtable hashtable;
L2Hashtable l2hashtable;
int contigs = 0;
std::string name;
std::vector<std::vector<char> * > genome_vector;
std::vector<std::string> genome_names;
uint64_t genome_size;
std::mutex mtx;


void init_table(uint64_t start, uint64_t end, std::vector<uint32_t> ** table);
void fill_table(uint64_t seed, std::vector<uint32_t> ** table, std::vector<char> *genome, uint64_t offset, std::mutex ** mtxs);
void construct_l2_table(uint32_t big_buckets, std::vector<uint64_t> * buckets, std::vector<std::vector<char>* > * genome, std::vector<uint32_t> ** l1table);
uint32_t get_genome_index(uint32_t location, std::vector<std::vector<char> * > * genome);
std::string check_genome(uint32_t location, std::vector<std::vector<char> * > * genome);
std::string reverse_hash(uint64_t hash);
void initialize_genome();
void create_name(std::string l1_hashtable_name, std::string l2_hashtable_name);

class L1Constructor {
  uint64_t my_seed;
  std::vector<uint32_t> ** my_table;
  std::vector<char> * const my_genome;
  uint64_t my_offset;
  std::mutex ** my_mutex;

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
      uint64_t hash = hashtable.get_hash(result.c_str());
      my_mutex[hash/100]->lock();
      if(my_table[hash] == 0){
	my_table[hash] = new std::vector<uint32_t>();
      }
      my_table[hash]->push_back(i+my_offset);
      my_mutex[hash/100]->unlock();
    }
  }

  L1Constructor(uint64_t seed, std::vector<uint32_t> ** table, std::vector<char> * const genome, uint64_t offset, std::mutex ** const mtxs):
    my_seed(seed), my_table(table), my_genome(genome), my_offset(offset), my_mutex(mtxs){}

};

#endif //L2_CONSTRUCT_TABLE_H_
