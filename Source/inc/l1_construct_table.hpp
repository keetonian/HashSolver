#ifndef L1_CONSTRUCT_TABLE_H_
#define L1_CONSTRUCT_TABLE_H_

#include <string>
#include <vector>
#include "stdint.h"
#include "genome.hpp"
#include "l1_hashtable.hpp"
#include "tbb/tbb.h"

#define PLACE(src, dest, mod) ((src<<(21*mod)) | dest)
#define GET_PLACE(src) ((src>>(21*mod))&0x1FFFFF)

char* genome_file=0;
char* directory = 0;
bool write_genome_to_file = true;
size_t seed=0;
uint32_t replace_n=10;
uint32_t l2_threshold=0;
Genome *genome;
Hashtable hashtable;
int contigs = 0;
std::string name;
std::vector<std::vector<char> * > genome_vector;
std::vector<std::string> genome_names;
uint64_t genome_size;
std::mutex mtx;


void init_table(uint64_t start, uint64_t end, std::vector<uint32_t> ** table);
void fill_table(uint64_t seed, std::vector<uint32_t> ** table, std::vector<char> *genome, uint64_t offset, std::mutex ** mtxs);
uint32_t get_genome_index(uint32_t location, std::vector<std::vector<char> * > * genome);
std::string check_genome(uint32_t location, std::vector<std::vector<char> * > * genome);
std::string reverse_hash(uint64_t hash);
void initialize_genome();
void create_name(std::string hashtable_name);

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

#endif //L1_CONSTRUCT_TABLE_H_
