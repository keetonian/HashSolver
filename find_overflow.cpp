#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <set>

#include "hashtable.h"
#include "l2hashtable.h"
using namespace std;

int main(int argc, char** argv){
  if(argc != 3){
    cout << "usage: ./check <2-level prefix> <big bucket threshold>" << endl;
    exit(1);
  }
  string name = argv[1];
  uint32_t threshold = atoi(argv[2]);
  string h2 = name + ".l2.hashtable";
  if(threshold != 0){
    cout << "Reading l2 table" << endl;
    l2_read_hashtable_from_file(h2.c_str());
  }

  uint64_t overflow = 0;

  vector<uint32_t> overflows;

  //for(uint32_t i = 0; i < 1000000000; i++)
    //cout << get_seed(l2_get_location(i), &genome) << endl;

  for(uint32_t i = 0; i < l2_get_num_tables() * l2table_size - 1; i++){
    if(i % 10000000 == 0)
      cout << i << "/" << l2_get_num_tables() * l2table_size << " verified" << endl;
    uint64_t offset_s = l2_get_offset2(i);
    uint64_t offset_e = l2_get_offset2(i+1);

    // Need to write this to a file.
    if(offset_e < offset_s){
      overflows.push_back(i+1);
      cout << i+1 << endl;
      cout << "\nOffset incremented\n" << endl;
      overflow++;
    }

  }

  l2_init_overflow(overflows.size());
  cout << overflows.size() << endl;

  for(int i = 0; i < overflows.size(); i++){
    cout << overflows.at(i) << endl;
    l2_set_overflow(i, overflows.at(i));
  }

  string o2 = name + ".l2.overflow";
  l2_write_overflow_to_file(o2.c_str());

  l2_free_memory();

  return 0;

}
