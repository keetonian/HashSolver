#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <set>
#include <map>

#include "l2_hashtable.hpp"
using namespace std;

static uint64_t MAX = 0x100000000;

int main(int argc, char** argv){
  if(argc != 2){
    cout << "usage: ./frequencies <2-level prefix> " << endl;
    exit(1);
  }
  string name = argv[1];
  string h2 = name + ".l2.hashtable";
    cout << "Reading l2 hash table" << endl;
    l2_read_hashtable_from_file(h2.c_str());

  uint64_t overflow = 0;

  //for(uint32_t i = 0; i < 1000000000; i++)
    //cout << get_seed(l2_get_location(i), &genome) << endl;
    //
    map<uint32_t, uint32_t> frequencies;
    cout << l2num_tables << endl;
    cout << l2table_size << endl;

  for(uint32_t i = 0; i < 70816*l2table_size; i++){
    // Go over every element in the table
    frequencies[l2_get_frequency2(i)]++;
  }

  for(auto it = frequencies.begin(); it != frequencies.end(); it++){
    cout << it->first << '\t' << it->second << endl;
  }

  return 0;

}
