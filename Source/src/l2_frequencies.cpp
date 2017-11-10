#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <set>
#include <map>

#include "l2_hashtable.hpp"
using namespace std;

int main(int argc, char** argv){
  if(argc != 2){
    cout << "usage: ./frequencies <2-level prefix> " << endl;
    exit(1);
  }
  string name = argv[1];
  L2Hashtable l2hashtable = L2Hashtable(name);

  uint64_t overflow = 0;

  //for(uint32_t i = 0; i < 1000000000; i++)
  //cout << get_seed(l2_get_location(i), &genome) << endl;
  //
  map<uint32_t, uint32_t> frequencies;
  cout << l2hashtable.l2_get_num_tables() << endl;
  cout << l2hashtable.l2_get_table_size() << endl;

  for(uint64_t i = 0; i < l2hashtable.l2_get_num_tables()*l2hashtable.l2_get_table_size(); i++){
    // Go over every element in the table
    frequencies[l2hashtable.l2_get_frequency2(i)]++;
  }

  for(auto it = frequencies.begin(); it != frequencies.end(); it++){
    cout << it->first << '\t' << it->second << endl;
  }

  return 0;

}
