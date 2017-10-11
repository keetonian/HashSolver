#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <set>

#include "l1_hashtable.hpp"
#include "l2_hashtable.hpp"
using namespace std;

string reverse_hash(uint64_t hash);
string get_seed(uint32_t location, vector<char> * genome);
uint64_t seed_size;

int main(int argc, char** argv){
  if(argc != 2){
    cout << "usage: ./check <2-level prefix>" << endl;
    exit(1);
  }
  string name = argv[1];
  Hashtable hashtable = Hashtable(name);
  L2Hashtable l2hashtable = L2Hashtable(name);

  seed_size = hashtable.get_seed_size();

  vector<char> genome;
  string g = name + ".fasta";
  ifstream fasta(g.c_str());
  if(fasta.is_open()){
    string line;
    while(getline(fasta, line)){
      if(line.size() == 0)
	continue;
      if(line[0] == '>')
	continue;
      for(uint32_t i = 0; i < line.size(); i++){
	genome.push_back(toupper(line[i]));
      }
    }
  }
  else{
    cerr << "unable to open genome" << endl;
    return 1;
  }
  fasta.close();

  cout << "Genome read" << endl;

  for(uint64_t i = 0; i < hashtable.get_table_size(); i++){
    if(i % (hashtable.get_table_size()>>4) == 0)
      cout << i << "/" << hashtable.get_table_size() << " verified" << endl;
    // Go over every element in the table
    if(l2hashtable.l2_get_threshold() != 0 && hashtable.get_frequency(i) >= l2hashtable.l2_get_threshold()){
      set<uint32_t> loc;
      uint32_t index_start = l2hashtable.l2_get_index(hashtable.get_offset(i),0,0,0);
      uint32_t index_end = l2hashtable.l2_get_index(hashtable.get_offset(i),5,5,255);

      uint64_t offset_start = l2hashtable.l2_get_offset2(index_start);
      uint64_t offset_end = l2hashtable.l2_get_offset2(index_end);
      //cout << l2hashtable.l2_get_offset2(index_start) << '\t' << l2hashtable.l2_get_offset2(index_end) << endl;

      if(offset_end < offset_start){
	cout << "Problem with offsets" << endl;
      }

      offset_end += l2hashtable.l2_get_frequency2(index_end);

      for(uint64_t j = offset_start; j < offset_end; j++){
	//cout << j << endl;
	loc.insert(l2hashtable.l2_get_location(j));
      }

      if(loc.size() != hashtable.get_frequency(i)){
	cout << "Size mismatch: " << i << '\t' << reverse_hash(i) << '\t' << hashtable.get_frequency(i) << '\t' << loc.size() << '\t' << (1+offset_end-offset_start)/36 << endl;
	for(auto it = loc.begin(); it != loc.end(); it++){
	  if(get_seed(*it, &genome) == reverse_hash(i))
	    cout << "Location matched!" << endl; 
	}
      }
      else{
	//cout << i << '\t' << reverse_hash(i) << '\t' << hashtable.get_frequency(i) << endl;// << hashtable.get_frequency(i) << endl;
	for(auto it = loc.begin(); it != loc.end(); it++){
	  if(get_seed(*it, &genome) != reverse_hash(i))
	    cout << "L2 failure: " << '\t' << i << '\t' << reverse_hash(i) << '\t' << get_seed(*it, &genome) << endl;
	}
      }
    }
    else{ // Level 1 hashing
      uint32_t offset = hashtable.get_offset(i);
      for(uint32_t j = offset; j < offset + hashtable.get_frequency(i); j++){
	if(get_seed(hashtable.get_location(j), &genome) != reverse_hash(i))
	  cout << "L1 Table failure: " << '\t' << i << '\t' << reverse_hash(i) << '\t' << get_seed(hashtable.get_location(j), &genome) << endl;
      }
    }
  }

  return 0;

}

// Reconstructs the string from the hash number
string reverse_hash(uint64_t hash){
  string reverse;
  for(uint32_t i = 0; i < seed_size; i++){
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

string get_seed(uint32_t location, vector<char> * genome){
  string seed;
  for(uint32_t i = 0; i < seed_size; i++){
    seed += genome->at(location+i);
  }

  return seed;

}
