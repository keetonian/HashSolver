#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <set>

#include "l1_hashtable.hpp"
using namespace std;

string reverse_hash(uint64_t hash);
string get_seed(uint32_t location, vector<char> * genome);
uint64_t seed_size;

int main(int argc, char** argv){
  if(argc != 2){
    cout << "usage:<command name> <L1 table name>" << endl;
    exit(1);
  }
  string name = argv[1];
  Hashtable hashtable(name);
  cout << "Read l1 files" << endl;

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

  //for(uint32_t i = 0; i < 1000000000; i++)
  //cout << get_seed(l2_get_location(i), &genome) << endl;

  for(uint64_t i = 0; i < hashtable.get_table_size(); i++){
    if(i % (hashtable.get_table_size()>>4) == 0)
      cout << i << "/" << hashtable.get_table_size() << " verified" << endl;
    // Go over every element in the table
    uint64_t offset = hashtable.get_offset(i);
    //cout << hashtable.get_frequency(i) << endl;
    for(uint64_t j = offset; j < offset + hashtable.get_frequency(i); j++){
      if(get_seed(hashtable.get_location(j), &genome) != reverse_hash(i))
	cerr << "L1 Table failure: " << '\t' << i << '\t' << reverse_hash(i) << '\t' << get_seed(hashtable.get_location(j), &genome) << endl;
    }
  }

  cout << "All locations successfully verified." << endl;

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
