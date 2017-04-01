#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <set>

#include "hashtable.h"
#include "l2hashtable.h"
using namespace std;

static uint64_t MAX = 0x100000000;
string reverse_hash(uint64_t hash);
string get_seed(uint32_t location, vector<char> * genome);

int main(int argc, char** argv){
  if(argc != 3){
    cout << "usage: ./check <2-level prefix> <big bucket threshold>" << endl;
    exit(1);
  }
  string name = argv[1];
  uint32_t threshold = atoi(argv[2]);
  string h = name + ".hashtable";
  string l = name + ".locations";
  string h2 = name + ".l2.hashtable";
  string l2 = name + ".l2.locations";
  read_table_from_file(h.c_str());
  cout << "read hashtable" << endl;
  read_locations_from_file(l.c_str());
  cout << "read files" << endl;
  if(threshold != 0){
    cout << "Reading l2 files" << endl;
    l2_read_hashtable_from_file(h2.c_str());
    l2_read_locations_from_file(l2.c_str());
    cout << "Read l2 files" << endl;
  }

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
      for(int i = 0; i < line.size(); i++){
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

  uint64_t overflow = 0;

  //for(uint32_t i = 0; i < 1000000000; i++)
    //cout << get_seed(l2_get_location(i), &genome) << endl;

  for(uint32_t i = 0; i < table_size; i++){
    if(i % 10000000 == 0)
      cout << i << "/" << table_size << " verified" << endl;
    // Go over every element in the table
    if(threshold != 0 && get_frequency(i) >= threshold){
      set<uint32_t> loc;
      uint32_t start = l2_get_index(get_offset(i),0,0,0);
      uint32_t end = l2_get_index(get_offset(i),5,5,255);

      uint64_t offset_s = l2_get_offset2(start) + overflow*MAX;
      uint64_t offset_e = l2_get_offset2(end) + overflow*MAX;
      //cout << l2_get_offset2(start) << '\t' << l2_get_offset2(end) << endl;
      if(offset_e < offset_s){
	offset_e += MAX;
	cout << "\nOffset incremented\n" << endl;
	overflow++;
      }

      if(offset_e < offset_s){
	cout << "Problem with offsets" << endl;
      }

      for(uint64_t j = offset_s; j < offset_e; j++){
	loc.insert(l2_get_location(j));
      }

      if(loc.size() != get_frequency(i)){
	cout << "Size mismatch: " << i << '\t' << reverse_hash(i) << '\t' << get_frequency(i) << '\t' << loc.size() << '\t' << (1+offset_e-offset_s)/36 << endl;
	for(auto it = loc.begin(); it != loc.end(); it++){
	  if(get_seed(*it, &genome) == reverse_hash(i))
	    cout << "Location matched!" << endl; 
	}
      }
      else{
	//cout << i << '\t' << reverse_hash(i) << '\t' << get_frequency(i) << endl;// << get_frequency(i) << endl;
	for(auto it = loc.begin(); it != loc.end(); it++){
	  if(get_seed(*it, &genome) != reverse_hash(i))
	    cout << "L2 failure: " << '\t' << i << '\t' << reverse_hash(i) << '\t' << get_seed(*it, &genome) << endl;
	}
      }
    }
    else{ // Level 1 hashing
      uint32_t offset = get_offset(i);
      for(uint32_t j = offset; j < offset + get_frequency(i); j++){
	if(get_seed(get_location(j), &genome) != reverse_hash(i))
	  cout << "L1 Table failure: " << '\t' << i << '\t' << reverse_hash(i) << '\t' << get_seed(get_location(j), &genome) << endl;
      }
    }
  }

  return 0;

}

// Reconstructs the string from the hash number
string reverse_hash(uint64_t hash){
  string reverse;
  for(uint32_t i = 0; i < 14; i++){
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
  for(uint32_t i = 0; i < 14; i++){
    seed += genome->at(location+i);
  }

  return seed;

}
