#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <set>

#include "hashtable.h"
#include "l2hashtable.h"
using namespace std;

int main(){
  string name = "seed14-1-10N-1000B";
  string h = name + ".hashtable";
  string l = name + ".locations";
  string h2 = name + ".l2.hashtable";
  string l2 = name + ".l2.locations";
  read_table_from_file(h.c_str());
  cout << "read hashtable" << endl;
  read_locations_from_file(l.c_str());
  cout << "read files" << endl;
  l2_read_hashtable_from_file(h2.c_str());
  l2_read_locations_from_file(l2.c_str());
  cout << "Read l2 files" << endl;

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

  uint64_t total = 0;
  cout << table_size << endl;
  for(uint32_t i = 0; i < table_size; i++){
    //cout << i << '\t' << get_frequency(i) << '\t' << get_offset(i);
    if(get_frequency(i) > 1023){
      uint32_t offset = get_offset(i);
      //cout << offset << endl;
      uint32_t index = l2_get_index(offset, 0, 0, 0);
      uint32_t index2 = l2_get_index(offset+1, 0, 0, 0);
      set<uint32_t> loc;
      uint32_t l2frequency2 = 0;
      for(uint32_t j = index; j < index + l2table_size; j++){
	l2frequency2 += l2_get_frequency2(j);
	//cout << l2_get_frequency2(j) << endl;
	if(l2_get_frequency2(j) > 0){
	  uint32_t l2offset = l2_get_offset2(j);
	}
      }

      for(uint32_t jj = l2_get_offset2(index); jj <= l2_get_offset(offset, 5, 5, 255)+1; jj++){
	loc.insert(l2_get_location(jj));
      }

      cout << index2 - l2_get_index(offset, 5, 5, 255) << endl;

      cout << get_frequency(i) << '\t' << loc.size() << endl;
      uint32_t l2frequency = l2_get_offset2(index2) - l2_get_offset2(index);
      if(get_frequency(i) != l2frequency/36 && l2frequency != l2frequency2)
	cout << get_frequency(i) << '\t' << l2frequency << '\t' << l2frequency2 << endl;
      else
	cout << "Match" << endl;
      // Now look into l2 hash
      //l2
    }
    //cout << endl;
  }

  cout << total << endl;

  return 0;

}
