#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <set>

#include "hashtable.h"
#include "l2hashtable.h"
using namespace std;

int main(){
  string name = "seed14-1-10N-1023B";
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
  string seed;
  while(1){
    cin >> seed;
    if(seed[0] == 'N')
      break;

    uint64_t hash = get_hash(seed.c_str());
    uint32_t frequency = get_frequency(hash);
    cout << frequency << endl;
    if(frequency > 1000){
      //l2
      uint32_t start = l2_get_index(get_offset(hash),0,0,0);
      uint32_t end = l2_get_index(get_offset(hash)+1,0,0,0);
      cout << "Large Bucket Encountered. Enter I: " << endl;
      uint32_t I;
      cin >> I;
      for(uint32_t J = 0; J < 6; J++){
	cout << "Enter 10bp seed for J" << J << endl;
	string l2seed;
	cin >> l2seed;
	uint8_t l2_hash = pearsonHash(l2seed.c_str());
	uint64_t index = l2_get_index(get_offset(hash), I, J, l2_hash);
	uint32_t l2freq = l2_get_frequency2(index);
	uint32_t l2offset = l2_get_offset2(index);

	for(uint32_t k = l2offset; k < l2offset+l2freq; k++){
	  cout << l2_get_location(k) << "\t";
	}
	cout << "\n" << endl;
      }

      //for(; start < end; start++){
	//cout << l2_get_location(start) << '\t';
      //}
      cout << "\n" << endl;
    }
    else{
      uint32_t offset = get_offset(hash);
      uint32_t end = offset + frequency;
      for(; offset < end; offset++){
	cout << get_location(offset) << '\t';
      }
      cout << "\n" << endl;
    }
  }

  cout << total << endl;

  return 0;

}
