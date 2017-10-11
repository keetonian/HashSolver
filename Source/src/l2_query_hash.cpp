#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <set>

#include "l1_hashtable.hpp"
#include "l2_hashtable.hpp"
using namespace std;

int main(){
  string name = "/tmp/DNA/seed14-10N-1023B";
  cout << "Reading l1 table" << endl;
  Hashtable hashtable = Hashtable(name);
  cout << "Reading l2 table" << endl;
  L2Hashtable l2hashtable = L2Hashtable(name);
  cout << "Tables read." << endl;

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

  uint64_t total = 0;
  cout << hashtable.get_table_size() << endl;
  string seed;
  while(1){
    cin >> seed;
    if(seed[0] == 'N')
      break;

    uint64_t hash = hashtable.get_hash(seed.c_str());
    uint32_t frequency = hashtable.get_frequency(hash);
    cout << frequency << endl;
    if(frequency > 1000){
      // l2
      //uint32_t start = l2hashtable.l2_get_index(hashtable.get_offset(hash),0,0,0);
      //uint32_t end = l2hashtable.l2_get_index(hashtable.get_offset(hash)+1,0,0,0);
      cout << "Large Bucket Encountered. Enter I: " << endl;
      uint32_t I;
      cin >> I;
      for(uint32_t J = 0; J < 6; J++){
	cout << "Enter 10bp seed for J" << J << endl;
	string l2seed;
	cin >> l2seed;
	uint8_t l2_hash = l2hashtable.pearsonHash(l2seed.c_str());
	uint64_t index = l2hashtable.l2_get_index(hashtable.get_offset(hash), I, J, l2_hash);
	uint32_t l2freq = l2hashtable.l2_get_frequency2(index);
	uint32_t l2offset = l2hashtable.l2_get_offset2(index);

	for(uint32_t k = l2offset; k < l2offset+l2freq; k++){
	  cout << l2hashtable.l2_get_location(k) << "\t";
	}
	cout << "\n" << endl;
      }

      //for(; start < end; start++){
	//cout << l2_get_location(start) << '\t';
      //}
      cout << "\n" << endl;
    }
    else{
      uint32_t offset = hashtable.get_offset(hash);
      uint32_t end = offset + frequency;
      for(; offset < end; offset++){
	cout << hashtable.get_location(offset) << '\t';
      }
      cout << "\n" << endl;
    }
  }

  cout << total << endl;

  return 0;

}
