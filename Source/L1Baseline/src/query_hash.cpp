#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <set>

#include "hashtable.hpp"
using namespace std;

int main(int argc, char** argv) {
  if (argc != 2) {
    cerr << "Usage: query_hash [name of L1Baseline hash table]" << endl;
    cerr << "Terminating." << endl;
  }
  string name = argv[1];
  Hashtable hashtable = Hashtable(name);

/*  vector<char> genome;
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
*/

  cerr << "Size of hash table: " << hashtable.get_table_size() << endl;
  string seed;
  while(1){
    cin >> seed;
    if(seed[0] == 'N')
      break;

    uint64_t hash = hashtable.get_hash(seed.c_str());
    uint32_t frequency = hashtable.get_frequency(hash);
    cout << frequency << endl;
    uint32_t offset = hashtable.get_offset(hash);
    uint32_t end = offset + frequency;
    for(; offset < end; offset++){
      cout << hashtable.get_location(offset) << '\t';
    }
    cout << "\n" << endl;
  }

  return 0;

}
