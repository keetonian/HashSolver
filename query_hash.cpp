#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <set>

#include "hashtable.h"
#include "l2hashtable.h"
using namespace std;

int main(){
  string name = "seed14-1-10N-0B";
  string h = name + ".hashtable";
  string l = name + ".locations";
  read_table_from_file(h.c_str());
  cout << "read hashtable" << endl;
  read_locations_from_file(l.c_str());
  cout << "read files" << endl;

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
    uint32_t offset = get_offset(hash);
    uint32_t end = offset + frequency;
    for(; offset < end; offset++){
      cout << get_location(offset) << '\t';
    }
    cout << "\n" << endl;
  }

  cout << total << endl;

  return 0;

}
