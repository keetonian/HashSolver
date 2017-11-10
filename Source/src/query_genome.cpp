#include <vector>
#include <limits.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <set>
#include <stdint.h>

using namespace std;

int main(int argc, char** argv) {
  if (argc != 3) {
    cerr << "Usage: query_genome [prefix of fasta file] [number of bases to print]" << endl;
    cerr << "Terminating." << endl;
    return 1;
  }
  string name = argv[1];
  uint32_t bases = atoi(argv[2]);

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

  uint32_t location;
  uint32_t prev_location = INT_MAX;
  while(1){
    cin >> location;
    if (prev_location == location)
      break;
    prev_location = location;

    for (uint32_t i = 0; i < bases; i++) {
      cout << genome.at(i+location);
    }
    cout << "\n" << endl;
  }

  return 0;

}
