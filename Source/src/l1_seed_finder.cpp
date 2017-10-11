#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <set>

using namespace std;

int main(){
  string name = "seed14-1-10N-1023B";

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
      for(unsigned int i = 0; i < line.size(); i++){
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

  cout << genome.size() << endl;
  uint32_t location;
  while(1){
    cin >> location;
    if(location > genome.size())
      break;

    for(uint32_t i = 0; i < 14; i++){
      cout << genome.at(location+i);
    }
    cout << '\n' << endl;
  }

  return 0;

}
