#include "hashtable.h"
#include <fstream>
#include <sstream>
#include <array>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <thread>
#include <mutex>
#include <stdlib.h>
//#include <boost/thread.hpp>

using namespace std;

mutex mtx;
char random_base();
void init_table(uint64_t start, uint64_t end, vector<uint32_t> ** table);
void fill_table(uint64_t step, uint64_t seed, vector<uint32_t> ** table, vector<char> *genome);

int main(int argc, char** argv){
  cout << "Starting program" << endl;
  if(argc != 4){
    cout << "Usage: ./construct_table <reference genome path> <step size> <seed size> "<< endl;
    return 0;
  }

  // Future: use the bwa index.
  // Construct from scratch to compare this with
  // future versions using the bwa index.
  char * file = argv[1];
  uint64_t step = atoi(argv[2]);

  size_t seed = atoi(argv[3]);
  const size_t table_s = 4ULL<<((seed-1ULL)*2ULL);
  seed_size = seed;
  table_size = table_s;
  set_seed_size(seed);

  stringstream name;
  name << "seed";
  name << seed;
  name << "-";
  name << step;

  //cout << seed_size << endl;
  //cout << table_s << endl;
  //cout << table_size << endl;

  cout << "seed size: " << seed << endl;

  cout << "Table size in Gb: " << ((sizeof(vector<uint32_t>*)*table_s)>>30) << endl;

  //array<vector<uint32_t> * , table_s> table;
  //fill(table.begin(), table.end(), new vector<uint32_t>());
  //vector<uint32_t> * hashtable[table];
  //map<string, uint32_t> names;
  vector<uint32_t> ** table = (vector<uint32_t> **)malloc(sizeof(vector<uint32_t>*)*table_s);
  if(table == 0){
    cerr << "Unable to allocate sufficient space" << endl;
    return 1;
  }

  // Initialize to 0
  for(uint32_t i = 0; i < table_s; i++){
    table[i] = 0;
  }

  cout << "Initialized table" << endl;

  string line;
  ifstream genome(file);
  string fastafile = name.str() + ".fasta";
  ofstream fasta(fastafile.c_str(), ofstream::out);
  vector<vector<char> * > genome_vector;


  if(genome.is_open()) {
    char previous = 0;
    uint32_t i = 0;
    string last;
    while(getline(genome, line)) {
      //cout << line << endl;
      if(line.size() == 0){
	fasta << line << endl;
	continue;
      }
      if(line[0] == '>'){
	//if(last != "\0")
	  //fasta<<last<<endl;
	cout << line << '\t' << genome_vector.size() << endl;
	genome_vector.push_back(new vector<char>());
	i = genome_vector.size()-1;
	fasta << line << endl;
	previous = 0;
	continue;
      }

      uint32_t j = 0;
      for(; j < line.size(); j++) {
	if(j < 2 && previous){
	  if(previous == 'N' && j == 0 && toupper(line[j]) != 'N'){
	    char base = random_base();
	    genome_vector.at(i)->at(genome_vector.at(i)->size()-1) = base;
	    last[last.size()-1] = base;
	    fasta << last << endl;
	  }
	  else if(previous == 'N' && j == 0 && toupper(line[j]) == 'N'){
	    fasta << last << endl;
	  }
	  else if(previous != 'N' && j == 1 && toupper(line[0]) == 'N' && toupper(line[j]) != 'N')
	  {
	    char base = random_base();
	    genome_vector.at(i)->at(genome_vector.at(i)->size()-1) = base;
	    line[0] = base;
	  }
	}

	if(j >= 2 && j < line.size()-1){
	  if(toupper(line[j]) == 'N' && toupper(line[j+1]) != 'N' && toupper(line[j-1]) != 'N')
	    line[j] = random_base();
	}

	genome_vector.at(i)->push_back(toupper(line[j]));
      }
      //cout << line << endl;
      previous = ((j>=2 && toupper(line[j-2]) == 'N') && (j>=1 && toupper(line[j-1]) == 'N')) ? 0 : toupper(line[j-1]);
      //cout << '\t' << previous << endl;
      if(previous != 'N' || line.size() == 1){
	fasta << line << endl;
	line = "\0";
      }
      else{
	last = line;
      }
    }
  }
  else{
    cerr << "Unable to open genome" << endl;
    return 1;
  }
  genome.close();
  fasta.close();

  cout << "Read reference genome" << endl;

  uint64_t genome_size = genome_vector.size();

  for(uint32_t i = 0; i < genome_size; i++){
    cout << "Chromosome " << i << "/" << genome_size << " started" << endl;
    fill_table(step, seed, table, genome_vector.at(i));
  }

  // Threading removed for now.
  //cout << "Combining tables" << endl;

  /*for(uint32_t i = 1; i < num_threads; i++){
    for(uint64_t j = 0; j < table_s; j++){
    if(table[i][j] == 0 || table[i][j]->size() == 0)
    continue;
    if(table[0][j] == 0)
    table[0][j] = new vector<uint32_t>();
    table[0][j]->insert(table[0][j]->end(), table[i][j]->begin(), table[i][j]->end());
    free(table[i][j]);
    }
    free(table[i]);
    }*/

  cout << "Calculating number of locations" << endl;

  uint64_t total_size = 0;
  for(uint64_t i = 0; i < table_s; i++){
    if(table[i] != 0)
      total_size += table[i]->size();
  }

  cout << "Completed initial table construction" << endl;

  cout << "initializing hashtable" << endl;
  initialize_hashtable();

  cout << "initializing locations array" << endl;
  initialize_location(total_size);

  cout << "populating hashtable and locations" << endl;
  total_size = 0;
  for(uint64_t i = 0; i < table_s; i++){
    uint32_t frequency = table[i] == 0 ? 0 : table[i]->size();
    set_frequency(i, frequency);
    set_offset(i, total_size);
    for(uint64_t j = 0; j < frequency; j++){
      set_location(j+total_size, table[i]->at(j));
    }
    total_size += frequency;
  }

  cout << "Completed internal table and list construction." << endl;

  string tablename = name.str() + ".hashtable";
  string locname = name.str() + ".locations";
  write_table_to_file(tablename.c_str());
  write_locations_to_file(locname.c_str());

  free_memory();
  free(table);

  return 0;
}

void init_table(uint64_t start, uint64_t end, vector<uint32_t> ** table){
  // No possible collisions, so no locks needed
  //cout << "Thread started" << endl;
  for(; start < end; start++){
    //cout << start << endl;
    table[start] = 0;
  }
  //cout << "Thread finished" << endl;
}

void fill_table(uint64_t step, uint64_t seed, vector<uint32_t> ** table, vector<char> *genome){
  //cout << "Starting thread" << endl;
  //cout << seed << endl;
  uint64_t index = 0;
  double e = genome->size();
  for(; index < e-seed; index += step){
    if((index) % 10000000 == 0)
      cout << (index)/e << endl;
    stringstream st;
    for(uint64_t j = index; j < seed+index; j++){
      st << genome->at(j);
    }
    string result = st.str();
    //cout << result << '\t' << result.size() << endl;
    size_t num = result.find("N");
    if(num < seed)
      continue;
    uint64_t hash = get_hash(result.c_str());
    //cout << result << endl;
    //cout << "hash: " << hash << endl;
    //cout << "tabl: " << table_size << endl;
    //mtx.lock();
    if(table[hash] == 0)
      table[hash] = new vector<uint32_t>();
    table[hash]->push_back(index);
    //cout << "hashed" << endl;
    //mtx.unlock();
  }
}

char random_base(){
  uint32_t r = rand() % 4;
  switch(r){
    case 0: return 'A';
    case 1: return 'C';
    case 2: return 'G';
    case 3: return 'T';
  }
  cout << "\n\n\n\n\n########################LOOK#################################\n\n\n\n\n" << endl;
}
