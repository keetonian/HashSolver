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
  if(argc != 5){
    cout << "Usage: ./construct_table <reference genome path> <step size> <seed size> <Replace # consecutive N's>"<< endl;
    return 0;
  }

  // Construct from scratch to compare this with
  // future versions 
  char * file = argv[1];
  uint64_t step = atoi(argv[2]);
  size_t seed = atoi(argv[3]);
  uint32_t replace_n = atoi(argv[4]);
  const size_t table_s = 4ULL<<((seed-1ULL)*2ULL);
  set_seed_size(seed);

  // Set up file name
  stringstream name;
  name << "seed";
  name << seed;
  name << "-";
  name << step;

  cout << "seed size: " << seed << endl;
  cout << "Table size in Gb: " << ((sizeof(vector<uint32_t>*)*table_s)>>30) << endl;

  vector<uint32_t> ** table = (vector<uint32_t> **)malloc(sizeof(vector<uint32_t>*)*table_s);
  if(table == 0){
    cerr << "Unable to allocate sufficient space" << endl;
    return 1;
  }

  // Initialize table buckets to 0
  for(uint32_t i = 0; i < table_s; i++){
    table[i] = 0;
  }

  cout << "Initialized table" << endl;

  // Prepare to parse reference genome
  string line;
  ifstream genome(file);
  vector<vector<char> * > genome_vector;
  vector<string> genome_names;
  if(genome.is_open()) {
    uint32_t i = 0;
    while(getline(genome, line)) {

      // Sometimes contains blank lines in between sections
      if(line.size() == 0){
	continue;
      }

      // Fasta format: name lines start with a > character
      if(line[0] == '>'){
	cout << line << endl;
	genome_vector.push_back(new vector<char>());
	genome_names.push_back(line);
	i = genome_vector.size()-1;
	continue;
      }

      // Read all the DNA into a vector for future use
      uint32_t j = 0;
      for(; j < line.size(); j++) {
	genome_vector.at(i)->push_back(toupper(line[j]));
      }
    }
  }
  else{
    cerr << "Unable to open genome" << endl;
    return 1;
  }
  genome.close();

  cout << "Reference genome read." << endl;
  cout << "Changing strings of N's less than " << replace_n << " characters." << endl;

  string fastafile = name.str() + ".fasta";
  ofstream fasta(fastafile.c_str(), ofstream::out);
  for(uint32_t i = 0; i < genome_vector.size(); i++){
    vector<char> * chromosome = genome_vector.at(i);
    if(i != 0)
      fasta << "\n";
    fasta << genome_names.at(i) << endl;
    for(uint32_t j = 0; j < chromosome->size(); j++){

      // Find N's
      if(chromosome->at(j) == 'N'){
	uint32_t k = j;

	// Find length of string. This does not deal with long strings of N's well.
	while(k < chromosome->size() && chromosome->at(k) == 'N' && k-j != replace_n)
	  k++;

	// Only change N if it is a spurious string of them
	if(k < chromosome->size() && k-j < replace_n && chromosome->at(j-1) != 'N'){
	  for(uint32_t l = j; l < k; l++)
	    chromosome->at(l) = random_base();
	}
      }
      
      // Add bases to fasta file, with only 60 bases on a line.
      if(j != 0 && j % 60 == 0)
	fasta << "\n";
      fasta << chromosome->at(j); 
    }
  }
  fasta.close();

  cout << "N's changed" << endl;

  uint64_t genome_size = genome_vector.size();

  for(uint32_t i = 0; i < genome_size; i++){
    cout << "Chromosome " << i << "/" << genome_size << " started" << endl;
    fill_table(step, seed, table, genome_vector.at(i));
  }

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
  cout << "Writing table and list to file" << endl;

  string tablename = name.str() + ".hashtable";
  string locname = name.str() + ".locations";
  write_table_to_file(tablename.c_str());
  write_locations_to_file(locname.c_str());

  cout << "Freeing memory" << endl;
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

    size_t num = result.find("N");
    if(num < seed)
      continue;
    uint64_t hash = get_hash(result.c_str());

    if(table[hash] == 0)
      table[hash] = new vector<uint32_t>();
    table[hash]->push_back(index);
  }
}

// G++ gives a warning about this function. It shouldn't.
char random_base(){
  uint32_t r = rand() % 4;
  switch(r){
    case 0: return 'A';
    case 1: return 'C';
    case 2: return 'G';
    case 3: return 'T';
  }
}
