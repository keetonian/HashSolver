#include <fstream>
#include <sstream>
#include <array>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <stdlib.h>

#include "hashtable.h"
#include "l2hashtable.h"

using namespace std;

char random_base();
void init_table(uint64_t start, uint64_t end, vector<uint32_t> ** table);
void fill_table(uint64_t step, uint64_t seed, vector<uint32_t> ** table, vector<char> *genome, uint64_t offset);
void construct_l2_table(uint32_t big_buckets, vector<uint64_t> * buckets, vector<vector<char> * > * genome, vector<uint32_t> ** l1table);
uint32_t get_genome_index(uint32_t location, vector<vector<char> * > * genome);
string check_genome(uint32_t location, vector<vector<char> * > * genome);
string reverse_hash(uint64_t hash);

int main(int argc, char** argv){
  cout << "Starting program" << endl;
  if(argc != 6){
    cout << "Usage: ./construct_table <reference genome path> <step size> <seed size> <Replace # consecutive N's> <L2 bucket threshold>"<< endl;
    cout << "Set <L2 bucket threshold> to 0 if you do not want any L2 hashing" << endl;
    return 0;
  }

  // Construct from scratch to compare this with
  // future versions 
  char * file = argv[1];
  uint64_t step = atoi(argv[2]);
  size_t seed = atoi(argv[3]);
  uint32_t replace_n = atoi(argv[4]);
  uint32_t bucket_threshold = atoi(argv[5]);
  const size_t table_s = 4ULL<<((seed-1ULL)*2ULL);
  set_seed_size(seed);

  // Set up file name
  stringstream name;
  name << "seed";

  name << seed;
  name << "-";
  name << step;
  name << "-";
  name << replace_n;
  name << "N-";
  name << bucket_threshold;
  name << "B";

  cout << "seed size: " << seed << endl;
  cout << "Table size in Gb: " << ((sizeof(vector<uint32_t>*)*table_s)>>30) << endl;

  vector<uint32_t> ** table = (vector<uint32_t> **)malloc(sizeof(vector<uint32_t>*)*table_s);
  if(table == 0){
    cerr << "Unable to allocate sufficient space" << endl;
    return 1;
  }

  // Initialize table buckets to 0
  for(uint64_t i = 0; i < table_s; i++){
    table[i] = 0;
  }

  cout << "Initialized table" << endl;
  cout << "Reading reference genome" << endl;

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
	char cc = toupper(line[j]);
	if(cc != 'A' && cc != 'C' && cc != 'G' && cc != 'T' && cc != 'N'){
	  cout << "Changing character " << cc << " to N" << endl;
	  cc = 'N';
	}
	genome_vector.at(i)->push_back(cc);
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
	if(k < chromosome->size() && k-j < replace_n && j > 0 && chromosome->at(j-1) != 'N'){
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

  uint64_t offset = 0;
  for(uint32_t i = 0; i < genome_size; i++){
    cout << "Chromosome " << (i+1) << "/" << genome_size << " started" << endl;
    fill_table(step, seed, table, genome_vector.at(i), offset);
    offset+=genome_vector.at(i)->size();
  }

  cout << "Calculating number of locations" << endl;

  uint64_t total_size = 0;
  for(uint64_t i = 0; i < table_s; i++){
    if(table[i] != 0)
      if(bucket_threshold == 0 || table[i]->size() < bucket_threshold)
	total_size += table[i]->size();
  }

  cout << "Completed initial table construction" << endl;

  cout << "initializing hashtable" << endl;
  initialize_hashtable();

  cout << "initializing locations array" << endl;
  initialize_location(total_size);

  cout << "populating hashtable and locations" << endl;
  total_size = 0;
  uint32_t big_buckets = 0;
  vector<uint64_t> buckets;
  for(uint64_t i = 0; i < table_s; i++){
    uint32_t frequency = table[i] == 0 ? 0 : table[i]->size();
    set_frequency(i, frequency);

    // Save big buckets, keep track of how many there are.
    if(bucket_threshold != 0 && frequency >= bucket_threshold){
      big_buckets++;
      set_offset(i, buckets.size());
      buckets.push_back(i);
      continue;
    }

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

  cout << "Freeing l1 tables" << endl;
  free_memory();

  if(buckets.size()){
    cout << "L2 Hashing" << endl;
    construct_l2_table(big_buckets, &buckets, &genome_vector, table);
    string l2hash_name = name.str();
    l2hash_name += ".l2.hashtable";
    string l2loc_name = name.str();
    l2loc_name += ".l2.locations";
    cout << "Writing l2 hash table to file" << endl;
    l2_write_hashtable_to_file(l2hash_name.c_str());
    cout << "Writing l2 locations array to file" << endl;
    l2_write_locations_to_file(l2loc_name.c_str());
  }

  cout << "Freeing hashtable memory" << endl;
  free(table);

  return 0;
}

/*
 * Initializes a table with 0's
 * */
void init_table(uint64_t start, uint64_t end, vector<uint32_t> ** table){
  // No possible collisions, so no locks needed
  //cout << "Thread started" << endl;
  for(; start < end; start++){
    //cout << start << endl;
    table[start] = 0;
  }
  //cout << "Thread finished" << endl;
}

/*
 * Populates the tables with information from a section of the parsed reference.
 * */
void fill_table(uint64_t step, uint64_t seed, vector<uint32_t> ** table, vector<char> *genome, uint64_t offset){
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
    table[hash]->push_back(index+offset);
  }
}

/*
 * Returns a random base
 * */
char random_base(){
  uint32_t r = rand() % 4;
  switch(r){
    case 0: return 'A';
    case 1: return 'C';
    case 2: return 'G';
    case 3: return 'T';
  }

  // Should never reach this, but included here for the compiler
  return 0;
}

/*
 * Constructs the L2 tables
 * */
void construct_l2_table(uint32_t big_buckets, vector<uint64_t> * buckets, vector<vector<char>* > * genome, vector<uint32_t> ** l1table){
  cout << "Starting L2 hashing" << endl;
  cout << "Elements in L2 hashing: " << big_buckets << endl;
  uint32_t l2seed_size = 10;
  uint32_t l1seed_size = 14;

  // Repurpose table for 2nd level hashing.
  // i: 6
  // j: 6
  // buckets: 256
  // big_buckets: number of these tables for L2
  // Each bucket contains a pointer to a vector where the locations will be stored.
  size_t tsize = 6 * 6 * 256;
  vector<uint32_t> ** table = (vector<uint32_t> **)malloc(tsize * sizeof(vector<uint32_t> * ));

  for(uint64_t i = 0; i < tsize; i++){
    table[i] = 0;
  }

  vector<uint32_t> l2temploc;
  l2_init_hashtable(big_buckets);


  cout << "Calculating hash tables for each seed" << endl;
  // Calculate hash tables for each location of each seed
  for(uint32_t i = 0; i < buckets->size(); i++){
    if(i%1000 == 0){
      cout << (i) << "/" << buckets->size() << endl;
    }
    // Initialize all locations to 0
    for(uint64_t ii = 0; ii < tsize; ii++){
      if(table[ii] != 0)
	free(table[ii]);
      table[ii] = 0;
    }

    // Get the current seed, frequency, offset
    uint64_t seed = buckets->at(i);
    uint32_t frequency = l1table[seed]->size(); 

    // Calculations for position in reference genome
    uint32_t index = 0;
    uint32_t max = 0; 
    uint32_t min = 0;
    max += genome->at(index)->size();
    for(uint32_t j = 0; j < frequency; j++){

      // Select a location to use
      uint32_t location = l1table[seed]->at(j);

      // Calculate reference genome location
      for(; index < genome->size(); index++){
	if(min <= location && max > location)
	  break;
	if(index == genome->size() -1){
	  cerr << "ERROR: Index out of bounds in reference genome" << endl;
	  break;
	}

	// Update the min, max
	min += genome->at(index)->size();
	max += genome->at(index+1)->size();
      }

      // Sanity check: seeds match
      string l1string;
      for(uint32_t k = location-min; k < (location-min) + l1seed_size; k++){
	if(k >= genome->at(index)->size())
	  cout << "here" << endl;
	l1string += genome->at(index)->at(k);
      }
      if(seed != get_hash(l1string.c_str()))
	cout << "Seeds don't match: " << reverse_hash(seed) << '\t' << l1string << '\t' << location << endl;

      // For each location, populate all 6 I/J combinations
      uint32_t Ioffset = 0;
      for(uint32_t I = 0; I < 7; I++){
	if(I==3){
	  Ioffset = 1;
	  continue;
	}
	for(uint32_t J = 0; J < 6; J++){
	  string l2string;
	  uint32_t start; 

	  if(J < I){ // Seeds before I
	    // Check boundaries
	    if(location - min < (I-J) * l2seed_size)
	      continue;

	    start = (location - min) - (I-J) * l2seed_size;
	  }
	  else { // Seeds after I
	    // Check boundaries
	    if(location-min + l1seed_size + (J-I)*l2seed_size + l2seed_size >= genome->at(index)->size())
	      continue;

	    // Find the starting location
	    start = (location-min) + l1seed_size + (J-I)*l2seed_size;
	  }

	  // Get the seed
	  for(uint32_t k = 0; k < l2seed_size; k++){
	    l2string += genome->at(index)->at(start+k);
	  }

	  // Get the hash value and bucket
	  uint8_t hash = pearsonHash(l2string.c_str());
	  uint32_t bucket = l2_get_index(0, I-Ioffset, J, hash);
	  if(table[bucket] == 0)
	    table[bucket] = new vector<uint32_t>();

	  // Fill the bucket. 
	  table[bucket]->push_back(location);
	}
      }

    }

    uint32_t startindex = l2_get_index(i, 0, 0, 0);
    for(uint32_t j = 0; j < (36<<8); j++){
      uint32_t frequency = (table[j] == 0) ? 0 : table[j]->size();
      l2_set_frequency(j+startindex, frequency);
      l2_set_offset(j+startindex, l2temploc.size());
      if(!frequency)
	continue;
      for(auto it = table[j]->begin(); it != table[j]->end(); it++){
	l2temploc.push_back(*it);
	//cout << *it << endl;
      }
    }
  }

  // Calculate necessary location array size
  cout << "Location size: " << l2temploc.size() << endl;

  // Initialize l2 tables
  l2_init_locations(l2temploc.size());
  for(uint64_t i = 0; i < l2temploc.size(); i++){
    l2_set_location(i, l2temploc.at(i));
    //cout << check_genome(l2temploc.at(i), genome) << "\n";
  }

  cout << "Finished populating l2 tables" << endl;

  // Free l2 table used in this function
  free(table);
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

string check_genome(uint32_t location, vector<vector<char> * > * genome){
  string line;
  uint32_t min = 0;
  uint32_t max = genome->at(0)->size();
  uint32_t index = 0;
  for(uint32_t i = 0; i < genome->size(); i++){
    if(location < max){
      index = i;
      break;
    }
    min += genome->at(i)->size();
    max += genome->at(i+1)->size();
  }
  for(uint32_t i = location-min; i < location+14-min; i++){
    line += genome->at(index)->at(i);
  }
  return line;
}
