#include <fstream>
#include <sstream>
#include <array>
#include <iostream>
#include <map>
#include <stdlib.h>
#include <mutex>
#include <tbb/concurrent_vector.h>

#include "l1_hashtable.hpp"
#include "l2_hashtable.hpp"
#include "genome.hpp"
#include "l2_construct_table.hpp"
#include "l2_commandline.hpp"

using namespace std;

int main(int argc, char** argv){
  int a = parseCommands(argc, argv);
  if(a || !argc || !seed || !genome_file){
    print_options();
    exit(0);
  }

  //TODO: update this to be more automatic
  //tbb::task_scheduler_init init(4);

  hashtable = Hashtable(seed);
  l2hashtable = L2Hashtable(l2_threshold);
  const size_t table_size = hashtable.get_table_size();
  create_name(hashtable.get_name(), l2hashtable.get_name());

  // Set up the table
  cout << "seed size: " << seed << endl;
  cout << "Table size in Gb: " << ((sizeof(vector<uint32_t>*)*table_size)>>30) << endl;
  vector<uint32_t> ** table = (vector<uint32_t> **)malloc(sizeof(vector<uint32_t>*)*table_size);
  if(table == 0){
    cerr << "Unable to allocate sufficient space" << endl;
    return 1;
  }

  mutex **mtxs;
#if defined (USE_MULTITHREADING)
  //cout << "Making mutexes" << endl;
  //cout << table_size/100 << endl;
  mtxs = (mutex**)malloc(sizeof(mutex*) * table_size/100 + 1);
  for(uint32_t i = 0; i < 1+table_size/100; i++){
    mutex *m = new mutex();
    mtxs[i] = m;
  }
  //cout << "Done" << endl;
#endif

  // Initialize table buckets 
  if(seed <= 14){
    for(uint64_t i = 0; i < table_size; i++){
      table[i] = new vector<uint32_t>(0);
      if(i < 50 || i+50 > table_size){
	if(i == 1 || i == table_size-1){
	  table[i]->reserve(100000);
	}
	else
	  table[i]->reserve(400);
      }
    }
  }
  else{
    for(uint64_t i = 0; i < table_size; i++){
      table[i] = 0;
    }
  }

  // Load the reference genome
  initialize_genome();

  // Fill the L1 table
  uint64_t offset = 0;
  cout << "Constructing table" << endl;
  for(uint32_t i = 0; i < genome_size; i++){
    cout << "Chromosome " << (i+1) << "/" << genome_size << " started" << endl;
    fill_table(seed, table, genome_vector.at(i), offset, mtxs);
    offset+=genome_vector.at(i)->size();
  }

  cout << "Calculating number of locations" << endl;

  uint64_t total_size = 0;
  for(uint64_t i = 0; i < table_size; i++){
    if(table[i] != 0)
      if(l2_threshold == 0 || table[i]->size() < l2_threshold)
	total_size += table[i]->size();
  }

  cout << "Completed initial table construction" << endl;

  cout << "initializing hashtable" << endl;
  hashtable.initialize_hashtable();

  cout << "initializing locations array" << endl;
  hashtable.initialize_location(total_size);

  cout << "populating hashtable and locations" << endl;
  total_size = 0;
  uint32_t big_buckets = 0;
  vector<uint64_t> buckets;
  for(uint64_t i = 0; i < table_size; i++){
    uint32_t frequency = table[i] == 0 ? 0 : table[i]->size();
    hashtable.set_frequency(i, frequency);

    // Save big buckets, keep track of how many there are.
    if(l2_threshold != 0 && frequency >= l2_threshold){
      big_buckets++;
      hashtable.set_offset(i, buckets.size());
      buckets.push_back(i);
      continue;
    }

    hashtable.set_offset(i, total_size);
    for(uint64_t j = 0; j < frequency; j++){
      hashtable.set_location(j+total_size, table[i]->at(j));
    }
    total_size += frequency;
  }

  cout << "Completed internal table and list construction." << endl;
  cout << "Writing table and list to file" << endl;

  hashtable.write_to_file(name.c_str());

  cout << "Freeing l1 tables" << endl;
  hashtable.free_memory();

  if(buckets.size()){
    cout << "L2 Hashing" << endl;
    construct_l2_table(big_buckets, &buckets, &genome_vector, table);
    l2hashtable.l2_write_to_file(name.c_str());
    l2hashtable.l2_free_memory();
  }

  cout << "Freeing hashtable memory" << endl;
  for(uint64_t i = 0; i < table_size; i++) {
    if(table[i] != 0)
      delete table[i];
  }
  cout << "Freeing table" << endl;
  free(table);

  cout << "Finished" << endl;
  return 0;
}

void create_name(string l1_hashtable_name, string l2_hashtable_name){
  // Set up file name
  stringstream namestream;
  namestream << directory << "/L1" << l1_hashtable_name << "-L2" << l2_hashtable_name << "-";

  namestream << seed << "S-";
  namestream << replace_n << "N-";
  namestream << contigs << "C-";
  namestream << l2_threshold << "T";
  name = namestream.str();
}

void initialize_genome(){
  // Prepare reference genome
  cout << "Reading reference genome" << endl;
  genome = new Genome(genome_file, contigs);
  cout << "Reference genome read. Changing strings of N's less than " << replace_n << " characters." << endl;
  int changed = genome->change_all_N(replace_n);

  string fastaname = name;
  fastaname += ".fasta";
  if(write_genome_to_file && changed){
    cout << "Characters replaced. Writing genome to file: " << fastaname << endl;
    genome->write_genome(fastaname.c_str(), contigs);
  }

  genome_vector = *(genome->get_genome());
  genome_names = *(genome->get_genome_names());
  genome_size = genome_vector.size();

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
void fill_table(uint64_t seed, vector<uint32_t> ** table, vector<char> *genome, uint64_t offset, mutex ** mtxs){
#if defined (USE_MULTITHREADING)
  tbb::parallel_for(tbb::blocked_range<size_t>(0, genome->size()-seed), L1Constructor(seed, table, genome, offset, mtxs));
#else
  uint64_t index = 0;
  double e = genome->size();

  for(; index < e-seed; index += 1){
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
    uint64_t hash = hashtable.get_hash(result.c_str());

    if(table[hash] == 0){
      mtx.lock();
      table[hash] = new vector<uint32_t>();
      mtx.unlock();
    }
    table[hash]->push_back(index+offset);
  }
#endif
}


/*
 * Constructs the L2 tables
 * */
void construct_l2_table(uint32_t big_buckets, vector<uint64_t> * buckets, vector<vector<char>* > * genome, vector<uint32_t> ** l1table){
  cout << "Starting L2 hashing" << endl;
  cout << "Elements in L2 hashing: " << big_buckets << endl;
  uint32_t l2seed_size = 10;
  uint32_t l1seed_size = 14;
  uint64_t i = 0;

  // Repurpose table for 2nd level hashing.
  // i: 6
  // j: 6
  // buckets: 256
  // big_buckets: number of these tables for L2
  // Each bucket contains a pointer to a vector where the locations will be stored.
  size_t tsize = 6 * 6 * 256;
  vector<uint32_t> ** table = (vector<uint32_t> **)malloc(tsize * sizeof(vector<uint32_t> * ));

  for(i = 0; i < tsize; i++){
    table[i] = 0;
  }

  vector<uint32_t> l2temploc;
  l2hashtable.l2_init_hashtable(big_buckets);


  cout << "Calculating hash tables for each seed" << endl;
  // Calculate hash tables for each location of each seed
  for(i = 0; i < buckets->size(); i++){
    if(i%(buckets->size()>>3) == 0){
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
      if(seed != hashtable.get_hash(l1string.c_str()))
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
	  uint8_t hash = l2hashtable.pearsonHash(l2string.c_str());
	  uint32_t bucket = l2hashtable.l2_get_index(0, I-Ioffset, J, hash);
	  if(table[bucket] == 0)
	    table[bucket] = new vector<uint32_t>();

	  // Fill the bucket. 
	  table[bucket]->push_back(location);
	}
      }

    }

    uint32_t startindex = l2hashtable.l2_get_index(i, 0, 0, 0);
    for(uint64_t j = 0; j < (36<<8); j++){
      uint32_t frequency = (table[j] == 0) ? 0 : table[j]->size();
      l2hashtable.l2_set_offset(j+startindex, l2temploc.size());
      if(!frequency)
	continue;
      for(auto it = table[j]->begin(); it != table[j]->end(); it++){
	l2temploc.push_back(*it);
	//cout << *it << endl;
      }
    }
  }
  // Set very last index.
  l2hashtable.l2_set_offset(l2hashtable.l2_get_index(i, 0, 0, 0), l2temploc.size());

  // Calculate necessary location array size
  cout << "Location size: " << l2temploc.size() << endl;

  // Initialize l2 tables
  l2hashtable.l2_init_locations(l2temploc.size());
  for(i = 0; i < l2temploc.size(); i++){
    l2hashtable.l2_set_location(i, l2temploc.at(i));
    //cout << check_genome(l2temploc.at(i), genome) << "\n";
  }

  // Calculate the overflow values
  uint32_t previous_value = 0;
  vector<uint32_t> overflow_values;
  cout << l2hashtable.l2_get_num_tables() << " " << buckets->size() << endl;
  for (i = 0; i < l2hashtable.l2_get_num_tables() * l2hashtable.l2_get_table_size() + 1; i++) {
    uint32_t current_value = l2hashtable.l2_get_offset2(i);
    if (current_value < previous_value)
      overflow_values.push_back(i);
    previous_value = current_value;
  }

  l2hashtable.l2_init_overflow(overflow_values.size());
  for (uint32_t i = 0; i < overflow_values.size(); i++) {
    l2hashtable.l2_set_overflow(i, overflow_values.at(i));
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

