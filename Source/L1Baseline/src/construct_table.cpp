#include <fstream>
#include <sstream>
#include <array>
#include <iostream>
#include <map>
#include <stdlib.h>
#include <mutex>

#include "hashtable.hpp"
#include "genome.hpp"
#include "construct_table.hpp"
#include "commandline.hpp"

using namespace std;

int main(int argc, char** argv){
  int a = parseCommands(argc, argv);
  if(a || !argc || !seed || !file){
    print_options();
    exit(0);
  }

  //TODO: update this to be more automatic
  //tbb::task_scheduler_init init(4);

  hashtable = Hashtable(seed);
  const size_t table_size = hashtable.get_table_size();
  create_name();

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
  vector<uint64_t> buckets;
  for(uint64_t i = 0; i < table_size; i++){
    uint32_t frequency = table[i] == 0 ? 0 : table[i]->size();

    hashtable.set_frequency(i, frequency);
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

void create_name(){
  // Set up file name
  stringstream namestream;
  namestream << "/tmp/DNA/seed";

  namestream << seed;
  namestream << "-";
  namestream << replace_n;
  namestream << "N-";
  namestream << l2_threshold;
  namestream << "B";
  name = namestream.str();
}

void initialize_genome(){
  // Prepare reference genome
  cout << "Reading reference genome" << endl;
  genome = new Genome(file, contigs);
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

