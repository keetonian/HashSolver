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
void init_table(uint64_t start, uint64_t end, vector<uint32_t> ** table);
void fill_table(uint64_t start, uint64_t end, uint64_t step, uint64_t seed, vector<uint32_t> ** table, vector<char> *genome);

int main(int argc, char** argv){
  cout << "Starting program" << endl;
  if(argc != 5){
    cout << "Usage: ./construct_table <reference genome path> <step size> <seed size> <thread count>" << endl;
    return 0;
  }

  uint64_t num_threads = atoi(argv[4]);

  thread threads[num_threads];

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

  //cout << seed_size << endl;
  //cout << table_s << endl;
  //cout << table_size << endl;

  cout << "seed size: " << seed << "\tthreads: " << num_threads << endl;

  cout << "Table size in Gb: " << ((sizeof(vector<uint32_t>*)*table_s)>>30) << endl;

  //array<vector<uint32_t> * , table_s> table;
  //fill(table.begin(), table.end(), new vector<uint32_t>());
  //vector<uint32_t> * hashtable[table];
  //map<string, uint32_t> names;
  vector<uint32_t> ** table[num_threads];
  for(uint64_t i = 0; i < num_threads; i++){
    table[i] = (vector<uint32_t>**)malloc(sizeof(vector<uint32_t>*)*table_s);
    if(table[i] == 0){
      cerr << "Unable to allocate sufficient space" << endl;
      return 1;
    }
    //cout << table[i] << endl;
    //cout << "Allocated table space" << endl;
    uint64_t increment = table_s/num_threads;

    for(uint64_t j = 0; j < num_threads; j++){
      if(j == num_threads-1)
	threads[j] = thread(init_table, j*increment, (j+1)*increment + table_s%num_threads, table[i]);
      else
	threads[j] = thread(init_table, j*increment, (j+1)*increment, table[i]);
    }

    for(auto& th: threads) th.join();
  }

  cout << "Initialized table" << endl;

  string line;
  ifstream genome(file);
  vector<char> genome_vector;
  genome_vector.reserve(3500000000);

  if(genome.is_open()) {
    while(getline(genome, line)) {
      if(line.size() == 0)
	continue;
      if(line[0] == '>'){
	cout << line << '\t' << genome_vector.size() << endl;
	continue;
      }

      for(uint32_t j = 0; j < line.size(); j++) {
	genome_vector.push_back(toupper(line[j]));
      }
    }
  }
  else{
    cerr << "Unable to open genome" << endl;
    return 1;
  }

  cout << "Read reference genome" << endl;

  uint64_t genome_size = genome_vector.size();
  uint64_t step_size = genome_size / num_threads;
  uint64_t step_mod = genome_size % num_threads;

  for(uint32_t i = 0; i < num_threads; i++){
    cout << "Thread " << i << " started" << endl;
    if(i == num_threads-1){
      threads[i] = thread(fill_table, i*step_size, (i+1)*step_size + step_mod - seed, step, seed, table[i], &genome_vector);
    }
    else{
      threads[i] = thread(fill_table, i*step_size, (i+1)*step_size, step, seed, table[i], &genome_vector);
    }
  }

  cout << "Joining threads" << endl;
  for(auto& th: threads) th.join();
  

  cout << "Combining tables" << endl;

  for(uint32_t i = 1; i < num_threads; i++){
    for(uint64_t j = 0; j < table_s; j++){
      if(table[i][j] == 0 || table[i][j]->size() == 0)
	continue;
      if(table[0][j] == 0)
	table[0][j] = new vector<uint32_t>();
      table[0][j]->insert(table[0][j]->end(), table[i][j]->begin(), table[i][j]->end());
      free(table[i][j]);
    }
    free(table[i]);
  }

  cout << "Calculating number of locations" << endl;

  uint64_t total_size = 0;
  for(uint64_t i = 0; i < table_s; i++){
    if(table[0][i] != 0)
      total_size += table[0][i]->size();
  }

  cout << "Completed initial table construction" << endl;

  cout << "initializing hashtable" << endl;
  initialize_hashtable();

  cout << "initializing locations array" << endl;
  initialize_location(total_size);

  cout << "populating hashtable and locations" << endl;
  total_size = 0;
  for(uint64_t i = 0; i < table_s; i++){
    uint32_t frequency = table[0][i] == 0? 0 : table[0][i]->size();
    set_frequency(i, frequency);
    set_offset(i, total_size);
    for(uint64_t j = 0; j < frequency; j++){
      set_location(j+total_size, table[0][i]->at(j));
    }
    total_size += frequency;
  }

  cout << "Completed internal table and list construction." << endl;

  stringstream name;
  name << "seed";
  name << seed;
  name << "-";
  name << step;
  string tablename = name.str() + ".hashtable";
  string locname = name.str() + ".locations";
  write_table_to_file(tablename.c_str());
  write_locations_to_file(locname.c_str());

  free_memory();
  free(table[0]);

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

void fill_table(uint64_t start, uint64_t end, uint64_t step, uint64_t seed, vector<uint32_t> ** table, vector<char> *genome){
  //cout << "Starting thread" << endl;
  //cout << seed << endl;
  uint64_t s = start;
  double e = end-start;
  for(; start < end; start += step){
    if((start - s) % 1000000 == 0 && s!=start)
      cout << (start-s)/e << endl;
    stringstream st;
    for(uint64_t j = start; j < seed+start; j++){
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
    table[hash]->push_back(start);
    //cout << "hashed" << endl;
    //mtx.unlock();
  }
}

