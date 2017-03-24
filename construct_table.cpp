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
//#include <boost/thread.hpp>

using namespace std;

mutex mtx;
void init_table(unsigned int start, unsigned int end, vector<uint32_t> ** table);
void fill_table(unsigned int start, unsigned int end, unsigned int step, unsigned int seed, vector<uint32_t> ** table, vector<char> *genome);

int main(int argc, char** argv){
  cout << "Starting program" << endl;
  if(argc != 3){
    cout << "Usage: ./construct_table <reference genome path> <step size>" << endl;
    return 0;
  }

  thread threads[32];

  // Future: use the bwa index.
  // Construct from scratch to compare this with
  // future versions using the bwa index.
  char * file = argv[1];
  unsigned int step = atoi(argv[2]);

  const size_t table_s = 268435456;
  size_t seed = 14;


  //array<vector<uint32_t> * , table_s> table;
  //fill(table.begin(), table.end(), new vector<uint32_t>());
  //vector<uint32_t> * hashtable[table];
  //map<string, uint32_t> names;
  vector<uint32_t> ** table = (vector<uint32_t>**)malloc(table_s*sizeof(vector<uint32_t>*));

  unsigned int increment = table_s/32;

  for(int i = 0; i < 32; i++)
    threads[i] = thread(init_table, i*increment, (i+1)*increment, table);

  for(auto& th: threads) th.join();

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

      for(int j = 0; j < line.size(); j++) {
        genome_vector.push_back(toupper(line[j]));
      }
    }
  }
  else{
    cerr << "Unable to open genome" << endl;
    return 1;
  }

  cout << "Read reference genome" << endl;

  /*unsigned int genome_size = genome_vector.size();
  unsigned int step_size = genome_size / 32;
  unsigned int step_mod = genome_size % 32;

  for(int i = 0; i < 32; i++){
    if(i == 31){
      threads[i] = thread(fill_table, i*step_size, (i+1)*step_size + step_mod - seed, step, seed, table, &genome_vector);
    }
    else{
      threads[i] = thread(fill_table, i*step_size, (i+1)*step_size, step, seed, table, &genome_vector);
    }
  }

  for(auto& th: threads) th.join();*/

  unsigned long total_size = 0;

  for(unsigned int i = 0; i < genome_vector.size() - seed; i+=step){
    if(i%10000000 == 0)
      cout << i << "/" << genome_vector.size() << endl;
    stringstream s;
    for(unsigned int j = 0; j < seed; j++){
      s<<genome_vector.at(i+j);
    }
    string result = s.str();
    size_t num = result.find("N");
    if(num < seed)
      continue;
    //cout << num << endl;
    //cout << '\t' << result << endl;
    uint32_t hash = get_hash(result.c_str());
    //cout << hash << endl;
    table[hash]->push_back(i);
    total_size++;
  }

  cout << "Completed initial table construction" << endl;

  initialize_hashtable();
  initialize_location(total_size);

  total_size = 0;
  for(unsigned int i = 0; i < table_s; i++){
    set_frequency(i, table[i]->size());
    set_offset(i, total_size);
    for(unsigned int j = 0; j < table[i]->size(); j++){
      set_location(j+total_size, table[i]->at(j));
    total_size += table[i]->size();
    }
  }

  cout << "Completed internal table and list construction." << endl;

  write_table_to_file("seed14.hashtable");
  write_locations_to_file("seed14.locations");

  free_memory();
  free(table);

  return 0;
}

void init_table(unsigned int start, unsigned int end, vector<uint32_t> ** table){
  // No possible collisions, so no locks needed
  cout << "Thread started" << endl;
  for(; start < end; start++)
    table[start] = new vector<uint32_t>();
}

void fill_table(unsigned int start, unsigned int end, unsigned int step, unsigned int seed, vector<uint32_t> ** table, vector<char> *genome){
  int s = start;
  for(; start < end; start += step){
    if((start - s) % 100000000 == 0 && s!=start)
      cout << "100 million" << endl;
    stringstream s;
    for(unsigned int j = start; j < seed+start; j++){
      s << genome->at(j);
    }
    string result = s.str();
    size_t num = result.find("N");
    if(num < seed)
      continue;
    uint32_t hash = get_hash(result.c_str());
    mtx.lock();
    table[hash]->push_back(start);
    mtx.unlock();
  }
}



