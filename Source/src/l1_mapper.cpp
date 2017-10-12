#include <vector>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <cstring>
#include <zlib.h>
#include <errno.h>
#include <assert.h>
#include <map>
#include "kseq.h"
#include "l1_mapper_commandline.hpp"
#include "l1_mapper.hpp"
#include "l1_hobbes_solver.hpp"
#include "l1_basic_solver.hpp"
#include "seed_solver.hpp"


using namespace std;

int main(int argc, char** argv) {

  // Parse Command Line Parameters
  int command_line = parseCommands(argc, argv);
  if(command_line || !hashtable_filename || !reads_filename || !seed_selection_algorithm) {
    print_options();
    exit(1);
  }

  // Load hash table
  hashtable.read_from_file(hashtable_filename);

  // Load genome
  read_genome();
  cout << "genome read" << endl;

  // Select Solver
  if (seed_selection_algorithm == 0x1) {
    solver = new BasicSolver();
    solver->loadHashTables(&hashtable);
  } else if (seed_selection_algorithm == 0x2) {
    solver = new HobbesSolver();
    solver->loadHashTables(&hashtable);
  } else if (seed_selection_algorithm == 0x4) {
    // Optimal solver.
    // Uses bwa instead of hash table.
  }

  // Load Read File
  gzFile fp;
  kseq_t *ks;
  fp = strcmp(reads_filename, "-")? gzopen(reads_filename, "r") : gzdopen(fileno(stdin), "r");
  if (NULL == fp) {
    fprintf(stderr, "Couldn't open %s : %s\n",
	strcmp(reads_filename, "-") ? reads_filename : "stdin",
	errno ? strerror(errno) : "Out of memory");
    exit(EXIT_FAILURE);
  }
  ks = kseq_init(fp); // initialize the FASTA/Q parser

  kseq_read(ks);

  /*
   * Consider here checking the parameters.
   *	ks-seq.l is the length of the reads.
   *	Check:
   *	  limit: max value is len(read) - number_of_seeds * seed_size
   *	  number_of_seeds: max value is len(read) / seed_size
   *	Don't make everything automatically assume 100 bp reads.
   */

  read_length = ks->seq.l;
  solver->init(read_length, number_of_seeds, hashtable.get_seed_size(), limit);

  // Initialize the read information for simple batching.
  ReadInformation * reads = (ReadInformation*)malloc(sizeof(ReadInformation) * group);
  for (uint32_t i = 0; i < group; i++) {
    ReadInformation read;
    read.seeds = new uint8_t[number_of_seeds];
    reads[i] = read;
  }


  // PIPELINE
  do {

    // Get batch of reads
    uint32_t i = 0;
    do {
      reads[i].read = new std::string(ks->seq.s);
      i++;
    } while(i < group && kseq_read(ks) >= 0);

    // What about doing it all backwards?
    // seed selection
    get_seeds(reads, i);

    // location selection
    get_locations(reads, i);

    // filtering
    filter_reads(reads, i);

    // SWA
    //finalize_read_locations(reads, i);

    // free memory
    free_read_memory(reads, i);

  } while (kseq_read(ks) >= 0);


  cout << "Done" << endl;
  /*
  //int wait;
  //cin >> wait;
  char benchFileName[80];
  strcpy(benchFileName, reads_filename);
  char* benchName = strtok(benchFileName, ".");

  //ofstream output( (string(benchName) + "." + to_string(seedNum) + "-" + to_string(seedLength) + string(".hobbes")) , ofstream::out);

  cout << "seedNum: " << number_of_seeds << " | seedLength: " << hashtable.get_seed_size() << endl;
  //output << "seedNum: " << number_of_seeds << " | seedLength: " << hashtable.get_seed_size() << endl;
  for (map<int, unsigned long long>::iterator iter = frequencyCounter.begin(); iter != frequencyCounter.end(); iter++)
    cout << iter->first << " " << iter->second << endl;


  //output.close();
  */
  kseq_destroy(ks);
  gzclose(fp);
  free(genome);

  return 0;
}

void get_seeds(ReadInformation * reads, uint32_t number_of_reads) {
  for (uint32_t i = 0; i < number_of_reads; i++) {

    int frequency = solver->solveDNA(*(reads[i].read), reads[i].seeds);
    reads[i].frequency = frequency;
  }
}

void get_locations(ReadInformation * reads, uint32_t number_of_reads) {
  // Set up variables
  uint32_t * locations;
  uint8_t * seeds;
  uint32_t index;
  uint32_t seed;
  uint32_t hash;
  uint32_t offset;
  uint32_t frequency;

  // Loop through every read
  for(uint32_t i = 0; i < number_of_reads; i++) {
    reads[i].locations = (uint32_t *)malloc(sizeof(uint32_t) * reads[i].frequency);
    locations = reads[i].locations;
    seeds = reads[i].seeds;
    std::string read_string = *(reads[i].read);
    index = 0;

    // Loop through every seed
    for(uint32_t j = 0; j < number_of_seeds; j++) {
      seed = seeds[j];
      hash = hashtable.get_hash(&read_string[seed]);
      offset = hashtable.get_offset(hash);
      frequency = hashtable.get_frequency(hash);

      for (uint32_t k = 0; k < frequency; k++) {
	// Subtract seed to make sure locations start at beginning of read
	locations[index] = hashtable.get_location(offset+k);// - seed;
	index++;
      }
    }
  }
}

void filter_reads(ReadInformation * reads, uint32_t number_of_reads) {
  char reference[read_length];
  // Implement filters, based on flags set.
  uint32_t * locations;
  for(uint32_t i = 0; i < number_of_reads; i++) {
    cout << i << ", " << reads[i].frequency << endl;
    locations = reads[i].locations;
    for (uint32_t j = 0; j < reads[i].frequency; j++) {
      decompress_2bit_dna(reference, locations[j]);
      cout << reference << "\n";
      cout << *(reads[i].read) << "\n\n";
    }
  }
}

void finalize_read_locations(ReadInformation * reads, uint32_t number_of_reads) {
  //SWA
}

// Read up: see if this is strictly necessary.
void free_read_memory(ReadInformation * reads, uint32_t number_of_reads) {
  for (uint32_t i = 0; i < number_of_reads; i++) {
    free(reads[i].locations);
    free(reads[i].read);
  }
}

void read_genome() {
  // size: 3.2 billion / 4 Bytes (approx 0xc0000000 >> 5 => 0x6000000)
  // WRONG! FIX ME.
  genome = (uint64_t *)malloc(0x6000000);
  string genome_filename = hashtable_filename;
  genome_filename += ".fasta";
  ifstream fasta(genome_filename.c_str());
  uint32_t index = 0;
  if (fasta.is_open()) {
    string line;
    while(getline(fasta, line)) {
      if(line.size() == 0)
	continue;
      if(line[0] == '>')
	continue;
      for (uint32_t i = 0; i < line.size(); i++) {
	// If this is the first for the number, make sure it is 0.
	if (index % 32 == 0)
	  genome[index>>5] = 0;
	// 32 bases can fit in a 64-bit number.
	// Index>>5 is Index/32
	// the other math converts the character into a 2-bit base,
	//  then places it in the correct spot inside of the number.
	// Or or Xor? The debate still rages.
	genome[index>>5] |= (char_values[line[i]]<<((index%32)*2));
      }
    }
  }
  else {
    cerr << "Unable to open genome." << endl;
  }
  fasta.close();
}

void decompress_2bit_dna(char * destination, uint32_t starting_index) {
  for (uint32_t i = 0; i < read_length; i++) {
    uint64_t dna = genome[(starting_index+i)>>5];
    uint8_t char_code = (dna >> (2*((starting_index + i)%32))) & 0x03;
    destination[i] = reverse_values[char_code];
  }
}
