#include <vector>
#include <set>
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
#include "sw_aligner.hpp"
#include "shd_filter.hpp"
#include <chrono>
#include "edlib.h"
#include "mapper_common.hpp"
#include "opal_aligner.hpp"


using namespace std;

int main(int argc, char** argv) {

  // Parse Command Line Parameters
  int command_line = parseCommands(argc, argv);
  if(command_line || !hashtable_filename || !reads_filename ) {
    print_options();
    exit(1);
  }

  // Load hash table
  hashtable.read_from_file(hashtable_filename);

  // Load genome
  read_genome_2bit();
  read_genome_char();
  cout << "genome read" << endl;

  // Select Solver
  switch(seed_selection) {
    case SeedSelection::naive: solver = new BasicSolver();
			       solver->loadHashTables(&hashtable);
			       break;
    case SeedSelection::hobbes: 
			       solver = new HobbesSolver();
			       solver->loadHashTables(&hashtable);
			       break;
    case SeedSelection::optimal:
			       // Optimal solver.
			       // Uses bwa instead of hash table.
			       break;
  }

  // Select SWA implementation
  switch(swa_function) {
    case SWAFunction::noswa: finalize_read_locations = &NoSWA;
			     break;
    case SWAFunction::seqalign: finalize_read_locations = &SWA_Seqalign;
				break;
    case SWAFunction::edlib: finalize_read_locations = &Meyers_Edlib;
			     break;
    case SWAFunction::opal: finalize_read_locations = &Opal;
			    break;
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

  // Prepare Smith Waterman step
  swa_threshold = read_length - error_threshold;

  // Initialize the read information for simple batching.
  ReadInformation * reads = (ReadInformation*)malloc(sizeof(ReadInformation) * group);
  for (uint32_t i = 0; i < group; i++) {
    ReadInformation read;
    read.seeds = new uint8_t[number_of_seeds];
    reads[i] = read;
  }

  int number_completed = 0;

  double time_seeds=0, time_locations=0, time_filter=0, time_swa=0;
  // PIPELINE
  do {

    // Get batch of reads
    uint32_t i = 0;
    do {
      if (pre_filter(ks->seq.s)) {
	reads[i].read = new std::string(ks->seq.s);
	i++;
      }
    } while(i < group && kseq_read(ks) >= 0);

    cout << number_completed << endl;
    // What about doing it all backwards?
    // seed selection
    auto start = chrono::steady_clock::now();
    get_seeds(reads, i);
    auto end = chrono::steady_clock::now();
    time_seeds += (end - start).count();

    // location selection
    start = chrono::steady_clock::now();
    get_locations(reads, i);
    end = chrono::steady_clock::now();
    time_locations += (end - start).count();

    // filtering
    start = chrono::steady_clock::now();
    filter_reads(reads, i);
    end = chrono::steady_clock::now();
    time_filter += (end - start).count();

    // SWA
    start = chrono::steady_clock::now();
    if(do_swa)
      finalize_read_locations(reads, i);
    end = chrono::steady_clock::now();
    time_swa += (end - start).count();

    // free memory
    free_read_memory(reads, i);
    number_completed++;

  } while (kseq_read(ks) >= 0);


  cerr << "Done" << endl;
  cerr << "Timing: " << endl;
  cerr << "Seeds: " << time_seeds << endl;
  cerr << "Locations: " << time_locations << endl;
  cerr << "Filters: " << time_filter << endl;
  cerr << "SWA: " << time_swa << endl;
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
  free(genome_char);
  hashtable.free_memory();

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
  uint32_t location;

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
	// Subtract seed position to make sure locations start at beginning of read
	location = hashtable.get_location(offset+k);
	locations[index] = location - seed;
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
    locations = reads[i].locations;
    //cout << reads[i].frequency << endl;
    for (uint32_t j = 0; j < reads[i].frequency; j++) {
      // filter. If something does not pass the filter, 
      // then make the location 0, decrement frequency.
      decompress_2bit_dna(reference, locations[j]);
      //cout << "\t" << locations[j] << "\n";
      //cout << "\t" << reference << "\n";
      //cout << "\t" << *(reads[i].read) << "\n\n";
      if(!shd_filter.filter(reads[i].read->c_str(), reference, error_threshold, read_length)){
	locations[j] = 0;
	reads[i].frequency--;
      }

    }
  }
}

// Returns false if read is
bool pre_filter(string read) {
  int count = 0;
  for(uint32_t i = 0; i < read_length; i++) {
    if (read[i] == 'N'){
      count++;
      if (count > error_threshold)
	return false;
    }
  }
  return true;
}

void NoSWA(ReadInformation * reads, uint32_t number_of_reads) {
  char reference[read_length];
  uint32_t * locations;
  uint32_t index;
  for (uint32_t i = 0; i < number_of_reads; i++) {
    locations = reads[i].locations;
    cout << *(reads[i].read) << " " << reads[i].frequency << "\n";
    index = 0;
    for (uint32_t j = 0; j < reads[i].frequency; j++) {
      while (!locations[index]) { index++; }
      decompress_2bit_dna(reference, locations[index]);
      cout << locations[index] << " ";
      index++;
    }
    cout << "\n";
  }
}

void SWA_Seqalign(ReadInformation * reads, uint32_t number_of_reads) {
  char reference[read_length];
  uint32_t * locations;
  uint32_t index;
  for (uint32_t i = 0; i < number_of_reads; i++) {
    locations = reads[i].locations;
    cout << *(reads[i].read) << " " << reads[i].frequency << "\n";
    index = 0;
    for (uint32_t j = 0; j < reads[i].frequency; j++) {
      while (!locations[index]) { index++; }
      decompress_2bit_dna(reference, locations[index]);
      if (swaligner.sw_align(reads[i].read->c_str(), reference, swa_threshold))
	cout << locations[index] << " ";
      index++;
    }
    cout << "\n";
  }
}

void Meyers_Edlib(ReadInformation * reads, uint32_t number_of_reads) {
  char reference[read_length];
  uint32_t * locations;
  uint32_t index;
  for (uint32_t i = 0; i < number_of_reads; i++) {
    locations = reads[i].locations;
    cout << *(reads[i].read) << " " << reads[i].frequency << "\n";
    index = 0;
    for (uint32_t j = 0; j < reads[i].frequency; j++) {
      while (!locations[index]) { 
	index++;
      }
      decompress_2bit_dna(reference, locations[index]);
      EdlibAlignResult result = edlibAlign(reads[i].read->c_str(), read_length, reference, read_length, edlibNewAlignConfig(error_threshold, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
      if (result.editDistance != -1 ){
	cout << locations[index] << " ";
      }
      edlibFreeAlignResult(result);
      index++;
    }
    cout << "\n";
  }
}

void convert_read(string read, unsigned char * destination, uint32_t read_length) {
  for (uint32_t i = 0; i < read_length; i++) {
    destination[i] = (unsigned char)char_values[read[i]];
  }
}

void Opal(ReadInformation * reads, uint32_t number_of_reads) {
  // Need a function to convert read to 0,1,2,3
  // Need a function to convert the reference to 0,1,2,3 (unsigned char)
  int reference_length = 0;
  unsigned char read[read_length];
  uint32_t * locations;
  uint32_t index;
  for (uint32_t i = 0; i < number_of_reads; i++) {
    unsigned char ** reference = (unsigned char **)malloc(reads[i].frequency * sizeof(unsigned char *));
    int * reference_lengths = (int*)malloc(reads[i].frequency * sizeof(int));
    locations = reads[i].locations;
    cout << *(reads[i].read) << " " << reads[i].frequency << "\n";
    index = 0;

    convert_read(*(reads[i].read), read, read_length);
    reference_length = reads[i].frequency;

    for (uint32_t j = 0; j < reads[i].frequency; j++) {
      while (!locations[index]) { index++; }
      //decompress_2bit_dna(reference, locations[index]);
      reference[j] = genome_char + locations[index];
      reference_lengths[j] = read_length;
      index++;
    }
    opal_aligner.opal_swa(read, read_length, reference, reference_length, reference_lengths);
    free(reference);
    free(reference_lengths);
  }
}

// Read up: see if this is strictly necessary.
void free_read_memory(ReadInformation * reads, uint32_t number_of_reads) {
  for (uint32_t i = 0; i < number_of_reads; i++) {
    free(reads[i].locations);
    free(reads[i].read);
  }
}

void read_genome_char() {
  // size: 3.2 billion / 4 Bytes (approx 0xc0000000 >> 5 => 0x6000000)
  // WRONG! FIX ME.
  genome_char = (unsigned char *)malloc(0xc0000000);
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
	//if (index % 32 == 0)
	//genome[index>>5] = 0;
	// 32 bases can fit in a 64-bit number.
	// Index>>5 is Index/32
	// the other math converts the character into a 2-bit base,
	//  then places it in the correct spot inside of the number.
	// Or or Xor? The debate still rages.
	genome_char[index] = (unsigned char)char_values[line[i]];
	index++;
      }
    }
  }
  else {
    cerr << "Unable to open genome." << endl;
  }
  //cout << index << endl;
  fasta.close();
}

void read_genome_2bit() {
  // size: 3.2 billion / 4 Bytes (approx 0xc0000000 >> 5 => 0x6000000)
  // WRONG! FIX ME.
  genome = (uint64_t *)malloc(0x7000000*sizeof(uint64_t));
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
	//if (index % 32 == 0)
	//genome[index>>5] = 0;
	// 32 bases can fit in a 64-bit number.
	// Index>>5 is Index/32
	// the other math converts the character into a 2-bit base,
	//  then places it in the correct spot inside of the number.
	// Or or Xor? The debate still rages.
	genome[index>>5] |= (((uint64_t)char_values[line[i]])<<((index%32)*2));
	index++;
      }
    }
  }
  else {
    cerr << "Unable to open genome." << endl;
  }
  cout << index << endl;
  //cout << index << endl;
  fasta.close();
}

void decompress_2bit_dna(char * destination, uint32_t starting_index) {
  for (uint32_t i = 0; i < read_length; i++) {
    uint64_t dna = genome[(starting_index+i)>>5];
    uint8_t char_code = (dna >> (2*((starting_index + i)%32))) & 0x03;
    destination[i] = reverse_values[char_code];
  }
}
