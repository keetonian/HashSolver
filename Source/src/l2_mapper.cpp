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
#include "mapper_commandline.hpp"
#include "l2_mapper.hpp"
#include "l1_hobbes_solver.hpp"
#include "l1_basic_solver.hpp"
#include "seed_solver.hpp"
#include "sw_aligner.hpp"
#include "shd_filter.hpp"
#include <chrono>
#include "edlib.h"
#include "mapper_common.hpp"
#include "opal_aligner.hpp"
#include "ssw_cpp.h"


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
  l2hashtable.l2_read_from_file(hashtable_filename);

  // Load genome
  read_genome_2bit();
  //read_genome_char();
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
    case SWAFunction::ssw: finalize_read_locations = &SSW;
			    break;
  }

  switch(filter_algorithm) {
    case FilterAlgorithm::none: filter_read_locations = &NoFilter;
				break;
    case FilterAlgorithm::SHD: filter_read_locations = &SHD;
			       break;
    case FilterAlgorithm::MAGNET: filter_read_locations = &MAGNET;
				  break;
    case FilterAlgorithm::QGRAM: filter_read_locations = &QGRAM;
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

  // Find read length
  read_length = ks->seq.l;

  // Init seed solver
  solver->init(read_length, number_of_seeds, hashtable.get_seed_size(), limit);

  // Prepare Smith Waterman step
  swa_threshold = read_length - error_threshold;

  // Initialize the read information for simple batching.

  time_seeds=0; time_locations=0; time_filter=0; time_swa=0;

  // PIPELINE
  do {
    // Pre filter out erroneous reads
    if (!pre_filter(ks->seq.s))
      continue;

    // Map the read
    map_read(ks->seq.s);

  } while (kseq_read(ks) >= 0);


  cerr << "Timing: " << endl;
  cerr << "Seeds: " << time_seeds << endl;
  cerr << "Locations: " << time_locations << endl;
  cerr << "Filters: " << time_filter << endl;
  cerr << "SWA: " << time_swa << endl;

  // Free memory that was used.
  kseq_destroy(ks);
  gzclose(fp);
  free(genome);
  //free(genome_char);
  hashtable.free_memory();

  return 0;
}

// Maps a single read.
void map_read(string read) {

  // Get the reverse read to map as well.
  string reverse = reverse_read(read);

  // seed selection
  auto start = chrono::steady_clock::now();
  uint8_t seeds[number_of_seeds];
  uint8_t reverse_seeds[number_of_seeds];
  solver->solveDNA(read, seeds);
  solver->solveDNA(reverse, reverse_seeds);
  auto end = chrono::steady_clock::now();
  time_seeds += (end - start).count();

  // location selection
  start = chrono::steady_clock::now();
  // Keep locations on the stack unless I have to move them.
  set<uint32_t> locations;
  set<uint32_t> reverse_locations;
  get_locations(&read, seeds, &locations);
  get_locations(&reverse, reverse_seeds, &reverse_locations);
  end = chrono::steady_clock::now();
  time_locations += (end - start).count();

  // filtering and SWA

  cout << read << "\n";
  filter_and_finalize_reads(&read, &locations);
  cout << "Reverse:" << "\n";
  filter_and_finalize_reads(&reverse, &reverse_locations);
}

// Reverses and inverts the string. O(N) complexity.
string reverse_read(string read) {
  string reverse = read;
  for (uint32_t i = 0; i < read_length; i++) {
    reverse[read_length-1-i] = reverse_values[char_values[read[i]] ^ 0x03];
  }
  return reverse;
}

// Uses one of the seed solvers to get the seeds
void get_seeds(string * read, uint8_t * seeds) {
  solver->solveDNA(*read, seeds);
}

// Uses the seeds found using the seed selection step
// Gets all of the locations, puts into a set (no duplicates)
void get_locations(string * read, uint8_t * seeds, set<uint32_t> * locations) {
  // Set up variables
  uint32_t seed;
  uint32_t hash;
  uint32_t offset;
  uint32_t frequency;
  uint64_t l2_offset;
  uint32_t l2_frequency;
  uint32_t location;
  uint32_t l2_seed_size = l2hashtable.l2_get_l2_seed_size();
  const char * rd = read->c_str();

  // Loop through every seed
  for(uint32_t i = 0; i < number_of_seeds; i++) {
    seed = seeds[i];
    hash = hashtable.get_hash(&((*read)[seed]));
    offset = hashtable.get_offset(hash);
    frequency = hashtable.get_frequency(hash);

    if ( frequency < l2hashtable.l2_get_threshold() ){
      // Find all locations for a seed
      for (uint32_t k = 0; k < frequency; k++) {
	location = hashtable.get_location(offset+k);
	// Subtract seed position to make sure locations start at beginning of read
	locations->insert(location - seed);
      }
    }
    else { // FIGURE OUT HOW MANY LOCATIONS L2 TABLE FILTERS
      // 2nd level hash
      int order_fix = 0;
      if(number_of_seeds%2 == 0 && i >= number_of_seeds/2)
	order_fix = 1;
      int idx = (number_of_seeds%2 == 0) ? (number_of_seeds/2 - 1) : number_of_seeds / 2;
      uint32_t j = 0;
      uint32_t jj = 0;
      // Elements before seed: i + order_fix
      for (; j < i + order_fix; j++) {
	uint8_t pearson = l2hashtable.pearsonHash(rd+seed-((j+1)*l2_seed_size));
	hash = l2hashtable.l2_get_index(offset, i, j, pearson);
	l2_offset = l2hashtable.l2_get_offset2(hash);
	l2_frequency = l2hashtable.l2_get_frequency2(hash);
	for (uint32_t l = 0; l < l2_frequency; l++) {
	  location = l2hashtable.l2_get_location(l + l2_offset);
	  locations->insert(location - seed);
	}
      }
      // Elements after seed 6 - i - order_fix
      for (j = 0; j < 6 -i - order_fix; j++) {
	uint8_t pearson = l2hashtable.pearsonHash(rd+seed+hashtable.get_seed_size()+(j*l2_seed_size));
	hash = l2hashtable.l2_get_index(offset, i, j, pearson);
	l2_offset = l2hashtable.l2_get_offset2(hash);
	l2_frequency = l2hashtable.l2_get_frequency2(hash);
	for (uint32_t l = 0; l < l2_frequency; l++) {
	  location = l2hashtable.l2_get_location(l + l2_offset);
	  locations->insert(location-seed);
	}
      }
    }
  }
}

// Filters and finalizes the reads using the specified filters, SWA implementations
void filter_and_finalize_reads(string * read, set<uint32_t> * locations) {
  char reference[read_length];
  chrono::time_point<chrono::steady_clock> start, end;
  // Implement filters, based on flags set.
  for (set<uint32_t>::iterator it=locations->begin(); it != locations->end(); ++it) {
    // Get the reference DNA from the genome
    decompress_2bit_dna(reference, *it);
    // Filter
    start = chrono::steady_clock::now();
    bool pass_filter = filter_read_locations(read, reference);
    end = chrono::steady_clock::now();
    time_filter += (end - start).count();
    if (pass_filter && do_swa) {
      // SWA
      start = chrono::steady_clock::now();
      bool pass_swa = finalize_read_locations(read, reference);
      end = chrono::steady_clock::now();
      time_swa += (end - start).count();
      if (pass_swa) {
	cout << *it << " ";
      }
    }
  }
  cout << "\n";
}

// Returns false if read contains too many N's
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

// No SWA alignment
bool NoSWA(string * read, char * reference) {
  return true;
}

// Aligns using seqalign
bool SWA_Seqalign(string * read, char * reference) {
  return swaligner.sw_align(read->c_str(), reference, swa_threshold);
}

// Aligns using edlib
// This can be optimized further.
bool Meyers_Edlib(string * read, char * reference) {
  bool res = false;
  EdlibAlignResult result = edlibAlign(read->c_str(), read_length, reference, read_length, edlibNewAlignConfig(error_threshold, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
  if (result.editDistance != -1 ){
    res = true;
  }
  edlibFreeAlignResult(result);
  return res;
}

// Converts the read from ACTG to 0 1 2 3
void convert_for_opal(string read, unsigned char * destination, uint32_t read_length) {
  for (uint32_t i = 0; i < read_length; i++) {
    destination[i] = (unsigned char)char_values[read[i]];
  }
}

// Aligns using opal
bool Opal(string * read, char * reference) {
  // Need a function to convert read to 0,1,2,3
  // Need a function to convert the reference to 0,1,2,3 (unsigned char)
  unsigned char rd[read_length];
  unsigned char rf[read_length];
  unsigned char * ref[1];
  int reference_lengths[1];
  reference_lengths[0] = read_length;
  int reference_length = 1;

  convert_for_opal(*read, rd, read_length);
  convert_for_opal(reference, rf, read_length);
  ref[0] = rf;

  return opal_aligner.opal_swa(rd, read_length, ref, reference_length, reference_lengths);
}

bool SSW(string * read, char * reference) {
  StripedSmithWaterman::Aligner aligner;
  StripedSmithWaterman::Filter filter;
  StripedSmithWaterman::Alignment alignment;
  aligner.Align(read->c_str(), reference, read_length, filter, &alignment, read_length);
  if (alignment.sw_score >= 190)
   return true;
  return false;
}

bool NoFilter(string * read, char * reference) {
  return true;
}

bool SHD(string * read, char * reference) {
  return shd_filter.filter(read->c_str(), reference, error_threshold, read_length);
}

bool MAGNET(string * read, char * reference) {
  return shd_filter.magnet(read->c_str(), reference, error_threshold, read_length);
}

bool QGRAM(string * read, char * reference) {
  return true;
}

// Read the genome into 8 bit representation, using 0 1 2 3 instead of ACTG
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

// Read the genome into 2bit representation
void read_genome_2bit() {
  // size: 3.2 billion / 4 Bytes (approx 0xc0000000 >> 5 => 0x6000000)
  // WRONG! FIX ME.
  genome = (uint64_t *)malloc(0x7000000*sizeof(uint64_t));
  printf("Genome: %p-%p\n", genome, genome + 0x7000000);
  printf("Char_to_2bit: %p-%p\n", char_values, char_values+128);
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

// Decompresses the DNA back into ACTG characters
void decompress_2bit_dna(char * destination, uint32_t starting_index) {
  for (uint32_t i = 0; i < read_length; i++) {
    uint64_t dna = genome[(starting_index+i)>>5];
    uint8_t char_code = (dna >> (2*((starting_index + i)%32))) & 0x03;
    destination[i] = reverse_values[char_code];
  }
}

// Decompresses the DNA into characters containing 0, 1, 2, 3 (uint8, not char).
void decompress_2bit_dna_number(char * destination, uint32_t starting_index) {
  for (uint32_t i = 0; i < read_length; i++) {
    uint64_t dna = genome[(starting_index+i)>>5];
    uint8_t char_code = (dna >> (2*((starting_index + i)%32))) & 0x03;
    destination[i] = char_code;
  }
}

