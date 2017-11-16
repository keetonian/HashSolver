#include <vector>
#include <cstdlib>
#include <algorithm>
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
#include "l1_mapper.hpp"
#include "l1_hobbes_solver.hpp"
#include "l1_basic_solver.hpp"
#include "fastHASHSolver.h"
#include "seed_solver.hpp"
#include "sw_aligner.hpp"
#include "shd_filter.hpp"
#include <chrono>
#include "edlib.h"
#include "mapper_common.hpp"
#include "opal_aligner.hpp"
#include "ssw_cpp.h"
#include "optimalSolverLN.h"
#include "fast_hash_solver.hpp"
#include "bwt.h"


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

  // Load BWT objects
  bwt_t * bwt = load_bwt(bwt_filename);
  if (NULL == bwt) {
    cerr << "bwt failed to load" << endl;
    return 1;
  }
  cerr << bwt << endl;

  // Load genome
  read_genome_2bit();
  //read_genome_char();
  //cout << "genome read" << endl;

  do_filter = true;

  uint32_t reserve_space = 5000000;
  locations.reserve(reserve_space);
  reverse_locations.reserve(reserve_space);
  locations.push_back(0);
  reverse_locations.push_back(0);
  fprintf(stderr, "Locations:%p-%p\n", &(locations.at(0)), &(locations.at(0)) + reserve_space);
  fprintf(stderr, "Locations2:%p-%p\n", &(reverse_locations.at(0)), &(reverse_locations.at(0)) + reserve_space);

  decompress_dna = & decompress_2bit_dna;

  // Select Solver
  switch(seed_selection) {
    case SeedSelection::naive: 
      solver = new BasicSolver();
      solver->loadTables(&hashtable);
      break;
    case SeedSelection::hobbes: 
      solver = new HobbesSolver();
      solver->loadTables(&hashtable);
      break;
    case SeedSelection::optimal:
      solver = new OptimalSolverLN();
      solver->loadTables(bwt);
      optimal_seed_selection = true;
      decompress_dna = &decompress_2bit_dna_pac;
      break;
    case SeedSelection::fasthash:
      solver = new FastHASHSolver();
      solver->loadTables(&hashtable);
      fasthash_seed_selection = true;
      break;
  }

  // Select SWA implementation
  switch(swa_function) {
    case SWAFunction::noswa: 
      finalize_read_locations = &NoSWA;
      break;
    case SWAFunction::seqalign: 
      finalize_read_locations = &SWA_Seqalign;
      break;
    case SWAFunction::edlib: 
      finalize_read_locations = &Meyers_Edlib;
      break;
    case SWAFunction::opal: 
      finalize_read_locations = &Opal;
      break;
    case SWAFunction::ssw: 
      finalize_read_locations = &SSW;
      break;
  }

  switch(filter_algorithm) {
    case FilterAlgorithm::none: 
      filter_read_locations = &NoFilter;
      do_filter = false;
      break;
    case FilterAlgorithm::SHD: 
      filter_read_locations = &SHD;
      break;
    case FilterAlgorithm::MAGNET: 
      filter_read_locations = &MAGNET;
      break;
    case FilterAlgorithm::QGRAM: 
      filter_read_locations = &QGRAM;
      break;
  }
  switch(second_filter_algorithm) {
    case FilterAlgorithm::none: 
      filter_read_locations2 = &NoFilter;
      do_filter2 = false;
      break;
    case FilterAlgorithm::SHD: 
      filter_read_locations2 = &SHD;
      break;
    case FilterAlgorithm::MAGNET: 
      filter_read_locations2 = &MAGNET;
      break;
    case FilterAlgorithm::QGRAM: 
      filter_read_locations2 = &QGRAM;
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
   *	ks->seq.l is the length of the reads.
   *	Check:
   *	  limit: max value is len(read) - number_of_seeds * seed_size
   *	  number_of_seeds: max value is len(read) / seed_size
   *	Don't make everything automatically assume 100 bp reads.
   */

  // Find read length
  read_length = ks->seq.l;

  // Init seed solver
  solver->init(read_length, number_of_seeds+fasthash_seeds, hashtable.get_seed_size(), limit);
  fasthash.init(read_length, number_of_seeds, hashtable.get_seed_size(), number_of_seeds + fasthash_seeds, seed_selection == SeedSelection::optimal ? 1 : 0);
  fasthash.load_hashtable(&hashtable);

  // Prepare Smith Waterman step
  swa_threshold = read_length - error_threshold;

  // Initialize the read information for simple batching.

  time_seeds=0; time_locations=0; time_filter=0; time_swa=0;

  cout << "start" << endl;
  // PIPELINE
  do {
    //cout << ks->name.s << "\n";
    // Pre filter out erroneous reads
    if (!pre_filter(ks->seq.s))
      continue;

    // Map the read
    map_read(ks->seq.s);

  } while (kseq_read(ks) >= 0);

  cout << "end" << endl;

  double period = 1000000000;
  double total_time = time_seeds + time_locations + time_filter + time_swa;
  cerr << "Timing: " << endl;
  cerr << "Seeds: " << time_seeds/period << "\t" << time_seeds / total_time << endl;
  //cerr << "Locations: " << time_locations/period << endl;
  cerr << "Filters: " << time_filter/period << "\t" << time_filter / total_time << endl;
  cerr << "SWA: " << time_swa/period << " " << time_swa / total_time << endl;
  cerr << "Total" << total_time/period << endl;

  // Free memory that was used.
  kseq_destroy(ks);
  gzclose(fp);
  free(genome);
  //free(genome_char);
  hashtable.free_memory();

  return 0;
}

// Maps a single read.
void map_read(const string &read) {


  // Get the reverse read to map as well.
  string reverse = reverse_read(read);
  uint8_t seeds[number_of_seeds + fasthash_seeds];
  uint8_t reverse_seeds[number_of_seeds + fasthash_seeds];

  auto start = chrono::steady_clock::now();
  // seed selection
  if (!fasthash_seeds && !optimal_seed_selection) {
    solver->solveDNA(read, seeds);
    solver->solveDNA(reverse, reverse_seeds);
    locations.clear();
    reverse_locations.clear();

    // location selection
    // Keep locations on the stack unless I have to move them.
    get_locations(read, seeds, locations);
    get_locations(reverse, reverse_seeds, reverse_locations);
  }
  else if (optimal_seed_selection && !fasthash_seed_selection) {
    locations.clear();
    reverse_locations.clear();
    solver->solveDNA(read, seeds, locations);
    solver->solveDNA(reverse, reverse_seeds, reverse_locations);
  }
  else if (optimal_seed_selection && fasthash_seeds) {
    solver->solveDNA(read, seeds, locations);
    solver->solveDNA(reverse, reverse_seeds, reverse_locations);
    locations.clear();
    reverse_locations.clear();
    
    fasthash.solveDNA(seeds, locations);
    fasthash.solveDNA(reverse_seeds, reverse_locations);
  }
  else if (!optimal_seed_selection && fasthash_seeds) {
    solver->solveDNA(read, seeds);
    solver->solveDNA(reverse, reverse_seeds);
    locations.clear();
    reverse_locations.clear();

    // Fast hash get locations
    fasthash.solveDNA(seeds, locations);
    fasthash.solveDNA(reverse_seeds, reverse_locations);
  }

  // Sort locations, remove duplicates
  std::sort(locations.begin(), locations.end());
  locations.erase(std::unique(locations.begin(), locations.end()), locations.end());
  std::sort(reverse_locations.begin(), reverse_locations.end());
  reverse_locations.erase(std::unique(reverse_locations.begin(), reverse_locations.end()), reverse_locations.end());
  // filtering and SWA

  auto end = chrono::steady_clock::now();
  time_seeds += (end - start).count();

  //cout << read << "\n";
  filter_and_finalize_reads(read, locations);
  //cout << "Reverse:" << "\n";
  filter_and_finalize_reads(reverse, reverse_locations);
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
void get_seeds(const string &read, uint8_t * seeds) {
  solver->solveDNA(read, seeds);
}

// Uses the seeds found using the seed selection step
// Gets all of the locations, puts into a fibonacci heap(no duplicates)
void get_locations(const string &read, uint8_t * seeds, std::vector<uint32_t> &locations) {
  // Set up variables
  uint32_t seed;
  uint32_t hash;
  uint32_t offset;
  uint32_t frequency;
  uint32_t location;

  // Loop through every seed
  for(uint32_t j = 0; j < number_of_seeds; j++) {
    seed = seeds[j];
    hash = hashtable.get_hash(&(read[seed]));
    offset = hashtable.get_offset(hash);
    frequency = hashtable.get_frequency(hash);
    //cout << frequency << endl;

    // Find all locations for a seed
    for (uint32_t k = 0; k < frequency; k++) {
      location = hashtable.get_location(offset+k);
      // Subtract seed position to make sure locations start at beginning of read
      locations.push_back(location - seed);
    }
  }
}

// Filters and finalizes the reads using the specified filters, SWA implementations
void filter_and_finalize_reads(const string &read, std::vector<uint32_t> & locations) {
  char reference[read_length];
  chrono::time_point<chrono::steady_clock> start, end;
  // Implement filters, based on flags set.
  for (auto it=locations.begin(); it != locations.end(); ++it) {
    // Get the reference DNA from the genome
    decompress_dna(reference, *it);
    // Filter
    start = chrono::steady_clock::now();
    int pass_filter = do_filter ? 0 : 1;
    if (do_filter)
      pass_filter = filter_read_locations(read, reference);
    if (do_filter2 && pass_filter == 1)
      pass_filter = filter_read_locations2(read, reference);
    end = chrono::steady_clock::now();
    time_filter += (end - start).count();
    if (pass_filter == 1 && do_swa) {
      // SWA
      start = chrono::steady_clock::now();
      bool pass_swa = finalize_read_locations(read, reference);
      end = chrono::steady_clock::now();
      time_swa += (end - start).count();
      if (pass_swa) {
	//cout << *it << " ";
      }
    }
    else if(pass_filter == 2 && do_swa) {
      //cout << *it << " ";
    }
  }
  //cout << "\n";
}

// Returns false if read contains too many N's
bool pre_filter(string read) {
  uint32_t count = 0;
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
bool NoSWA(const string &read, char * reference) {
  return true;
}

// Aligns using seqalign
bool SWA_Seqalign(const string &read, char * reference) {
  return swaligner.sw_align(&(read[0]), reference, swa_threshold);
}

// Aligns using edlib
// This can be optimized further.
bool Meyers_Edlib(const string &read, char * reference) {
  bool res = false;
  EdlibAlignResult result = edlibAlign(&(read[0]), read_length, reference, read_length, edlibNewAlignConfig(error_threshold, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
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
bool Opal(const string &read, char * reference) {
  // Need a function to convert read to 0,1,2,3
  // Need a function to convert the reference to 0,1,2,3 (unsigned char)
  unsigned char rd[read_length];
  unsigned char rf[read_length];
  unsigned char * ref[1];
  int reference_lengths[1];
  reference_lengths[0] = read_length;
  int reference_length = 1;

  convert_for_opal(read, rd, read_length);
  convert_for_opal(reference, rf, read_length);
  ref[0] = rf;

  return opal_aligner.opal_swa(rd, read_length, ref, reference_length, reference_lengths);
}

bool SSW(const string &read, char * reference) {
  StripedSmithWaterman::Aligner aligner(1, 1, 1, 1);
  StripedSmithWaterman::Filter filter(false, false, read_length-error_threshold, read_length + error_threshold);
  StripedSmithWaterman::Alignment alignment;
  aligner.Align(&(read[0]), reference, read_length, filter, &alignment, read_length);
  if (alignment.sw_score >= 95)
    return true;
  return false;
}

int NoFilter(const string &read, char * reference) {
  return 1;
}

int SHD(const string &read, char * reference) {
  return shd_filter.filter(&(read[0]), reference, error_threshold, read_length);
}

int MAGNET(const string &read, char * reference) {
  return shd_filter.magnet(&(read[0]), reference, error_threshold, read_length);
}

int QGRAM(const string &read, char * reference) {
  return 1;
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
  size_t genome_size = 0x7000000;
  genome = (uint64_t *)malloc(genome_size*sizeof(uint64_t));
  fprintf(stderr, "Genome:%p-%p\n", genome, genome + genome_size);
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
  //cout << index << endl;
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

void decompress_2bit_dna_pac(char * destination, uint32_t starting_index) {
  for(uint32_t k = 0; k < read_length; k++) {
    destination[k] = reverse_values[bwa_genome[((starting_index + k)>>2)] >>((~(k+starting_index)&3)<<1) & 3];
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

bwt_t *load_bwt(const char *hint)
{
  string tmp;
  string prefix;
  bwt_t *bwt;
  prefix = hint;
  if (prefix == "") {
    return 0;
  }

  // Load BWT
  tmp = prefix + ".bwt";
  bwt = bwt_restore_bwt(tmp.c_str());

  // Load SA
  tmp = prefix + ".sa";
  bwt_restore_sa(tmp.c_str(), bwt);

  // Load size of genome (from ANN)
  tmp = prefix + ".ann";
  ifstream ann(tmp.c_str());
  uint64_t genome_size;
  if (ann.is_open()) {
    // Genome size should be first entry in ann file
    ann >> genome_size;
  } else {
    cerr << "Unable to open bwt genome file (.pac)." << endl;
  }
  ann.close();

  // Load BWT Genome
  tmp = prefix + ".pac";
  load_bwt_genome(tmp.c_str(), genome_size); 
  return bwt;
}

void load_bwt_genome(const char * filename, uint64_t size) {
  // Get the right number of bytes (size/4 + 3)
  size = (size>>2) + 3;
  bwa_genome = (uint8_t *)malloc(size * sizeof(uint8_t));
  fprintf(stderr, "BWTGenome:%p-%p\n", bwa_genome, bwa_genome + size);
  FILE * f = fopen(filename, "rb");
  size_t data = fread(bwa_genome,  sizeof(uint8_t), size, f);
  if (data < size - 4) {
    cerr << "Expected " << size << " but read " << data << " from " << filename << endl;
  }
  fclose(f);
}
