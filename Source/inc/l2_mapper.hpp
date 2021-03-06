#ifndef L2_MAPPER_H_
#define L2_MAPPER_H_

#include <stdint.h>
#include <boost/heap/fibonacci_heap.hpp>
#include "l1_hashtable.hpp"
#include "l2_hashtable.hpp"
#include "seed_solver.hpp"
#include "sw_aligner.hpp"
#include "shd_filter.hpp"
#include "opal_aligner.hpp"
#include "mapper_common.hpp"
#include "bwt.h"
#include "fast_hash_solver.hpp"

typedef bool (*Finalize_Reads) (const std::string &read, char * reference);
typedef int (*Filter_Reads) (const std::string &read, char * reference);
typedef void (*Decompress_DNA) (char * destination, uint32_t starting_index);

char* reads_filename = 0;
char* directory_name = 0;
char* hashtable_filename = 0;
char* bwt_filename = 0;
uint8_t filters;
uint8_t seed_selection_algorithm;
uint32_t number_of_seeds = 1;
uint32_t error_threshold;
uint32_t swa_threshold;
uint32_t limit;
uint64_t * genome;
uint8_t * bwa_genome;
unsigned char * genome_char;
uint32_t read_length;
bool do_swa;
bool do_filter;
bool do_filter2;
double time_seeds, time_locations, time_filter, time_swa, time_decompress, total_shift_time;

bool fasthash_seed_selection = false;
bool optimal_seed_selection = false;
int fasthash_seeds;

SWAFunction swa_function;
SeedSelection seed_selection;
FilterAlgorithm filter_algorithm;
FilterAlgorithm second_filter_algorithm;

std::vector<uint32_t> locations;
std::vector<uint32_t> reverse_locations;

SeedSolver * solver;
Hashtable hashtable;
L2Hashtable l2hashtable;
SWAligner swaligner;
SHDFilter shd_filter;
OpalAligner opal_aligner;
FastHashSolver fasthash;

Finalize_Reads finalize_read_locations = NULL;
Filter_Reads filter_read_locations = NULL;
Filter_Reads filter_read_locations2 = NULL;
Decompress_DNA decompress_dna = NULL;

void map_read(const std::string &read);

void get_seeds(const std::string &read, uint8_t * seeds);

void get_locations(const std::string & read, uint8_t * seeds, std::vector<uint32_t> & locations);

void filter_and_finalize_reads(const std::string &read, std::vector<uint32_t> & locations);

bool pre_filter(const string &read);

std::string reverse_read(const std::string &read);

bool NoSWA(const std::string &read, char * reference);
bool SWA_Seqalign(const std::string &read, char * reference);
bool Meyers_Edlib(const std::string &read, char * reference);
bool Opal(const std::string &read, char * reference);
bool SSW(const std::string &read, char * reference);

int NoFilter(const std::string &read, char * reference);
int SHD(const std::string &read, char * reference);
int MAGNET(const std::string &read, char * reference);
int QGRAM(const std::string &read, char * reference);

void decompress_2bit_dna(char * destination, uint32_t starting_index);
void decompress_2bit_dna_pac(char * destination, uint32_t starting_index);
void convert_for_opal(const char* read, unsigned char * destination, uint32_t read_length);

void read_genome_2bit();

void read_genome_char();

bwt_t * load_bwt(const char *hint);

void load_bwt_genome(const char *filename, uint64_t size);

#endif //L2_MAPPER_H_
