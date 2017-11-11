#ifndef L1_MAPPER_H_
#define L1_MAPPER_H_

#include <stdint.h>
#include "l1_hashtable.hpp"
#include "seed_solver.hpp"
#include "sw_aligner.hpp"
#include "shd_filter.hpp"
#include "opal_aligner.hpp"
#include "mapper_common.hpp"
#include "bwt.h"

typedef bool (*Finalize_Reads) (std::string * read, char * reference);
typedef bool (*Filter_Reads) (std::string * read, char * reference);
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
double time_seeds, time_locations, time_filter, time_swa;

bool fasthash_seed_selection = false;
bool optimal_seed_selection = false;

std::vector<uint32_t> *locations;
std::vector<uint32_t> *reverse_locations;

SWAFunction swa_function;
SeedSelection seed_selection;
FilterAlgorithm filter_algorithm;

SeedSolver * solver;
Hashtable hashtable;
SWAligner swaligner;
SHDFilter shd_filter;
OpalAligner opal_aligner;

Finalize_Reads finalize_read_locations = NULL;
Filter_Reads filter_read_locations = NULL;
Decompress_DNA decompress_dna = NULL;

void map_read(std::string read);

void get_seeds(std::string * read, uint8_t * seeds);

void get_locations(std::string * read, uint8_t * seeds, std::vector<uint32_t> & locations);

void filter_and_finalize_reads(std::string * reads, std::vector<uint32_t> & locations);

bool pre_filter(string read);

std::string reverse_read(std::string read);

bool NoSWA(std::string * read, char * reference);
bool SWA_Seqalign(std::string * read, char * reference);
bool Meyers_Edlib(std::string * read, char * reference);
bool Opal(std::string * read, char * reference);
bool SSW(std::string * read, char * reference);

bool NoFilter(std::string * read, char * reference);
bool SHD(std::string * read, char * reference);
bool MAGNET(std::string * read, char * reference);
bool QGRAM(std::string * read, char * reference);

void decompress_2bit_dna(char * destination, uint32_t starting_index);
void decompress_2bit_dna_pac(char * destination, uint32_t starting_index);
void convert_for_opal(std::string read, unsigned char * destination, uint32_t read_length);

void read_genome_2bit();

void read_genome_char();

bwt_t * load_bwt(const char *hint);

void load_bwt_genome(const char *filename, uint64_t size);

#endif //L1_MAPPER_H_
