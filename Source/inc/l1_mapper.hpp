#ifndef L1_MAPPER_H_
#define L1_MAPPER_H_

#include <stdint.h>
#include <set>
#include "l1_hashtable.hpp"
#include "seed_solver.hpp"
#include "sw_aligner.hpp"
#include "shd_filter.hpp"
#include "opal_aligner.hpp"
#include "mapper_common.hpp"

typedef bool (*Finalize_Reads) (std::string * read, char * reference);
typedef bool (*Filter_Reads) (std::string * read, char * reference);

char* reads_filename = 0;
char* directory_name = 0;
char* hashtable_filename = 0;
uint8_t filters;
uint8_t seed_selection_algorithm;
uint32_t number_of_seeds = 1;
uint32_t error_threshold;
uint32_t swa_threshold;
uint32_t limit;
uint32_t group;
uint64_t * genome;
unsigned char * genome_char;
uint32_t read_length;
bool do_swa;
double time_seeds, time_locations, time_filter, time_swa;

SWAFunction swa_function;
SeedSelection seed_selection;

SeedSolver * solver;
Hashtable hashtable;
SWAligner swaligner;
SHDFilter shd_filter;
OpalAligner opal_aligner;

Finalize_Reads finalize_read_locations = NULL;
Filter_Reads filter_read_locations = NULL;

void map_read(std::string read);

void get_seeds(std::string * read, uint8_t * seeds);

void get_locations(std::string * read, uint8_t * seeds, std::set<uint32_t> * locations);

void filter_and_finalize_reads(std::string * reads, std::set<uint32_t> * locations);

bool pre_filter(string read);

std::string reverse_read(std::string read);

bool NoSWA(std::string * read, char * reference);
bool SWA_Seqalign(std::string * read, char * reference);
bool Meyers_Edlib(std::string * read, char * reference);
bool Opal(std::string * read, char * reference);

void decompress_2bit_dna(char * destination, uint32_t starting_index);
void convert_for_opal(std::string read, unsigned char * destination, uint32_t read_length);

void read_genome_2bit();

void read_genome_char();


#endif //L1_MAPPER_H_
