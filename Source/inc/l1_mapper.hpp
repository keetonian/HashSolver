#ifndef L1_MAPPER_H_
#define L1_MAPPER_H_

#include <stdint.h>
#include "l1_hashtable.hpp"
#include "seed_solver.hpp"
#include "sw_aligner.hpp"
#include "shd_filter.hpp"
#include "opal_aligner.hpp"
#include "mapper_common.hpp"

typedef struct {
  uint32_t frequency;
  uint8_t * seeds;
  std::string * read;
  uint32_t * locations;
} ReadInformation;

typedef void (*Finalize_Reads) (ReadInformation * reads, uint32_t number_of_reads);

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

SWAFunction swa_function;
SeedSelection seed_selection;

SeedSolver * solver;
Hashtable hashtable;
SWAligner swaligner;
SHDFilter shd_filter;
OpalAligner opal_aligner;

Finalize_Reads finalize_read_locations = NULL;

void get_seeds(ReadInformation * reads, uint32_t number_of_reads);

void get_locations(ReadInformation * reads, uint32_t number_of_reads);

void free_read_memory(ReadInformation * reads, uint32_t number_of_reads);

void filter_reads(ReadInformation * reads, uint32_t number_of_reads);

bool pre_filter(string read);

void NoSWA(ReadInformation * reads, uint32_t number_of_reads);
void SWA_Seqalign(ReadInformation * reads, uint32_t number_of_reads);
void Meyers_Edlib(ReadInformation * reads, uint32_t number_of_reads);
void Opal(ReadInformation * reads, uint32_t number_of_reads);

void decompress_2bit_dna(char * destination, uint32_t starting_index);
void convert_read(std::string read, unsigned char * destination, uint32_t read_length);

void read_genome_2bit();

void read_genome_char();


#endif //L1_MAPPER_H_
