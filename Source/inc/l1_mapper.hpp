#ifndef L1_MAPPER_H_
#define L1_MAPPER_H_

#include <stdint.h>
#include "l1_hashtable.hpp"
#include "seed_solver.hpp"

typedef struct {
  uint32_t frequency;
  uint8_t * seeds;
  std::string * read;
  uint32_t * locations;
} ReadInformation;

char* reads_filename = 0;
char* directory_name = 0;
char* hashtable_filename = 0;
uint8_t filters;
uint8_t seed_selection_algorithm;
uint32_t number_of_seeds = 1;
uint32_t error_threshold;
uint32_t limit;
uint32_t group;
uint64_t * genome;
uint32_t read_length;

SeedSolver * solver;
Hashtable hashtable;

void get_seeds(ReadInformation * reads, uint32_t number_of_reads);

void get_locations(ReadInformation * reads, uint32_t number_of_reads);

void free_read_memory(ReadInformation * reads, uint32_t number_of_reads);

void filter_reads(ReadInformation * reads, uint32_t number_of_reads);

void finalize_read_locations(ReadInformation * reads, uint32_t number_of_reads);

void decompress_2bit_dna(char * destination, uint32_t starting_index);

void read_genome();


#endif //L1_MAPPER_H_
