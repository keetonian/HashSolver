#ifndef L1_MAPPER_COMMANDLINE_H_
#define L1_MAPPER_COMMANDLINE_H_
#include <iostream>
#include <getopt.h>
#include <stdint.h>
#include "mapper_common.hpp"

using namespace std;
extern char* reads_filename;
extern char* directory_name;
extern char* hashtable_filename;
extern char* bwt_filename;
extern uint8_t filters;
extern SeedSelection seed_selection;
extern SWAFunction swa_function;
extern FilterAlgorithm filter_algorithm;
extern FilterAlgorithm second_filter_algorithm;
extern uint32_t number_of_seeds;
extern uint32_t error_threshold;
extern uint32_t limit;
extern bool do_swa;
extern int fasthash_seeds;
std::string version = "1.0";
std::string default_directory = ".";

const uint8_t PAIRED_END = 0x1;

void print_options(){
  std::cout << "Usage: " << std::endl;
  std::cout << "-t  --hashtable\tPrefix for the desired L1 hashtable files" << std::endl;
  std::cout << "-b  --bwt\tPrefix for the desired bwt files" << std::endl;
  std::cout << "-r  --reads\tName of read file" << std::endl;
  std::cout << "-c  --seeds\tNumber of seeds to use." << std::endl; 
  std::cout << "\t\t\tIf the number is above the maximum, then the maximum will be used." << std::endl;
  std::cout << "-l  --limit\tThe limit on how much space the seeds can be moved." << std::endl; 
  std::cout << "\t\t\tIf the number is greater than the maximum, then the maximum will be used." << std::endl;
  std::cout << "-e  --error\tPercent error" << std::endl;
  std::cout << "-s  --seed-select\tChoose seed selection algorithm" << std::endl;
  std::cout << "\t\t\t0: Naive." << std::endl;
  std::cout << "\t\t\t1: Hobbes." << std::endl;
  std::cout << "\t\t\t2: Optimal." << std::endl;
  std::cout << "\t\t\t2: Naive + FastHash." << std::endl;
  std::cout << "-p  --paired\tUse paired-end alignment (NOT IMPLEMENTED)" << std::endl;
  std::cout << "-x  --swa-option\tChoose what to do with the swa step." << std::endl;
  std::cout << "\t\t\t0: don't do SWA and don't print matches." << std::endl;
  std::cout << "\t\t\t1: pass all possible locations as matches." << std::endl;
  std::cout << "\t\t\t2: use seq-align's SWA implementation" << std::endl;
  std::cout << "\t\t\t3: use edlib's SWA implementation" << std::endl;
  std::cout << "\t\t\t4: use opal's SWA implementation" << std::endl;
  std::cout << "\t\t\t5: use Complete-Striped-Smith-Waterman-Library's SWA implementation" << std::endl;
  std::cout << "-f  --filter\tSelect which additional filtering to use:" << std::endl;
  std::cout << "\t\t\t0: No additional filtering" << std::endl;
  std::cout << "\t\t\t1: Use SHD filtering" << std::endl;
  std::cout << "\t\t\t2: Use MAGNET filtering" << std::endl;
  std::cout << "\t\t\t3: Use Q-gram filtering (NOT IMPLEMENTED)" << std::endl;
  std::cout << "-F  --fasthash\t Use fast hash filter:" << std::endl;
  std::cout << "\t\t\t#: # of extra seeds to use in fasthash filter" << std::endl;
  std::cout << "-d  --directory\tOutput to the specified directory." << std::endl;
  std::cout << "-h  --help\tPrint this message." << std::endl;
  std::cout << "-v  --version\tPrint the software version" << std::endl;
}

int parseCommands(int argc, char** argv){
  int index;
  int c;
  int required = 2; // genome and hashtable

  static struct option longOptions[] =
  {
    {"hashtable",   required_argument,  0,  't'},
    {"bwt",	    required_argument,  0,  'b'},
    {"reads",	    required_argument,  0,  'r'},
    {"seeds",	    required_argument,  0,  'c'},
    {"limit",	    required_argument,  0,  'l'},
    {"error",	    required_argument,  0,  'e'},
    {"seed-select", required_argument,  0,  's'},
    {"paired",	    required_argument,	0,  'p'},
    {"swa-option",  required_argument,  0,  'x'},
    {"filter",	    required_argument,  0,  'f'},
    {"fasthash",    required_argument,  0,  'F'},
    {"directory",   required_argument,  0,  'd'},
    {"help",	    no_argument,	0,  'h'},
    {"version",	    no_argument,	0,  'v'},
    {0,0,0,0}
  };

  // Write to the current directory unless otherwise directed.
  directory_name = (char*)default_directory.c_str();
  do_swa = true;
  seed_selection = SeedSelection::hobbes;
  swa_function = SWAFunction::edlib;
  fasthash_seeds = 0;
  second_filter_algorithm = FilterAlgorithm::none;
  FilterAlgorithm *filt = & filter_algorithm;

  // Check that only one seed selection algorithm is selected.
  // Load all options into variables.
  // Defaults: naive seed selection, no filters, write to std out.
  while( (c = getopt_long(argc, argv, "r:t:s:e:l:c:x:f:F:b:psqd:vh", longOptions, &index)) != -1){
    switch(c)
    {
      case 't':
	hashtable_filename = optarg;
	required--;
	break;
      case 'b':
	bwt_filename = optarg;
	break;
      case 'r':
	reads_filename = optarg;
	required--;
	break;
      case 'd':
	directory_name = optarg;
	break;
      case 'c':
	number_of_seeds = atoi(optarg);
	break;
      case 'l':
	limit = atoi(optarg);
	break;
      case 'e':
	error_threshold = atoi(optarg);
	break;
      case 's':
	switch(atoi(optarg)){
	  case 0: seed_selection = SeedSelection::naive;
		  break;
	  case 1: seed_selection = SeedSelection::hobbes;
		  break;
	  case 2: seed_selection = SeedSelection::optimal;
		  break;
	  case 3: seed_selection = SeedSelection::fasthash;
		  break;
	}
	break;
      case 'p':
	filters |= PAIRED_END;
	break;
      case 'x':
	switch(atoi(optarg)) {
	  case 0: swa_function = SWAFunction::noswa;
		  do_swa = false;
		  break;
	  case 1: swa_function = SWAFunction::noswa;
		  break;
	  case 2: swa_function = SWAFunction::seqalign;
		  break;
	  case 3: swa_function = SWAFunction::edlib;
		  break;
	  case 4: swa_function = SWAFunction::opal;
		  break;
	  case 5: swa_function = SWAFunction::ssw;
		  break;
	}
	break;
      case 'f':
	switch(atoi(optarg)) {
	  case 0: *filt = FilterAlgorithm::none;
		  break;
	  case 1: *filt = FilterAlgorithm::SHD;
		  break;
	  case 2: *filt = FilterAlgorithm::MAGNET;
		  break;
	  case 3: *filt = FilterAlgorithm::QGRAM;
		  break;
	}
	filt = &second_filter_algorithm;
	break;
      case 'F':
	fasthash_seeds = atoi(optarg);
	break;
      case 'h':
	print_options();
	exit(0);
      case 'v':
	std::cout << "Version " << version << std::endl;
	exit(0);
    }
  }

  if (required) {
    std::cerr << "Please provide a hashtable prefix and a file of reads." << std::endl;
    print_options();
    exit(1);
  }

  return 0;

}

#endif //L1_MAPPER_COMMANDLINE_H_
