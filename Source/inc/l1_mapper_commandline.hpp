#ifndef L1_MAPPER_COMMANDLINE_H_
#define L1_MAPPER_COMMANDLINE_H_
#include <iostream>
#include <getopt.h>
#include <stdint.h>

using namespace std;
extern char* reads_filename;
extern char* directory_name;
extern char* hashtable_filename;
extern uint8_t filters;
extern uint8_t seed_selection_algorithm;
extern uint32_t number_of_seeds;
extern uint32_t error_threshold;
extern uint32_t limit;
extern uint32_t group;
std::string version = "1.0";
std::string default_directory = ".";

const uint8_t NAIVE = 0x1;
const uint8_t HOBBES = 0x2;
const uint8_t OPTIMAL = 0x4;

const uint8_t PAIRED_END = 0x1;
const uint8_t SHD = 0x2;
const uint8_t QGRAM = 0x4;

void print_options(){
  std::cout << "Usage: " << std::endl;
  std::cout << "-t  --hashtable\tPrefix for the desired L1 hashtable files" << std::endl;
  std::cout << "-r  --reads\tPrefix for the read file" << std::endl;
  std::cout << "-c  --seeds\tNumber of seeds to use." << std::endl; 
  std::cout << "\t\t\tIf the number is above the maximum, then the maximum will be used." << std::endl;
  std::cout << "-l  --limit\tThe limit on how much space the seeds can be moved." << std::endl; 
  std::cout << "\t\t\tIf the number is greater than the maximum, then the maximum will be used." << std::endl;
  std::cout << "-e  --error\tPercent error" << std::endl;
  std::cout << "-o  --optimal\tUse optimal seed selection" << std::endl;
  std::cout << "-b  --hobbes\tUse hobbes (sub-optimal) seed selection" << std::endl;
  std::cout << "-n  --naive\tUse naive seed selection" << std::endl;
  std::cout << "-p  --paired-end\tUse paired-end filtering" << std::endl;
  std::cout << "-s  --shd\tUse shifted-hamming-distance filtering" << std::endl;
  std::cout << "-q  --qgram\tUse Q-gram filtering" << std::endl;
  std::cout << "-g  --group\tProcess the reads in groups of n numbers." << std::endl;
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
    {"hashtable",   required_argument,  0,  'h'},
    {"reads",	    required_argument,  0,  'r'},
    {"seeds",	    required_argument,  0,  'c'},
    {"limit",	    required_argument,  0,  'l'},
    {"error",	    required_argument,  0,  'e'},
    {"optimal",	    required_argument,  0,  'o'},
    {"hobbes",	    required_argument,  0,  'b'},
    {"naive",	    required_argument,  0,  'n'},
    {"paired-end",  required_argument,	0,  'p'},
    {"shd",	    required_argument,  0,  's'},
    {"qgram",	    required_argument,  0,  'q'},
    {"group",	    required_argument,  0,  'g'},
    {"directory",   required_argument,  0,  'd'},
    {"help",	    no_argument,	0,  'h'},
    {"version",	    no_argument,	0,  'v'},
    {0,0,0,0}
  };

  // Write to the current directory unless otherwise directed.
  directory_name = (char*)default_directory.c_str();
  group = 1;

  // Check that only one seed selection algorithm is selected.
  // Load all options into variables.
  // Defaults: naive seed selection, no filters, write to std out.
  while( (c = getopt_long(argc, argv, "r:t:s:e:l:c:g:obnpsqd:vh", longOptions, &index)) != -1){
    switch(c)
    {
      case 't':
	hashtable_filename = optarg;
	required--;
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
      case 'g':
	group = atoi(optarg);
	break;
      case 'e':
	error_threshold = atoi(optarg);
	break;
      case 'n':
	seed_selection_algorithm |= NAIVE;
	break;
      case 'b':
	seed_selection_algorithm |= HOBBES;
	break;
      case 'o':
	seed_selection_algorithm |= OPTIMAL;
	break;
      case 'p':
	filters |= PAIRED_END;
	break;
      case 's':
	filters |= SHD;
	break;
      case 'q':
	filters |= QGRAM;
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

  if (seed_selection_algorithm == 0x3 || seed_selection_algorithm > 0x4) {
    std::cerr << "Please only select one seed selection algorithm." << std::endl;
    print_options();
    exit(1);
  }

  if (!seed_selection_algorithm)
    seed_selection_algorithm |= NAIVE;

  return 0;

}

#endif //L1_MAPPER_COMMANDLINE_H_
