#ifndef L2_COMMANDLINE_H_
#define L2_COMMANDLINE_H_

#include <iostream>
#include <getopt.h>

using namespace std;
extern char* genome_file;
extern char* directory;
extern size_t seed_size;
extern uint32_t replace_n;
std::string version = "1.0";
extern uint32_t l2thresh;
std::string default_directory = ".";

void print_options(){
  std::cout << "Usage: " << std::endl;
  std::cout << "-g  --genome\tPath to the reference genome" << std::endl;
  std::cout << "-c  --contigs\tSpecify how many contigs to use." << std::endl;
  std::cout << "\t\t  0: Use all" << std::endl;
  std::cout << "\t\t  24: Use main contigs" << std::endl;
  std::cout << "-s  --seedsize\tSpecify seed size." << std::endl;
  std::cout << "\t\t  Sizes between 6 and 16 are preferred." << std::endl;
  std::cout << "-N  --replaceN\tReplace consecutive N's of specified length with random chars." << std::endl;
  std::cout << "-l  --l2thresh\tL2 bucket threshold." << std::endl;
  std::cout << "\t\t  0: Don't use L2 hashing"<< std::endl;
  std::cout << "\t\t  Recommended to use sizes greater than 1023." << std::endl;
  std::cout << "-f  --hash\tUse specified hash function." << std::endl;
  std::cout << "\t\t  Options:" << std::endl;
  std::cout << "\t\t    ICompress1" << endl;
  std::cout << "\t\t    ICompress2" << endl;
  std::cout << "\t\t    JCompress1" << endl;
  std::cout << "\t\t    JCompress2" << endl;
  std::cout << "\t\t    IJCompress1" << endl;
  std::cout << "\t\t    IJCompress2" << endl;
  std::cout << "\t\t    Hobbes" << endl;
  std::cout << "\t\t    LSH1" << endl;
  std::cout << "\t\t    LSH2" << endl;
  std::cout << "\t\t    Freq" << endl;
  std::cout << "\t\t    Mask" << endl;
  std::cout << "-d  --directory\tOutput to the specified directory." << std::endl;
  std::cout << "-w  --write_out\tDisable writing the modified genome to file" << std::endl;
  std::cout << "-h  --help\tPrint this message." << std::endl;
  std::cout << "-v  --version\tPrint the software version" << std::endl;
}

int parseCommands(int argc, char** argv){
  int index;
  int c;
  int required = 2;

  static struct option longOptions[] =
  {
    {"genome",	    required_argument,  0,  'g'},
    {"contigs",	    required_argument,  0,  'c'},
    {"seedsize",    required_argument,  0,  's'},
    {"replaceN",    required_argument,  0,  'N'},
    {"l2thresh",    required_argument,	0,  'l'},
    {"hash",	    required_argument,	0,  'f'},
    {"directory",   required_argument,	0,  'd'},
    {"write_out",   required_argument,	0,  'w'},
    {"help",	    no_argument,	0,  'h'},
    {"version",	    no_argument,	0,  'v'},
    {0,0,0,0}
  };

  // Write to the current directory unless otherwise directed.
  directory = (char*)default_directory.c_str();

  while( (c = getopt_long(argc, argv, "g:d:c:s:N:l:f:vh", longOptions, &index)) != -1){
    switch(c)
    {
      case 'g':
	genome_file = optarg;
	required--;
	break;
      case 'd':
	directory = optarg;
	break;
      case 'c':
	contigs = atoi(optarg);
	break;
      case 's':
	seed = atoi(optarg);
	required--;
	break;
      case 'N':
	replace_n = atoi(optarg);
	break;
      case 'l':
	l2_threshold= atoi(optarg);
	break;
      case 'w':
	write_genome_to_file = false;
	break;
      case 'h':
	print_options();
	exit(0);
      case 'v':
	std::cout << "Version " << version << std::endl;
	exit(0);
    }
  }

  return 0;

}

#endif //L2_COMMANDLINE_H_
