#include "hobbesSolver.h"
#include <vector>
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


using namespace std;

int main(int argc, const char* argv[]) {
  int frequency;
  map <int, unsigned long long> frequencyCounter;

  HobbesSolver solver;

  if (argc != 5) {
    cerr << "USAGE: ./testOptimalSolver <hash 2-level base> <Read File> <number of seeds> <seed length>" << endl;
    exit(1);
  }

  gzFile fp;
  kseq_t *ks;

  string tablename = argv[1];
  tablename += ".hashtable";
  string locname = argv[1];
  locname += ".locations";
  string l2tablename = argv[1];
  l2tablename += ".l2.hashtable";
  string l2locname = argv[1];
  l2locname += ".l2.locations";

  read_table_from_file(tablename.c_str());
  cout << "Table read" << endl;
  read_locations_from_file(locname.c_str());
  cout << "Locations read" << endl;

  l2_read_hashtable_from_file(l2tablename.c_str());
  l2_read_locations_from_file(l2locname.c_str());

  fp = strcmp(argv[2], "-")? gzopen(argv[2], "r") : gzdopen(fileno(stdin), "r");
  if (NULL == fp) {
    fprintf(stderr, "Couldn't open %s : %s\n",
	strcmp(argv[2], "-") ? argv[2] : "stdin",
	errno ? strerror(errno) : "Out of memory");
    exit(EXIT_FAILURE);
  }
  ks = kseq_init(fp); // initialize the FASTA/Q parser

  int seedNum, seedLength;
  seedNum = atoi(argv[3]);
  seedLength = atoi(argv[4]);


  while (kseq_read(ks) >= 0) {
    //cout << ks->seq.s << endl;
    solver.init(ks->seq.l, seedNum, seedLength);
    frequency = solver.solveDNA(ks->seq.s);
    //cout << frequency << endl;
    if(frequency < 0)
      continue;
    frequencyCounter[frequency]++;
    if(frequency > 1000000)
      cout << ks->seq.s << '\t' << frequency << endl;
  }


  cout << "Done" << endl;
  //int wait;
  //cin >> wait;
  char benchFileName[80];
  strcpy(benchFileName, argv[2]);
  char* benchName = strtok(benchFileName, ".");

  //ofstream output( (string(benchName) + "." + to_string(seedNum) + "-" + to_string(seedLength) + string(".hobbes")) , ofstream::out);

  cout << "seedNum: " << seedNum << " | seedLength: " << seedLength << endl;
  //output << "seedNum: " << seedNum << " | seedLength: " << seedLength << endl;
  for (map<int, unsigned long long>::iterator iter = frequencyCounter.begin(); iter != frequencyCounter.end(); iter++)
    cout << iter->first << " " << iter->second << endl;


  //output.close();

  kseq_destroy(ks);
  gzclose(fp);

  return 0;
}
