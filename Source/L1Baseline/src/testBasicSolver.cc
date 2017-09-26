#include "basicSolver.h"
#include <vector>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cstring>
#include <zlib.h>
#include <errno.h>
#include <assert.h>
#include <map>
#include "kseq.h"
#include "hashtable.h"


using namespace std;

int main(int argc, const char* argv[]) {
	int frequency;
	map <int, unsigned long long> frequencyCounter;

	BasicSolver solver;

	if (argc != 5) {
		cerr << "USAGE: ./testOptimalSolver <BWA index base> <Read File> <number of seeds> <seed length>" << endl;
		exit(1);
	}

	gzFile fp;
	kseq_t *ks;

	string name = argv[1];
	solver.loadHashTables(name);

	

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
	  solver.init(ks->seq.l, seedNum, seedLength);
	  frequency = solver.solveDNA(ks->seq.s);
	  frequencyCounter[frequency]++;
	  if(frequency > 1000000)
	    cout << ks->seq.s << '\t' << frequency << endl;
	}

	
	char benchFileName[80];
	strcpy(benchFileName, argv[2]);
	char* benchName = strtok(benchFileName, ".");
	
	ofstream output( (string(benchName) + "." + to_string(seedNum) + "-" + to_string(seedLength) + string(".basic")) , ofstream::out);

	cout << "seedNum: " << seedNum << " | seedLength: " << seedLength << endl;
	//output << "seedNum: " << seedNum << " | seedLength: " << seedLength << endl;
	for (map<int, unsigned long long>::iterator iter = frequencyCounter.begin(); iter != frequencyCounter.end(); iter++)
		output << iter->first << " " << iter->second << endl;

	
	//output.close();

	kseq_destroy(ks);
	gzclose(fp);
	
	return 0;
}
