CXX=g++
CFLAGS=-std=c++11 -O2 -Wall  #-Wno-unused-variable -I . #-pthread -lpthread
DFLAGS=-DUSE_DEBUGGING_OUTPUT -DUSE_MULTITHREADING
TBB=./tbb/build/linux_intel64_gcc_cc5.4.0_libc2.23_kernel4.4.0_release
LIB=-ltbb -lpthread -pthread -lrt

all: build query query2 check print_seed

build: construct_table.cpp hashtable.o l2hashtable.o genome.o #concurrent_vector.o
	$(CXX) $(CFLAGS) $(DFLAGS) -o $@ $^ -L$(TBB) -I./tbb/include -I . $(LIB)

hashtable.o: hashtable.c
	gcc -O2 -c $^

genome.o: genome.cpp
	$(CXX) $(CFLAGS) $(DFLAGS) -c $^ -I .

l2hashtable.o: l2hashtable.c
	gcc -O2 -c $^

query: query_hash.cpp hashtable.o l2hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^

query2: query_hash_l2.cpp hashtable.o l2hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^

check: check_hash.cpp hashtable.o l2hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^

print_seed: seed_finder.cpp hashtable.o l2hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^

frequencies: l2_frequencies.cpp l2hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^

clean:
	rm build hashtable.o query genome.o l2hashtable.o

run:
	./build -g /tmp/GRCh38_full_analysis_set_plus_decoy_hla.fa -s 14 -c 2

