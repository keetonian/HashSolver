CXX=g++
CFLAGS=-std=c++11 -O2 -Wall -Wno-unused-variable #-pthread -lpthread

all: build query query2 check print_seed

build: construct_table.cpp hashtable.o l2hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^

hashtable.o: hashtable.h hashtable.c
	gcc -O2 -c $^

l2hashtable.o: l2hashtable.h l2hashtable.c
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
	rm build hashtable.o query

run:
	./build ../human_g1k_v37_hs37d5.fasta 2 14 16

run2: 
	./build ../20.fa 1 14 16
