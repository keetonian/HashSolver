CXX=g++
CFLAGS=-std=c++11 -O2 -Wall -Wno-unused-variable -pthread -lpthread

all: build query

build: construct_table_thread.cpp hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^

hashtable.o: hashtable.h hashtable.c
	gcc -O2 -c $^

query: query_hash.cpp hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^

clean:
	rm build hashtable.o query

run:
	./build ../human_g1k_v37_hs37d5.fasta 2 14 16

run2: 
	./build ../20.fa 1 14 16
