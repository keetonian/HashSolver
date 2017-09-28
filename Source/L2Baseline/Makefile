CXX=g++
CFLAGS=-std=c++17 -O2 -Wall -g
DFLAGS=-DUSE_DEBUGGING_OUTPUT -DUSE_MULTITHREADING
LIBRARIES=../lib
TBB=$(LIBRARIES)/tbb/build/linux_intel64_gcc_cc5.4.0_libc2.23_kernel4.4.0_release
LIB=-ltbb -lpthread -pthread -lrt
INC=./inc
SRC=./src
OBJ=./build
TST=./tests
BIN=./bin
OBJECTS= $(OBJ)/genome.o $(OBJ)/hashtable.o $(OBJ)/l2hashtable.o
TEST_CPP_FILES=$(TST)/genome_test.cpp $(TST)/hashtable_test.cpp $(TST)/main_test.cpp $(TST)/l2hashtable_test.cpp

EXE= $(BIN)/L2ConstructTable $(BIN)/L2Query $(BIN)/L2Check $(BIN)/L2PrintSeed $(TST)/mapper_tests #testHobbesSolver testBasicSolver

all: $(EXE)

$(BIN)/L2ConstructTable: $(SRC)/construct_table.cpp $(OBJ)/hashtable.o $(OBJ)/l2hashtable.o $(OBJ)/genome.o
	$(CXX) $(CFLAGS) $(DFLAGS) -o $@ $^ -L$(TBB) -I $(LIBRARIES)/tbb/include -I $(INC) $(LIB)

$(OBJ)/hashtable.o: $(SRC)/hashtable.cpp
	$(CXX) $(CFLAGS) $(DFLAGS) -c $^ -I $(INC)
	mv hashtable.o $(OBJ)/

$(OBJ)/genome.o: $(SRC)/genome.cpp
	$(CXX) $(CFLAGS) $(DFLAGS) -c $^ -I $(INC)
	mv genome.o $(OBJ)/

$(OBJ)/l2hashtable.o: $(SRC)/l2hashtable.cpp
	$(CXX) $(CFLAGS) $(DFLAGS) -c $^ -I $(INC)
	mv l2hashtable.o $(OBJ)/

$(BIN)/L2Query: $(SRC)/query_hash_l2.cpp $(OBJ)/hashtable.o $(OBJ)/l2hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^ -I $(INC)

$(BIN)/L2Check: $(SRC)/check_hash.cpp $(OBJ)/hashtable.o $(OBJ)/l2hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^ -I $(INC)

$(BIN)/L2PrintSeed: $(SRC)/seed_finder.cpp $(OBJ)/hashtable.o $(OBJ)/l2hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^ -I $(INC)

$(BIN)/L2Frequencies: $(SRC)/l2_frequencies.cpp $(OBJ)/l2hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^ -I $(INC)

$(TST)/mapper_tests: $(TEST_CPP_FILES)
	$(CXX) $(CFLAGS) -o $@ $^ -I $(INC) $(OBJECTS) -lgtest -lpthread

clean:
	rm $(EXE) $(OBJ)/*.o

test:
	sh $(TST)/run_tests.sh

run:
	./build -g /tmp/DNA/GRCh38_full_analysis_set_plus_decoy_hla.fa -s 14 -c 2

run2:
	./build -g /tmp/DNA/GRCh38_full_analysis_set_plus_decoy_hla.fa -s 14 -c 2 -l 2047

runall:
	./build -g /tmp/DNA/GRCh38_full_analysis_set_plus_decoy_hla.fa -s 14 -c 0


#XX = g++
##CPPFLAGS=-O3
##CPPFLAGS=-g -DDEBUG
#CPPFLAGS=-O3 -DTEST -DHAVE_PTHREAD -std=c++11 -g
#LIB= -lz -lrt -lm -lpthread
#EXE= testHobbesSolver testBasicSolver 

#all: $(EXE)

#KmerHash.o: KmerHash.cc KmerHash.h
#	$(CXX) $(CPPFLAGS) -c $<
	
#RefDB.o: RefDB.cc RefDB.h
#	$(CXX) $(CPPFLAGS) -c $<

#hobbesSolver.o: hobbesSolver.cc hobbesSolver.h hashtable.hpp hashtable.o l2hashtable.o
#	$(CXX) $(CPPFLAGS) -c $< -o $@ 

#basicSolver.o: basicSolver.cc basicSolver.h
#	$(CXX) $(CPPFLAGS) $(INC) -c $< -o $@ 

#testHobbesSolver: testHobbesSolver.cc hobbesSolver.o hashtable.o l2hashtable.o utils.o malloc_wrap.o
#	$(CXX) $(CPPFLAGS) -o $@ $^ $(LIB)

#testBasicSolver: testBasicSolver.cc basicSolver.o hashtable.o utils.o malloc_wrap.o
#	$(CXX) $(CPPFLAGS) $(INC) -o $@ $^ $(LIB) 
