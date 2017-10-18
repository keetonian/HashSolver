CXX=g++
CFLAGS=-std=c++14 -O2 -Wall -g -Wextra -Wno-char-subscripts
DFLAGS=-DUSE_DEBUGGING_OUTPUT -DUSE_MULTITHREADING
LIBRARIES=./lib
TBB=$(LIBRARIES)/tbb/build/linux_intel64_gcc_cc6.3.0_libc2.24_kernel4.10.0_release
LIB=-ltbb -lpthread -pthread -lrt
INC=./inc
SRC=./src
OBJ=./build
TST=./tests
BIN=./bin
SWA=./lib/seq-align
SWA_ITEMS=-I $(SWA)/libs -I $(SWA)/src -L $(SWA)/libs/string_buffer -L $(SWA)/src -lalign -lstrbuf -lz
SHD=./lib/Shifted-Hamming-Distance
SHD_ITEMS=-I $(SHD) -mbmi -msse4.2
SHD_OBJECTS=$(SHD)/mask.o $(SHD)/print.o $(SHD)/bit_convert.o $(SHD)/popcount.o $(SHD)/vector_filter.o
EDLIB=./lib/edlib
EDLIB_ITEMS=-I $(EDLIB)/edlib/include -L $(EDLIB)/build/lib -L $(EDLIB) -ledlib_static
OPAL=./lib/opal
OPAL_ITEMS=-I $(OPAL)/src
SSW=./lib/Complete-Striped-Smith-Waterman-Library
SSW_ITEMS=-I $(SSW)/src
SSW_FILES=$(SSW)/src/ssw_cpp.cpp $(SSW)/src/ssw.c
TEST_CPP_FILES=$(TST)/genome_test.cpp $(TST)/l1_hashtable_test.cpp $(TST)/main_test.cpp $(TST)/l2_hashtable_test.cpp

L1Baseline=$(BIN)/L1BaselineConstructor $(BIN)/L1BaselineQuery $(BIN)/L1BaselineCheck $(BIN)/L1BaselinePrintSeed $(BIN)/L1BaselineMapper
L1Compressed=$(BIN)/L1CompressedConstructor $(BIN)/L1CompressedQuery $(BIN)/L1CompressedCheck $(BIN)/L1CompressedPrintSeed #$(BIN)/L1CompressedMapper
L2Baseline=$(BIN)/L2BaselineConstructor $(BIN)/L2BaselineQuery $(BIN)/L2BaselineCheck $(BIN)/L2BaselinePrintSeed
L2CompressedBaseline=$(BIN)/L2CompressedBaselineConstructor $(BIN)/L2CompressedBaselineQuery $(BIN)/L2CompressedBaselineCheck $(BIN)/L2CompressedBaselinePrintSeed

EXE=$(L1Baseline) $(L1Compressed) $(L2Baseline) $(L2CompressedBaseline)

all: $(EXE)

$(OBJ)/l1_hashtable.o: $(SRC)/l1_hashtable.cpp
	$(CXX) $(CFLAGS) $(DFLAGS) -c $^ -o $@ -I $(INC)

$(OBJ)/l1_hashtable_compressed.o: $(SRC)/l1_hashtable_compressed.cpp
	$(CXX) $(CFLAGS) $(DFLAGS) -c $^ -o $@ -I $(INC)

$(OBJ)/genome.o: $(SRC)/genome.cpp
	$(CXX) $(CFLAGS) $(DFLAGS) -c $^ -o $@ -I $(INC)

$(OBJ)/l2_hashtable.o: $(SRC)/l2_hashtable.cpp
	$(CXX) $(CFLAGS) $(DFLAGS) -c $^ -o $@ -I $(INC)

$(OBJ)/l1_hobbes_solver.o: $(SRC)/l1_hobbes_solver.cpp $(OBJ)/l1_hashtable.o
	$(CXX) $(CFLAGS) -c $^ -o $@ -I $(INC)

$(OBJ)/l1_basic_solver.o: $(SRC)/l1_basic_solver.cpp $(OBJ)/l1_hashtable.o
	$(CXX) $(CFLAGS) -c $^ -o $@ -I $(INC)

$(OBJ)/utils.o: $(SRC)/utils.c 
	gcc -O2 -Wall -g -DHAVE_PTHREAD -DUSE_MALLOC_WRAPPERS -c $^ -o $@ -I $(INC) #-fPIC

$(OBJ)/malloc_wrap.o: $(SRC)/malloc_wrap.c
	gcc -O2 -Wall -g -DHAVE_PTHREAD -DUSE_MALLOC_WRAPPERS -c $^ -o $@ -I $(INC) #-fPIC

$(OBJ)/l2_hobbes_solver.o: $(SRC)/l2_hobbes_solver.cc $(OBJ)/l1_hashtable.o $(OBJ)/l2_hashtable.o
	$(CXX) $(CFLAGS) -c $^ -o $@ -I $(INC)

$(OBJ)/sw_aligner.o: $(SRC)/sw_aligner.cpp 
	$(CXX) $(CFLAGS) -c $^ -o $@ -I $(INC) $(SWA_ITEMS)

$(OBJ)/shd_filter.o: $(SRC)/shd_filter.cpp $(SHD_OBJECTS)
	$(CXX) $(CFLAGS) -c $^ -o $@ -I $(INC) $(SHD_ITEMS)

$(OBJ)/opal_aligner.o: $(SRC)/opal_aligner.cpp $(OPAL)/src/opal.o
	$(CXX) $(CFLAGS) -c $^ -o $@ -I $(INC) $(OPAL_ITEMS)



$(BIN)/L1BaselineConstructor: $(SRC)/l1_construct_table.cpp  $(OBJ)/genome.o $(OBJ)/l1_hashtable.o
	$(CXX) $(CFLAGS) $(DFLAGS) -o $@ $^ -L$(TBB) -I $(LIBRARIES)/tbb/include -I $(INC) $(LIB)

$(BIN)/L1BaselineMapper: $(SRC)/l1_mapper.cpp  $(OBJ)/genome.o $(OBJ)/l1_hashtable.o $(OBJ)/l1_hobbes_solver.o $(OBJ)/utils.o $(OBJ)/malloc_wrap.o $(OBJ)/l1_basic_solver.o $(OBJ)/sw_aligner.o $(OBJ)/shd_filter.o $(SHD_OBJECTS) $(OBJ)/opal_aligner.o $(OPAL)/src/opal.o $(SSW_FILES)
	$(CXX) $(CFLAGS) $(DFLAGS) -o $@ $^ -I $(INC) $(SWA_ITEMS) $(SHD_ITEMS) $(EDLIB_ITEMS) $(OPAL_ITEMS) $(SSW_ITEMS) -lz -lm -lrt -lpthread

$(BIN)/L1BaselineQuery: $(SRC)/l1_query_hash.cpp $(OBJ)/l1_hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^ -I $(INC)

$(BIN)/L1BaselineCheck: $(SRC)/l1_check_hash.cpp $(OBJ)/l1_hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^ -I $(INC)

$(BIN)/L1BaselinePrintSeed: $(SRC)/l1_seed_finder.cpp $(OBJ)/l1_hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^ -I $(INC)


$(BIN)/L1CompressedConstructor: $(SRC)/l1_construct_table.cpp  $(OBJ)/genome.o $(OBJ)/l1_hashtable_compressed.o
	$(CXX) $(CFLAGS) $(DFLAGS) -o $@ $^ -L$(TBB) -I $(LIBRARIES)/tbb/include -I $(INC) $(LIB)

$(BIN)/L1CompressedMapper: $(SRC)/l1_mapper.cpp  $(OBJ)/genome.o $(OBJ)/l1_hashtable_compressed.o $(OBJ)/l1_hobbes_solver.o $(OBJ)/utils.o $(OBJ)/malloc_wrap.o $(OBJ)/l1_basic_solver.o $(OBJ)/sw_aligner.o $(OBJ)/shd_filter.o $(SHD_OBJECTS)
	$(CXX) $(CFLAGS) $(DFLAGS) -o $@ $^ -I $(INC) $(SWA_ITEMS) $(SHD_ITEMS) $(EDLIB_ITEMS) -lz -lm -lrt -lpthread

$(BIN)/L1CompressedQuery: $(SRC)/l1_query_hash.cpp $(OBJ)/l1_hashtable_compressed.o
	$(CXX) $(CFLAGS) -o $@ $^ -I $(INC)

$(BIN)/L1CompressedCheck: $(SRC)/l1_check_hash.cpp $(OBJ)/l1_hashtable_compressed.o
	$(CXX) $(CFLAGS) -o $@ $^ -I $(INC)

$(BIN)/L1CompressedPrintSeed: $(SRC)/l1_seed_finder.cpp $(OBJ)/l1_hashtable_compressed.o
	$(CXX) $(CFLAGS) -o $@ $^ -I $(INC)



$(BIN)/L2BaselineConstructor: $(SRC)/l2_construct_table.cpp $(OBJ)/l1_hashtable.o $(OBJ)/l2_hashtable.o $(OBJ)/genome.o
	$(CXX) $(CFLAGS) $(DFLAGS) -o $@ $^ -L$(TBB) -I $(LIBRARIES)/tbb/include -I $(INC) $(LIB)

$(BIN)/L2BaselineQuery: $(SRC)/l2_query_hash.cpp $(OBJ)/l1_hashtable.o $(OBJ)/l2_hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^ -I $(INC)

$(BIN)/L2BaselineCheck: $(SRC)/l2_check_hash.cpp $(OBJ)/l1_hashtable.o $(OBJ)/l2_hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^ -I $(INC)

$(BIN)/L2BaselinePrintSeed: $(SRC)/l2_seed_finder.cpp $(OBJ)/l1_hashtable.o $(OBJ)/l2_hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^ -I $(INC)


$(BIN)/L2CompressedBaselineConstructor: $(SRC)/l2_construct_table.cpp $(OBJ)/l1_hashtable_compressed.o $(OBJ)/l2_hashtable.o $(OBJ)/genome.o
	$(CXX) $(CFLAGS) $(DFLAGS) -o $@ $^ -L$(TBB) -I $(LIBRARIES)/tbb/include -I $(INC) $(LIB)

$(BIN)/L2CompressedBaselineQuery: $(SRC)/l2_query_hash.cpp $(OBJ)/l1_hashtable_compressed.o $(OBJ)/l2_hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^ -I $(INC)

$(BIN)/L2CompressedBaselineCheck: $(SRC)/l2_check_hash.cpp $(OBJ)/l1_hashtable_compressed.o $(OBJ)/l2_hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^ -I $(INC)

$(BIN)/L2CompressedBaselinePrintSeed: $(SRC)/l2_seed_finder.cpp $(OBJ)/l1_hashtable_compressed.o $(OBJ)/l2_hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^ -I $(INC)


$(BIN)/L2BaselineFrequencies: $(SRC)/l2_frequencies.cpp $(OBJ)/l2_hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^ -I $(INC)

$(TST)/mapper_tests: $(TEST_CPP_FILES) $(OBJ)/genome.o $(OBJ)/l1_hashtable.o $(OBJ)/l2_hashtable.o
	$(CXX) $(CFLAGS) -o $@ $^ -I $(INC) -lgtest -lpthread -pthread

clean:
	rm $(EXE) $(OBJ)/*.o

test:
	sh $(TST)/run_tests.sh


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

#hobbesSolver.o: hobbesSolver.cc hobbesSolver.h l1_hashtable.hpp hashtable.o l2_hashtable.o
#	$(CXX) $(CPPFLAGS) -c $< -o $@ 

#basicSolver.o: basicSolver.cc basicSolver.h
#	$(CXX) $(CPPFLAGS) $(INC) -c $< -o $@ 

#testHobbesSolver: testHobbesSolver.cc hobbesSolver.o l1_hashtable.o l2_hashtable.o utils.o malloc_wrap.o
#	$(CXX) $(CPPFLAGS) -o $@ $^ $(LIB)

#testBasicSolver: testBasicSolver.cc basicSolver.o l1_hashtable.o utils.o malloc_wrap.o
#	$(CXX) $(CPPFLAGS) $(INC) -o $@ $^ $(LIB) 
