TEST_CPP_FILES=genome_test.cpp hashtable_test.cpp main_test.cpp l2hashtable_test.cpp
DFLAGS=-DUSE_DEBUGGING_OUTPUT
CXX=g++
INC=..
SRC=..
CPP_SOURCES= $(SRC)/genome.cpp $(SRC)/hashtable.cpp $(SRC)/l2hashtable.cpp

mapper_tests: $(TEST_CPP_FILES)
	$(CXX) -std=c++11 -I $(INC) $(CPP_SOURCES) $^ -o $@ -lgtest -lpthread -g

clean:
	rm mapper_tests

run:
	sh ./run_tests.sh
