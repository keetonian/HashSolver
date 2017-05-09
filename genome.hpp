#ifndef GENOME_H_
#define GENOME_H_

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include "stdint.h"

class Genome {
  private:
    std::vector<std::vector<char> * > genome_vector;
    //vector<string* > genome_vector;
    //string ** genome;

    std::vector<std::string> genome_names;
#if defined (USE_DEBUGGING_OUTPUT)
    bool debug = true;
#else
    bool debug = false;
#endif

  public:
    Genome();
    Genome(char* filename);
    Genome(char* filename, int options);
    ~Genome();


    char* getstring(uint32_t index);
    char* getstring(std::string chromosome, uint32_t index);

    void change(uint32_t index, char c);
    void change(std::string chromosome, uint32_t index, char c);
    int change_all_N(uint32_t consecutive = 10);

    void read_genome(const char* filename, int options);
    void write_genome(const char* filename, int options);

    std::vector<char> * get_chromosome(std::string chromosome);
    std::vector<std::vector<char> * > * get_genome();
    std::vector<std::string> * get_genome_names();

    void delete_genome();

  protected:
    int error_check();
    char random_base();
    char change_base(char cc);

    char* getstring(uint32_t chromosome_index, uint32_t index);
    uint32_t get_chromosome_index(std::string chromosome);
};

#endif //GENOME_H_
