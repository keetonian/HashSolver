#include "testing_common.hpp"
#include "genome.hpp"
#include <iostream>

TEST(genome_test, creates_empty_genome) {
  Genome genome;
  EXPECT_EQ(genome.get_genome()->size(), 0);
}

TEST(genome_test, loads_genome_no_options) {
  Genome genome = Genome(test_genome_file_location.c_str());
  EXPECT_TRUE(genome.get_genome()->size() > 0);
  EXPECT_EQ(genome.get_genome()->size(), genome.get_genome_names()->size());
}

TEST(genome_test, stores_genome) {
  Genome genome = Genome(test_genome_file_location.c_str());
  std::string genome_output_file = temporary_file_directory + "/generated_genome.fa";
  genome.write_genome(genome_output_file.c_str(), 0);
  Genome written_genome = Genome(genome_output_file.c_str());
  EXPECT_EQ(genome.get_genome()->size(), written_genome.get_genome()->size());
  for(uint32_t i = 0; i < genome.get_genome()->size(); i++) {
    EXPECT_EQ(genome.get_genome()->at(i)->size(), written_genome.get_genome()->at(i)->size());
    for(uint32_t j = 0; j < genome.get_genome()->at(i)->size(); j++) {
      EXPECT_EQ(genome.get_genome()->at(i)->at(j), written_genome.get_genome()->at(i)->at(j));
    }
  }
}

bool removedConsecutiveNs(std::vector<std::vector<char> *> * chromosomes, uint32_t number_consecutive) {
  for(uint32_t i = 0; i < chromosomes->size(); i++) {
    std::vector<char>* chromosome = chromosomes->at(i);
    uint32_t counter = 0;
    for(uint32_t j = 0; j < chromosome->size(); j++) {
      if (chromosome->at(j) == 'N') {
	counter++;
      }
      else if(counter > 0) {
	if (counter <= number_consecutive) {
	  std::cout << "Failed on chromosome: " << i 
	    << "\nDid not change " << number_consecutive << " Ns.\n" << "Printing previous characters." << std::endl;
	  int k = (j - number_consecutive - 1) < 0 ? 0 : j - number_consecutive - 1;
	  for (; k <= j; k++ ) {
	    std::cout << chromosome->at(k);
	  }
	  std::cout << std::endl;
	  return false;
	}
	counter = 0;
      }
    }
  }
  return true;
}


TEST(genome_test, remove_entries_of_N) {
  Genome genome = Genome(test_genome_file_location.c_str());
  for(uint32_t i = 0; i < 10; i++) {
    genome.change_all_N(i);
    EXPECT_TRUE(removedConsecutiveNs(genome.get_genome(), i));
  }
}
