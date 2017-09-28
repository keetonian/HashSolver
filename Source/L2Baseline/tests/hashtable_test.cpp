#include "testing_common.hpp"
#include "hashtable.hpp"
#include <iostream>

TEST(hashtable_test, creates_empty_hashtable) {
  Hashtable hashtable = Hashtable();
  EXPECT_EQ(hashtable.get_seed_size(), 14);
  EXPECT_EQ(hashtable.get_table_size(), 4ULL << (2 * (14-1)));
}

TEST(hashtable_test, creates_empty_hashtable_for_size) {
  Hashtable hashtable = Hashtable(16);
  EXPECT_EQ(hashtable.get_seed_size(), 16);
  EXPECT_EQ(hashtable.get_table_size(), 4ULL << (2 * (16-1)));
}

// Make a test for something too small or too big?

TEST(hashtable_test, set_seed_size) {
  Hashtable hashtable = Hashtable();
  hashtable.set_seed_size(13);
  EXPECT_EQ(hashtable.get_seed_size(), 13);
  EXPECT_EQ(hashtable.get_table_size(), 4ULL << (2 * (13-1)));
  hashtable.set_seed_size(10);
  EXPECT_EQ(hashtable.get_seed_size(), 10);
  EXPECT_EQ(hashtable.get_table_size(), 4ULL << (2 * (10-1)));
}

TEST(hashtable_test, initialize_table_has_0_entries) {
  Hashtable hashtable = Hashtable(10);
  hashtable.initialize_hashtable();
  for (uint64_t i = 0; i < hashtable.get_table_size(); i++) {
    EXPECT_EQ(hashtable.get_frequency(i), 0);
    EXPECT_EQ(hashtable.get_offset(i), 0);
  }
}

TEST(hashtable_test, get_hash) {
  Hashtable hashtable = Hashtable(10);
  ASSERT_EQ(hashtable.get_hash("AAAAAAAAAA"), 0);
  ASSERT_EQ(hashtable.get_hash("TTTTTTTTTT"), hashtable.get_table_size()-1);
}

uint32_t populate_hashtable(Hashtable * hashtable) {
  uint32_t offset = 0;
  for(uint64_t i = 0; i < hashtable->get_table_size(); i++) {
    uint32_t frequency = rand() % 100;
    hashtable->set_frequency(i, frequency);
    hashtable->set_offset(i, offset);
    offset += frequency;
  }
  return offset;
}

TEST(hashtable_test, add_table_values) {
  Hashtable hashtable = Hashtable(10);
  hashtable.initialize_hashtable();
  populate_hashtable(&hashtable);
  uint64_t offset = 0;
  for(uint64_t i = 0; i < hashtable.get_table_size(); i++) {
    EXPECT_EQ(hashtable.get_offset(i), offset);
    offset += hashtable.get_frequency(i);
  }
}

TEST(hashtable_test, store_and_load_hash_table) {
  Hashtable hashtable = Hashtable(10);
  hashtable.initialize_hashtable();
  populate_hashtable(&hashtable);
  std::string temp_file_name = temporary_file_directory + "/hash_table";
  hashtable.write_table_to_file(temp_file_name.c_str());
  Hashtable table_from_file = Hashtable();
  table_from_file.read_table_from_file(temp_file_name.c_str());
  EXPECT_EQ(hashtable.get_seed_size(), table_from_file.get_seed_size());
  EXPECT_EQ(hashtable.get_table_size(), table_from_file.get_table_size());
  for (uint32_t i = 0; i < hashtable.get_table_size(); i++) {
    EXPECT_EQ(hashtable.get_frequency(i), table_from_file.get_frequency(i));
    EXPECT_EQ(hashtable.get_offset(i), table_from_file.get_offset(i));
  }
}

/**
 * Sets up the hashtable, locations array
 * Populates the locations array with the frequency of the seed.
 */
void populate_locations(Hashtable* hashtable) {
  uint64_t locations_size = populate_hashtable(hashtable);
  hashtable->initialize_location(locations_size);
  for (uint64_t i = 0; i < hashtable->get_table_size(); i++) {
    for (uint64_t j = 0; j < hashtable->get_frequency(i); j++) {
      hashtable->set_location(hashtable->get_offset(i) + j, hashtable->get_frequency(i));
    }
  }
}

TEST(hashtable_test, add_location_values) {
  Hashtable hashtable = Hashtable(10);
  hashtable.initialize_hashtable();
  populate_locations(&hashtable);
  for (uint64_t i = 0; i < hashtable.get_table_size(); i++) {
    for (uint64_t j = 0; j < hashtable.get_frequency(i); j++) {
      EXPECT_EQ(hashtable.get_location(hashtable.get_offset(i) + j), hashtable.get_frequency(i));
    }
  }
}

TEST(hashtable_test, store_and_load_locations) {
  Hashtable hashtable = Hashtable(10);
  hashtable.initialize_hashtable();
  populate_locations(&hashtable);
  std::string temp_file_name = temporary_file_directory + "hash_table";
  hashtable.write_to_file(temp_file_name.c_str());
  Hashtable table_from_file = Hashtable(temp_file_name);
  EXPECT_EQ(hashtable.get_locations_size(), table_from_file.get_locations_size());
  for (uint64_t i = 0; i < hashtable.get_locations_size(); i++) {
    EXPECT_EQ(hashtable.get_location(i), table_from_file.get_location(i));
  }
}


