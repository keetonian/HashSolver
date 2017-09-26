#! /bin/bash

# Exit if any command fails.
set -e

# Make sure files and directories are where they should be.
# Put everything in /tmp for speed, save SSD life.

TEST_DIR=/tmp/DNA/testing

mkdir -p -v $TEST_DIR

if [ ! -f $TEST_DIR/test_genome.fa ]; then
  echo "Copying test genome to $TEST_DIR"
  cp ./test_files/test_genome.fa $TEST_DIR/
fi

# Run tests
./mapper_tests

# Remove testing files
rm -rf /tmp/DNA/testing
