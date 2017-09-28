#! /bin/bash

# Exit if any command fails.
set -e

# Make sure files and directories are where they should be.
# Put everything in /tmp for speed, save SSD life.

TEST_DIR=/tmp/DNA/testing

mkdir -p -v $TEST_DIR

if [ ! -f $TEST_DIR/test_genome.fa ]; then
  echo "Copying test genome to $TEST_DIR"
  find . -name "test_genome.fa" -exec cp {} $TEST_DIR/ \;
fi

# Run tests
$(find . -name "mapper_tests")

echo $?

# Remove testing files
rm -rf /tmp/DNA/testing
