#include "shd_filter.hpp"
#include "vector_filter.h"
#include <string.h>

SHDFilter::SHDFilter() {
  for (int i = 0; i < 128; i++) {
    blank_string[i] = 0x0;
  }
  strncpy(shd_read, blank_string, 128);
  strncpy(shd_ref, blank_string, 128);
}

void SHDFilter::init(uint32_t error_threshold, uint32_t read_length) {
  hamming_masks = (char*)malloc(sizeof(char) * (error_threshold*2+1) * read_length);
}

int SHDFilter::magnet(const char * read, const char * ref, int error, int read_length) {
  strncpy(shd_read, read, read_length);
  strncpy(shd_ref, ref, read_length);
  return bit_magnet_filter_sse1(shd_read, shd_ref, read_length, error);
}

int SHDFilter::filter(const char * read, const char * ref, int error, int read_length) {
  strncpy(shd_read, read, read_length);
  strncpy(shd_ref, ref, read_length);
  return bit_vec_filter_sse1(shd_read, shd_ref, read_length, error);
}

int SHDFilter::magnet_scalar(const char * read, const char* ref, int error, int read_length) {
  int edit_distance = 0;
  int i = 0;
  int magnet_result = 0;


  // Create normal hamming mask
  for (; i < read_length; i++) {
    char c = !(read[i] == ref[i]);
    hamming_masks[i] = c;
    edit_distance += c;
  }
  if (edit_distance <= error) {
    return 2;
  }

  uint32_t offset = 0;
  // Create shifted hamming masks
  for (int j = 1; j <= error; j++) {
    // Shift right
    edit_distance = 0;
    offset += read_length;

    // Shift read right
    for(i = 0; i < read_length; i++) {
      if (i < j)
	hamming_masks[i+offset] = 0;
      else {
	char c = !(read[i] == ref[j+i]);
	hamming_masks[i+offset] = c;
	edit_distance += c;
      }
    }

    if (edit_distance <= error) {
      return 2;
    }

    offset += read_length;
    edit_distance = 0;
    // Shift read left
    for (i = 0; i < read_length; i++) {
      if(i >= read_length - j)
	hamming_masks[i+offset] = 0;
      else {
	char c = !(read[i+j] == ref[i]);
	hamming_masks[i + offset] = c;
	edit_distance += c;
      }
    }

    if (edit_distance <= error) {
      return 2;
    }
  }

  // TODO do i repeat this error or error+1 times?
  for (int k = 0; k <= error+1; k++) {
    uint32_t longest_subsequence = 0;
    uint32_t starting_index = 0;
    offset = 0;
    // Find longest subsequence of 0's
    for (i = 0; i < 2*error + 1; i++) {
      uint32_t subsequence = 0;
      for (int j = 0; j < read_length; j++) {
	if(!hamming_masks[offset+j])
	  subsequence++;
	else {
	  if (subsequence > longest_subsequence) {
	    longest_subsequence = subsequence;
	    starting_index = j - subsequence;
	  }
	  subsequence = 0;
	  if (longest_subsequence > read_length - j)
	    break;
	}
      }
      offset += read_length;
    }

    // Update result
    // This can be a scalar, since we know exactly what we're adding
    // It needs to add up to read_length - error.
    magnet_result += longest_subsequence;
    if (magnet_result > read_length - error)
      return 1;

    // Update other vectors
    offset = 0;
    for (i = 0; i < 2*error+1; i++) {
      uint32_t j = starting_index ? starting_index - 1 : 0;
      for (; j < starting_index + longest_subsequence + 1; j++) {
	hamming_masks[j + offset] = 1;
      }
      offset += read_length;
    }

  }

  return 0;
}


