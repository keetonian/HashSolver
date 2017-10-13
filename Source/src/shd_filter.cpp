#include "shd_filter.hpp"
#include "vector_filter.h"
#include <string.h>

SHDFilter::SHDFilter() {
  strncpy(shd_read, blank_string, 128);
  strncpy(shd_ref, blank_string, 128);
}

bool SHDFilter::filter(const char * read, const char * ref, int error, int read_length) {
  strncpy(shd_read, read, read_length);
  strncpy(shd_ref, ref, read_length);
  if(bit_vec_filter_sse1(shd_read, shd_ref, read_length, error))
    return true;
  return false;
}
