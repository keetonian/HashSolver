#ifndef SHD_FILTER_H_
#define SHD_FILTER_H_

#include "vector_filter.h"

class SHDFilter {
  private:
    char shd_read[128] __aligned__;
    char shd_ref[128] __aligned__;
    char blank_string[128];
    char * hamming_masks;
  public:
    SHDFilter();
    void init(uint32_t error_threshold, uint32_t read_length);
    int filter(const char * read, const char * ref, int error, int read_length);
    int magnet(const char * read, const char * ref, int error, int read_length);
    int magnet_scalar(const char * read, const char * ref, int error, int read_length);
};
#endif //SHD_FILTER_H_
