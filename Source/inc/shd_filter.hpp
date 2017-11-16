#ifndef SHD_FILTER_H_
#define SHD_FILTER_H_

#include "vector_filter.h"

class SHDFilter {
  private:
    char shd_read[128] __aligned__;
    char shd_ref[128] __aligned__;
    char blank_string[128];
  public:
    SHDFilter();
    int filter(const char * read, const char * ref, int error, int read_length);
    int magnet(const char * read, const char * ref, int error, int read_length);
};
#endif //SHD_FILTER_H_
