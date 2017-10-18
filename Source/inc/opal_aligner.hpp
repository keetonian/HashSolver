#ifndef OPAL_ALIGNER_H_
#define OPAL_ALIGNER_H_

#include "opal.h"

class OpalAligner {
  private:
    // Does this need to be an int? Can it be smaller?
    int score_matrix_other[16] = { 
      0, 1, 1, 1,
      1, 0, 1, 1,
      1, 1, 0, 1,
      1, 1, 1, 0 };

    int score_matrix[16] = {
      1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1 };

    // A C G T
    int alphabet_length = 4;

    // Use smith waterman (SW) mode
    int mode_code = OPAL_MODE_SW;

    // Penalties
    int gap_open = 0;

    int gap_extend = 0;

    // Only search for the score.
    // Other search types return alignment information
    int search_type = OPAL_SEARCH_SCORE;

    // Uses characters for alignment scores until overflows.
    int overflow_method = OPAL_OVERFLOW_SIMPLE;

  public:
    //unsigned char * opal_swa_char(unsigned char * query, int query_length, 
    //unsigned char ** db, int db_length, int * db_seq_lengths);

    bool opal_swa(unsigned char * query, int query_length, 
	unsigned char ** db, int db_length, int * db_seq_lengths);
};
#endif // OPAL_ALIGNER_H_
