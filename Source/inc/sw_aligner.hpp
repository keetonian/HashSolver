#ifndef SW_ALIGNER_H_
#define SW_ALIGNER_H_

#include "smith_waterman.h"

class SWAligner {
  private:
  const int match = 1;
  const int mismatch = 0;
  const int gap_open = 0;
  const int gap_extend = 0;

  const bool no_start_gap_penalty = false;
  const bool no_end_gap_penalty = false;
  const bool no_gaps_in_read = false;
  const bool no_gaps_in_reference = false;
  const bool no_mismatches = false;
  const bool case_sensitive = true; // All should be upper case

  scoring_t scoring;
  sw_aligner_t * swa;
  alignment_t * alignment;

  public:
  SWAligner();
  void sw_init();
  bool sw_align(const char* read, const char* reference, int threshold);
};

#endif //SW_ALIGNER_H_
