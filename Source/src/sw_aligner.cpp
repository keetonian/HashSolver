#include "sw_aligner.hpp"
#include "smith_waterman.h"
#include <iostream>

SWAligner::SWAligner() {
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
      no_start_gap_penalty, no_end_gap_penalty, 
      no_gaps_in_read, no_gaps_in_reference,
      no_mismatches, case_sensitive);

  alignment = alignment_create(256);
  this->sw_init();
}

void SWAligner::sw_init() {
  swa = smith_waterman_new();
}

bool SWAligner::sw_align(const char* read, const char* reference, int threshold) {
  sw_init();
  smith_waterman_align(read, reference, &scoring, swa);

  // Only pull out the best score.
  smith_waterman_fetch(swa, alignment);
  int score = alignment->score;

  // Makes swa work right. However, I would much rather reset than reallocate.
  // Is there a way to do this?
  smith_waterman_free(swa);

  if (score >= threshold)
    return true;
  return false;
}
