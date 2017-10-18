#include <unistd.h>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <vector>
#include <ctime>
#include <string>
#include <climits>
#include "opal_aligner.hpp"



bool OpalAligner::opal_swa(unsigned char * query, int query_length, unsigned char ** db, int db_length, int * db_seq_lengths) {
  // ----------------------------- MAIN CALCULATION ----------------------------- //
  bool retval = false;
  OpalSearchResult** results = new OpalSearchResult*[db_length];
  for (int i = 0; i < db_length; i++) {
    results[i] = new OpalSearchResult;
    opalInitSearchResult(results[i]);
  }
  int resultCode = opalSearchDatabase(query, query_length, db, db_length, db_seq_lengths,
      gap_open, gap_extend, score_matrix, alphabet_length,
      results, search_type, mode_code, overflow_method);
  if (resultCode) {
    printf("\nOpal database search failed with error code: %d\n", resultCode);
  }

  for (int i = 0; i < db_length; i++) {
    if (results[i]->score > 94)
      retval = true;
  }

  for (int i = 0; i < db_length; i++) {
    if (results[i]->alignment) {
      free(results[i]->alignment);
    }
    delete (results[i]);
  }
  delete[] results;
  return retval;
}

