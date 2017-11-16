/*
 * vector_ed.c
 *
 *  Created on: Nov 8, 2013
 *      Author: hxin
 */

#include "print.h"
#include "vector_filter.h"
#include <nmmintrin.h>
#include <tmmintrin.h>
//#include <x86intrin.h>
#include <stdio.h>
#include <string.h>
#include "popcount.h"
#include "bit_convert.h"
#include "mask.h"
#include "time.h"

uint8_t MASK_01[32] __aligned__ = { 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55,
  0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55,
  0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55,
  0x55 };

/*
 * By little endians, left shift should actually be right shift in x86 convention
 */

__m128i right_alignr_helper(__m128i prev, __m128i curr, int shift_num) {
  switch (shift_num) {
    case 0:
      return _mm_alignr_epi8(curr, prev, 16);
      break;
    case 1:
      return _mm_alignr_epi8(curr, prev, 15);
      break;
    case 2:
      return _mm_alignr_epi8(curr, prev, 14);
      break;
    case 3:
      return _mm_alignr_epi8(curr, prev, 13);
      break;
    case 4:
      return _mm_alignr_epi8(curr, prev, 12);
      break;
    case 5:
      return _mm_alignr_epi8(curr, prev, 11);
      break;
    case 6:
      return _mm_alignr_epi8(curr, prev, 10);
      break;
    case 7:
      return _mm_alignr_epi8(curr, prev, 9);
      break;
    case 8:
      return _mm_alignr_epi8(curr, prev, 8);
      break;
    case 9:
      return _mm_alignr_epi8(curr, prev, 7);
      break;
    case 10:
      return _mm_alignr_epi8(curr, prev, 6);
      break;
    case 11:
      return _mm_alignr_epi8(curr, prev, 5);
      break;
    case 12:
      return _mm_alignr_epi8(curr, prev, 4);
      break;
    case 13:
      return _mm_alignr_epi8(curr, prev, 3);
      break;
    case 14:
      return _mm_alignr_epi8(curr, prev, 2);
      break;
    case 15:
      return _mm_alignr_epi8(curr, prev, 1);
      break;
    default:
      printf("Error! Invalid shift! From vector_fiter.c: shift_num: %d\n",
	  shift_num);
      exit(1);
      break;
  }
}

__m128i left_alignr_helper(__m128i curr, __m128i next, int shift_num) {
  switch (shift_num) {
    case 0:
      return _mm_alignr_epi8(next, curr, 0);
      break;
    case 1:
      return _mm_alignr_epi8(next, curr, 1);
      break;
    case 2:
      return _mm_alignr_epi8(next, curr, 2);
      break;
    case 3:
      return _mm_alignr_epi8(next, curr, 3);
      break;
    case 4:
      return _mm_alignr_epi8(next, curr, 4);
      break;
    case 5:
      return _mm_alignr_epi8(next, curr, 5);
      break;
    case 6:
      return _mm_alignr_epi8(next, curr, 6);
      break;
    case 7:
      return _mm_alignr_epi8(next, curr, 7);
      break;
    case 8:
      return _mm_alignr_epi8(next, curr, 8);
      break;
    case 9:
      return _mm_alignr_epi8(next, curr, 9);
      break;
    case 10:
      return _mm_alignr_epi8(next, curr, 10);
      break;
    case 11:
      return _mm_alignr_epi8(next, curr, 11);
      break;
    case 12:
      return _mm_alignr_epi8(next, curr, 12);
      break;
    case 13:
      return _mm_alignr_epi8(next, curr, 13);
      break;
    case 14:
      return _mm_alignr_epi8(next, curr, 14);
      break;
    case 15:
      return _mm_alignr_epi8(next, curr, 15);
      break;
    default:
      printf("Error! Invalid shift! From vector_fiter.c: shift_num: %d\n",
	  shift_num);
      exit(1);
      break;
  }
}

__m128i shift_right_sse1(__m128i vec, int shift_num) {
  if (shift_num >= 64) {
    vec = _mm_slli_si128(vec, 8);
    return _mm_slli_epi64(vec, shift_num % 64);
  }
  else {
    __m128i carryover = _mm_slli_si128(vec, 8);
    carryover = _mm_srli_epi64(carryover, 64 - shift_num);
    vec = _mm_slli_epi64(vec, shift_num);
    return _mm_or_si128(vec, carryover);
  }
}

__m128i shift_left_sse1(__m128i vec, int shift_num) {
  if (shift_num >= 64) {
    vec = _mm_srli_si128(vec, 8);
    return _mm_srli_epi64(vec, shift_num % 64);
  }
  else {
    __m128i carryover = _mm_srli_si128(vec, 8);
    carryover = _mm_slli_epi64(carryover, 64 - shift_num);
    vec = _mm_srli_epi64(vec, shift_num);
    return _mm_or_si128(vec, carryover);
  }
}

__m128i shift_right_sse11(__m128i pri_vec, __m128i vec, int shift_num) {
  if (shift_num % 4 == 0)
    return right_alignr_helper(pri_vec, vec, shift_num / 4);

  __m128i carryover;
  __m128i shiftee;
  __m128i mask;

  carryover = right_alignr_helper(pri_vec, vec, shift_num / 4 + 1);
  //	print128_bit(carryover);
  carryover = _mm_srli_epi64(carryover, (4 - (shift_num % 4)) * 2);
  //	print128_bit(carryover);

  if (shift_num > 4)
    shiftee = right_alignr_helper(pri_vec, vec, shift_num / 4);
  else
    shiftee = vec;

  shiftee = _mm_slli_epi64(shiftee, (shift_num % 4) * 2);

  return _mm_or_si128(shiftee, carryover);
}

__m128i shift_left_sse11(__m128i vec, __m128i next_vec, int shift_num) {
  if (shift_num % 4 == 0)
    return left_alignr_helper(vec, next_vec, shift_num / 4);

  __m128i carryover;
  __m128i shiftee;
  __m128i mask;

  carryover = left_alignr_helper(vec, next_vec, shift_num / 4 + 1);
  //	print128_bit(carryover);
  carryover = _mm_slli_epi64(carryover, (4 - (shift_num % 4)) * 2);
  //	print128_bit(carryover);

  if (shift_num > 4)
    shiftee = left_alignr_helper(vec, next_vec, shift_num / 4);
  else
    shiftee = vec;

  shiftee = _mm_srli_epi64(shiftee, (shift_num % 4) * 2);

  return _mm_or_si128(shiftee, carryover);
}

__m128i xor11complement_sse(__m128i input) {
  __m128i temp, result;
  __m128i mask = _mm_load_si128((__m128i *) MASK_01);

  temp = _mm_and_si128(input, mask);
  temp = _mm_slli_epi16(temp, 1);
  result = _mm_or_si128(temp, input);

  mask = _mm_slli_epi16(mask, 1);
  temp = _mm_and_si128(input, mask);
  temp = _mm_srli_epi16(temp, 1);
  result = _mm_or_si128(result, temp);

  return result;
}

void flip_false_zero(__m128i& vec) {

  //	printf("vec: \t\t");
  //	print128_bit(vec);

  //For not crossing bits
  __m128i *boundary= (__m128i *) MASK_7F;
  //	printf("MASK_7F: \t");
  //	print128_bit(*boundary);

  __m128i shift = _mm_and_si128(*boundary, vec);

  //	printf("After and: \t");
  //	print128_bit(shift);

  __m128i *mask = (__m128i *) MASK_0TO1;

  shift = _mm_shuffle_epi8(*mask, shift);
  vec = _mm_or_si128(vec, shift);

  //	printf("Last cases %d: \t", 0);
  //	print128_bit(vec);

  int i;
  for (i = 1; i < 4; i++) {
    shift = _mm_srli_epi16(vec, i);
    shift = _mm_and_si128(*boundary, shift);
    //		printf("shift %d: \t", i);
    //		print128_bit(shift);
    shift = _mm_shuffle_epi8(*mask, shift);
    //		printf("shuffle %d: \t", i);
    //		print128_bit(shift);
    shift = _mm_slli_epi16(shift, i);
    vec = _mm_or_si128(vec, shift);
    //		printf("Last cases %d: \t", i);
    //		print128_bit(vec);
  }

  //For the crossing bits
  __m128i shifted_vec = shift_right_sse1(vec, 4);
  //	printf("shifted_vec: \t");
  //	print128_bit(shifted_vec);

  shift = _mm_and_si128(*boundary, shifted_vec);
  //	printf("After and: \t");
  //	print128_bit(shift);

  shift = _mm_shuffle_epi8(*mask, shift);
  shifted_vec = _mm_or_si128(shifted_vec, shift);
  //	printf("Cross cases %d: \t", 0);
  //	print128_bit(shifted_vec);

  for (i = 1; i < 4; i++) {
    shift = _mm_srli_epi16(shifted_vec, i);
    shift = _mm_and_si128(*boundary, shift);
    shift = _mm_shuffle_epi8(*mask, shift);
    shift = _mm_slli_epi16(shift, i);
    shifted_vec = _mm_or_si128(shifted_vec, shift);
    //		printf("Cross cases %d: \t", i);
    //		print128_bit(shifted_vec);
  }

  shifted_vec = shift_left_sse1(shifted_vec, 4);
  vec = _mm_or_si128(shifted_vec, vec);
  //	printf("Final case: \t");
  //	print128_bit(vec);


}

int bit_vec_filter_m128_sse1(uint8_t *read_vec0, uint8_t *read_vec1, uint8_t
    *ref_vec0, uint8_t *ref_vec1, __m128i mask, int max_error) {

  int total_difference = 0;

  //Start iteration
  int j;
  //read data
  __m128i read_XMM0 = *((__m128i *) (read_vec0));
  __m128i read_XMM1 = *((__m128i *) (read_vec1));
  //ref data
  __m128i ref_XMM0 = *((__m128i *) (ref_vec0));
  __m128i ref_XMM1 = *((__m128i *) (ref_vec1));

  __m128i shift_XMM;
  __m128i diff_XMM;
  __m128i temp_diff_XMM;
  __m128i temp_shift_XMM;
  __m128i temp_mask;

  diff_XMM = _mm_xor_si128(read_XMM0, ref_XMM0);
  temp_diff_XMM = _mm_xor_si128(read_XMM1, ref_XMM1);
  //print128_bit(diff_XMM);
  //print128_bit(temp_diff_XMM);
  diff_XMM = _mm_or_si128(diff_XMM, temp_diff_XMM);
  //print128_bit(diff_XMM);
  //printf("\n");

  total_difference = popcount1_m128i_sse(diff_XMM);
  if (total_difference < max_error)
    return 2;

  flip_false_zero(diff_XMM);

  //	printf("diff_XMM: \t");
  //	print128_bit(diff_XMM);

  for (j = 1; j <= max_error; j++) {
    temp_mask = _mm_load_si128( (__m128i *) (MASK_SSE_BEG1 + (j - 1) *
	  SSE_BYTE_NUM));
    temp_mask = _mm_and_si128(temp_mask, mask);

    //Right shift read
    shift_XMM = shift_right_sse1(read_XMM0, j);
    temp_diff_XMM = _mm_xor_si128(shift_XMM, ref_XMM0);
    shift_XMM = shift_right_sse1(read_XMM1, j);
    temp_shift_XMM = _mm_xor_si128(shift_XMM, ref_XMM1);
    temp_diff_XMM = _mm_or_si128(temp_shift_XMM, temp_diff_XMM);
    temp_diff_XMM = _mm_and_si128(temp_diff_XMM, temp_mask);
    //		printf("Before flip: \t");
    //		print128_bit(temp_diff_XMM);
    flip_false_zero(temp_diff_XMM);
    //		printf("After flip: \t");
    //		print128_bit(temp_diff_XMM);
    diff_XMM = _mm_and_si128(diff_XMM, temp_diff_XMM);

    //		printf("diff_XMM: \t");
    //		print128_bit(diff_XMM);

    //Right shift ref
    shift_XMM = shift_right_sse1(ref_XMM0, j);
    temp_diff_XMM = _mm_xor_si128(shift_XMM, read_XMM0);
    shift_XMM = shift_right_sse1(ref_XMM1, j);
    temp_shift_XMM = _mm_xor_si128(shift_XMM, read_XMM1);
    temp_diff_XMM = _mm_or_si128(temp_shift_XMM, temp_diff_XMM);
    temp_diff_XMM = _mm_and_si128(temp_diff_XMM, temp_mask);
    //		printf("Before flip: \t");
    //		print128_bit(temp_diff_XMM);
    flip_false_zero(temp_diff_XMM);
    //		printf("After flip: \t");
    //		print128_bit(temp_diff_XMM);
    diff_XMM = _mm_and_si128(diff_XMM, temp_diff_XMM);

    //		printf("diff_XMM: \t");
    //		print128_bit(diff_XMM);
  }

  total_difference = popcount1_m128i_sse(diff_XMM);

  //	printf("total_difference: %d\n", total_difference);

  if (total_difference > (max_error) )
    return 0;
  else
    return 1;
}

// First bit, last 27 bits set high
uint8_t MAGNET_MASK[16] __aligned__ = { 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0xe0, 0xff, 0xff, 0xff };

uint8_t MAGNET_MASK_2[16] __aligned__ = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0xff, 0xff, 0xff };

uint8_t ALL_BITS_SET[16] __aligned__ = { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
  0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff };

int bit_magnet_filter_m128_sse1(uint8_t *read_vec0, uint8_t *read_vec1, uint8_t
    *ref_vec0, uint8_t *ref_vec1, __m128i mask, int max_error) {

  int total_difference = 0;
  int longest_subsequence = 0;
  uint32_t index_of_mask = 0;
  uint32_t read_length = 100;

  //Start iteration
  int j;
  int num_vectors = max_error*2 + 1;
  //read data
  __m128i read_XMM0 = *((__m128i *) (read_vec0));
  __m128i read_XMM1 = *((__m128i *) (read_vec1));
  //ref data
  __m128i ref_XMM0 = *((__m128i *) (ref_vec0));
  __m128i ref_XMM1 = *((__m128i *) (ref_vec1));

  __m128i shift_XMM;
  __m128i diff_XMM[num_vectors];
  __m128i temp_diff_XMM;
  __m128i temp_shift_XMM;
  __m128i temp_mask;

  diff_XMM[0] = _mm_xor_si128(read_XMM0, ref_XMM0);
  temp_diff_XMM = _mm_xor_si128(read_XMM1, ref_XMM1);
  diff_XMM[0] = _mm_or_si128(diff_XMM[0], temp_diff_XMM);

  total_difference = popcount1_m128i_sse(diff_XMM[0]);
  //printf("Total difference: %d\n", total_difference);

  if (total_difference <= max_error)
    return 2;
  else if(max_error == 0)
    return 0;

  //	printf("diff_XMM: \t");
  //	print128_bit(diff_XMM);

  for (j = 1; j <= max_error; j++) {
    temp_mask = _mm_load_si128( (__m128i *) (MASK_SSE_BEG1 + (j - 1) *
	  SSE_BYTE_NUM));
    temp_mask = _mm_and_si128(temp_mask, mask);

    //Right shift read
    shift_XMM = shift_right_sse1(read_XMM0, j);
    temp_diff_XMM = _mm_xor_si128(shift_XMM, ref_XMM0);
    shift_XMM = shift_right_sse1(read_XMM1, j);
    temp_shift_XMM = _mm_xor_si128(shift_XMM, ref_XMM1);
    temp_diff_XMM = _mm_or_si128(temp_shift_XMM, temp_diff_XMM);
    diff_XMM[j] = _mm_and_si128(temp_diff_XMM, temp_mask);
    // Shift back so that starts with leftmost bit
    diff_XMM[j] = shift_left_sse1(diff_XMM[j], j);

/*    total_difference = popcount1_m128i_sse(diff_XMM[j]);
    if (total_difference <= max_error)
      return 2;
*/
    //Right shift ref
    shift_XMM = shift_right_sse1(ref_XMM0, j);
    temp_diff_XMM = _mm_xor_si128(shift_XMM, read_XMM0);
    shift_XMM = shift_right_sse1(ref_XMM1, j);
    temp_shift_XMM = _mm_xor_si128(shift_XMM, read_XMM1);
    temp_diff_XMM = _mm_or_si128(temp_shift_XMM, temp_diff_XMM);
    diff_XMM[num_vectors - j] = _mm_and_si128(temp_diff_XMM, temp_mask);
    // Shift back so that starts with leftmost bit
    diff_XMM[num_vectors - j] = shift_left_sse1(diff_XMM[num_vectors - j], j);
/*
    total_difference = popcount1_m128i_sse(diff_XMM[num_vectors - j]);
    if (total_difference <= max_error)
      return 2;
*/
  }

  // Create result register.
  __m128i result_XMM = _mm_load_si128((__m128i *)ALL_BITS_SET);

  // Create register with all bits set high
  temp_diff_XMM = _mm_or_si128(result_XMM, result_XMM);

  // Load magnet mask into temp_mask register
  temp_mask = _mm_load_si128( (__m128i *) MAGNET_MASK_2);
  int i = 0, k = 0, res = 0;

  // Prepare bit vectors by shifting over 1, setting first and last bits.
  for (j = 0; j < num_vectors; j++) {
    //print128_bit(diff_XMM[j]);
    diff_XMM[j] = _mm_or_si128(temp_mask, diff_XMM[j]);
    //print128_bit(diff_XMM[j]);
  }
  uint32_t cumulative_subsequences = 0;

  for (j = 0; j <= max_error; j++) {
    // Find longest subsequence of 0's
    longest_subsequence = 0;
    for ( k = 0; k < num_vectors; k++ ) {
      // Shifting right works better with AND operations on the inverse,
      // Shiftling left works better with OR operations
      // Flip all the bits
      temp_shift_XMM = _mm_xor_si128(diff_XMM[k], temp_diff_XMM);
      //temp_shift_XMM = _mm_or_si128(shift_right_sse1(diff_XMM[k], 1), diff_XMM[k]);
      //print128_bit(temp_shift_XMM);
      int subsequence = 0;
      // DO BINARY SEARCH
      //clock_t cpu_time1 = clock();
      while(!_mm_test_all_zeros(temp_shift_XMM, temp_diff_XMM)) {
	// Shift until true.
	temp_shift_XMM = _mm_and_si128(shift_right_sse1(temp_shift_XMM, 1), temp_shift_XMM);
	subsequence++;
	//printf("S: %d\n", subsequence);
	//print128_bit(temp_shift_XMM);
      }
      //clock_t cpu_time2 = clock();
      //total_shift_time += (double)(cpu_time2 - cpu_time1);

      longest_subsequence = (subsequence>=longest_subsequence)*subsequence + (subsequence<longest_subsequence)*longest_subsequence;
      index_of_mask = (subsequence>=longest_subsequence)*k + (subsequence<longest_subsequence)*index_of_mask;

    }

    if (longest_subsequence <= 1) 
      break;
    if (longest_subsequence + max_error < (read_length - cumulative_subsequences)/(max_error + 1 - j))
      return 0;

    cumulative_subsequences += longest_subsequence;
    //printf("Longest subsequence: %d\n", longest_subsequence);
    //printf("Mask: \t");
    //print128_bit(diff_XMM[index_of_mask]);
    // Create mask for subsequence
    // USE BINARY METHOD

    // Invert bits
    temp_shift_XMM = _mm_xor_si128(diff_XMM[index_of_mask], temp_diff_XMM);
    // Want the last 0 to stay, not be overwritten to 1.
    for (k = 0; k < longest_subsequence-1; k++) {
      temp_shift_XMM = _mm_and_si128(shift_right_sse1(temp_shift_XMM, 1), temp_shift_XMM);
    }
    //print128_bit(temp_shift_XMM);
    for(k = 1; k < longest_subsequence; k++) {
      temp_shift_XMM = _mm_or_si128(shift_left_sse1(temp_shift_XMM, 1), temp_shift_XMM);
    }
    //print128_bit(temp_shift_XMM);

    // Invert back to normal
    temp_shift_XMM = _mm_xor_si128(temp_shift_XMM, temp_diff_XMM);
    //print128_bit(diff_XMM[index_of_mask]);
    //print128_bit(temp_shift_XMM);

    // Update result
    result_XMM = _mm_and_si128(temp_shift_XMM, result_XMM);

    // Set bits around longest subsequence to 0 as well
    temp_shift_XMM = _mm_xor_si128(temp_shift_XMM, temp_diff_XMM);
    temp_shift_XMM = _mm_or_si128(shift_right_sse1(temp_shift_XMM, 1), temp_shift_XMM);
    temp_shift_XMM = _mm_or_si128(shift_left_sse1(temp_shift_XMM, 1), temp_shift_XMM);
    // Make sure rightmost bits (and leftmost bit) are still set.
    //temp_shift_XMM = _mm_or_si128(temp_shift_XMM, temp_mask);
    // XOR to flip all bits
    //print128_bit(temp_shift_XMM);

    // Set bits for longest subsequence in all other masks high
    for (k = 0; k < num_vectors; k++) {
      //print128_bit(diff_XMM[k]);
      diff_XMM[k] = _mm_or_si128(temp_shift_XMM, diff_XMM[k]);
      //print128_bit(diff_XMM[k]);
      //printf("\n");
    }
  }

  // Remove the mask
  //printf("before mask removal:\n");
  //print128_bit(result_XMM);

  result_XMM = _mm_andnot_si128(temp_mask, result_XMM);

  //printf("After mask removal:\n");
  //print128_bit(result_XMM);
  // Shift back into correct position (shift left 1)
  //printf("Before shifting:\n");
  //print128_bit(result_XMM);
  //result_XMM = shift_left_sse1(result_XMM, 1);
  //printf("After shifting:\n");
  //print128_bit(result_XMM);

  // Find the difference
  total_difference = popcount1_m128i_sse(result_XMM);

  //printf("total_difference: %d\n", total_difference);

  if (total_difference > (max_error) )
    return 0;
  else
    return 1;
}

int bit_magnet_bin_filter_m128_sse1(uint8_t *read_vec0, uint8_t *read_vec1, uint8_t
    *ref_vec0, uint8_t *ref_vec1, __m128i mask, int max_error) {

  int total_difference = 0;
  int longest_subsequence = 0;
  uint32_t index_of_mask = 0;

  //Start iteration
  int j;
  int num_vectors = max_error*2 + 1;
  //read data
  __m128i read_XMM0 = *((__m128i *) (read_vec0));
  __m128i read_XMM1 = *((__m128i *) (read_vec1));
  //ref data
  __m128i ref_XMM0 = *((__m128i *) (ref_vec0));
  __m128i ref_XMM1 = *((__m128i *) (ref_vec1));

  __m128i shift_XMM;
  __m128i diff_XMM[num_vectors];
  __m128i temp_diff_XMM;
  __m128i temp_shift_XMM;
  __m128i temp_mask;

  diff_XMM[0] = _mm_xor_si128(read_XMM0, ref_XMM0);
  temp_diff_XMM = _mm_xor_si128(read_XMM1, ref_XMM1);
  diff_XMM[0] = _mm_or_si128(diff_XMM[0], temp_diff_XMM);

  total_difference = popcount1_m128i_sse(diff_XMM[0]);
  //printf("Total difference: %d\n", total_difference);

  if (total_difference <= max_error)
    return 1;
  else if(max_error == 0)
    return 0;

  	printf("diff_XMM: \t");
  //	print128_bit(diff_XMM);

  for (j = 1; j <= max_error; j++) {
    temp_mask = _mm_load_si128( (__m128i *) (MASK_SSE_BEG1 + (j - 1) *
	  SSE_BYTE_NUM));
    temp_mask = _mm_and_si128(temp_mask, mask);

    //Right shift read
    shift_XMM = shift_right_sse1(read_XMM0, j);
    temp_diff_XMM = _mm_xor_si128(shift_XMM, ref_XMM0);
    shift_XMM = shift_right_sse1(read_XMM1, j);
    temp_shift_XMM = _mm_xor_si128(shift_XMM, ref_XMM1);
    temp_diff_XMM = _mm_or_si128(temp_shift_XMM, temp_diff_XMM);
    diff_XMM[j] = _mm_and_si128(temp_diff_XMM, temp_mask);
    // Shift back so that starts with leftmost bit
    diff_XMM[j] = shift_left_sse1(diff_XMM[j], j);

    total_difference = popcount1_m128i_sse(diff_XMM[j]);
    if (total_difference <= max_error)
      return 1;

    //Right shift ref
    shift_XMM = shift_right_sse1(ref_XMM0, j);
    temp_diff_XMM = _mm_xor_si128(shift_XMM, read_XMM0);
    shift_XMM = shift_right_sse1(ref_XMM1, j);
    temp_shift_XMM = _mm_xor_si128(shift_XMM, read_XMM1);
    temp_diff_XMM = _mm_or_si128(temp_shift_XMM, temp_diff_XMM);
    diff_XMM[num_vectors - j] = _mm_and_si128(temp_diff_XMM, temp_mask);
    // Shift back so that starts with leftmost bit
    diff_XMM[num_vectors - j] = shift_left_sse1(diff_XMM[num_vectors - j], j);

    total_difference = popcount1_m128i_sse(diff_XMM[num_vectors - j]);
    if (total_difference <= max_error)
      return 1;
  }

  // Create result register.
  __m128i result_XMM = _mm_load_si128((__m128i *)ALL_BITS_SET);

  // Create register with all bits set high
  temp_diff_XMM = _mm_or_si128(result_XMM, result_XMM);

  // Load magnet mask into temp_mask register
  temp_mask = _mm_load_si128( (__m128i *) MAGNET_MASK_2);
  int i = 0, k = 0, res = 0;

  // Prepare bit vectors by shifting over 1, setting first and last bits.
  for (j = 0; j < num_vectors; j++) {
    //print128_bit(diff_XMM[j]);
    diff_XMM[j] = _mm_or_si128(temp_mask, diff_XMM[j]);
    //print128_bit(diff_XMM[j]);
  }

  for (j = 0; j <= max_error; j++) {
    // Find longest subsequence of 0's
    longest_subsequence = 0;
    for ( k = 0; k < num_vectors; k++ ) {
      // Shifting right works better with AND operations on the inverse,
      // Shiftling left works better with OR operations
      // Flip all the bits
      temp_shift_XMM = _mm_xor_si128(diff_XMM[k], temp_diff_XMM);
      //temp_shift_XMM = _mm_or_si128(shift_right_sse1(diff_XMM[k], 1), diff_XMM[k]);
      //print128_bit(temp_shift_XMM);
      int subsequence = 0;
      int shift_width = 1;
      int all_zero = 0;
      int max = 128;
      // DO BINARY SEARCH
      while(max-subsequence>1 /*&& max > longest_subsequence*/) {
	// Store previous result.
	shift_XMM = _mm_or_si128(temp_shift_XMM, temp_shift_XMM);
	// Shift until true.
	temp_shift_XMM = _mm_and_si128(shift_right_sse1(temp_shift_XMM, shift_width), temp_shift_XMM);
	if(_mm_test_all_zeros(temp_shift_XMM, temp_diff_XMM)){
	  temp_shift_XMM = _mm_or_si128(shift_XMM, shift_XMM);
	  max = subsequence + shift_width;
	  shift_width = shift_width>>1;
	  printf("%d\n", shift_width);
	  continue;
	}
	//printf("max: %d/%d\n", max, subsequence);

	subsequence += shift_width;
	shift_width = shift_width<<1;
	printf("%d\n", shift_width);
	//printf("S: %d\n", subsequence);
	//print128_bit(temp_shift_XMM);
      }
      printf("%d/%d\n", subsequence, max);
      subsequence = max;

      longest_subsequence = (subsequence>=longest_subsequence)*subsequence + (subsequence<longest_subsequence)*longest_subsequence;
      index_of_mask = (subsequence>=longest_subsequence)*k + (subsequence<longest_subsequence)*index_of_mask;

    }
    if (longest_subsequence <= 1)
      break;
    //printf("Longest subsequence: %d\n", longest_subsequence);
    //printf("Mask: \t");
    //print128_bit(diff_XMM[index_of_mask]);
    // Create mask for subsequence
    // USE BINARY METHOD

    // Invert bits
    temp_shift_XMM = _mm_xor_si128(diff_XMM[index_of_mask], temp_diff_XMM);
    // Want the last 0 to stay, not be overwritten to 1.
    for (k = 0; k < longest_subsequence-1; k++) {
      temp_shift_XMM = _mm_and_si128(shift_right_sse1(temp_shift_XMM, 1), temp_shift_XMM);
    }
    //print128_bit(temp_shift_XMM);
    for(k = 1; k < longest_subsequence; k++) {
      temp_shift_XMM = _mm_or_si128(shift_left_sse1(temp_shift_XMM, 1), temp_shift_XMM);
    }
    //print128_bit(temp_shift_XMM);

    // Invert back to normal
    temp_shift_XMM = _mm_xor_si128(temp_shift_XMM, temp_diff_XMM);
    //print128_bit(diff_XMM[index_of_mask]);
    //print128_bit(temp_shift_XMM);

    // Update result
    result_XMM = _mm_and_si128(temp_shift_XMM, result_XMM);

    // Set bits around longest subsequence to 0 as well
    temp_shift_XMM = _mm_xor_si128(temp_shift_XMM, temp_diff_XMM);
    temp_shift_XMM = _mm_or_si128(shift_right_sse1(temp_shift_XMM, 1), temp_shift_XMM);
    temp_shift_XMM = _mm_or_si128(shift_left_sse1(temp_shift_XMM, 1), temp_shift_XMM);
    // Make sure rightmost bits (and leftmost bit) are still set.
    //temp_shift_XMM = _mm_or_si128(temp_shift_XMM, temp_mask);
    // XOR to flip all bits
    //print128_bit(temp_shift_XMM);

    // Set bits for longest subsequence in all other masks high
    for (k = 0; k < num_vectors; k++) {
      //print128_bit(diff_XMM[k]);
      diff_XMM[k] = _mm_or_si128(temp_shift_XMM, diff_XMM[k]);
      //print128_bit(diff_XMM[k]);
      //printf("\n");
    }
  }

  // Remove the mask
  //printf("before mask removal:\n");
  //print128_bit(result_XMM);
  result_XMM = _mm_andnot_si128(temp_mask, result_XMM);
  //printf("After mask removal:\n");
  //print128_bit(result_XMM);
  // Shift back into correct position (shift left 1)
  //printf("Before shifting:\n");
  //print128_bit(result_XMM);
  //result_XMM = shift_left_sse1(result_XMM, 1);
  //printf("After shifting:\n");
  //print128_bit(result_XMM);
  // Find the difference
  total_difference = popcount1_m128i_sse(result_XMM);

  printf("total_difference: %d\n", total_difference);

  if (total_difference > (max_error) )
    return 0;
  else
    return 1;
}

int bit_vec_filter_no_flipping_m128_sse1(uint8_t *read_vec0, uint8_t *read_vec1, uint8_t
    *ref_vec0, uint8_t *ref_vec1, __m128i mask, int max_error) {

  int total_difference = 0;

  //Start iteration
  int j;
  //read data
  __m128i read_XMM0 = *((__m128i *) (read_vec0));
  __m128i read_XMM1 = *((__m128i *) (read_vec1));
  //ref data
  __m128i ref_XMM0 = *((__m128i *) (ref_vec0));
  __m128i ref_XMM1 = *((__m128i *) (ref_vec1));

  __m128i shift_XMM;
  __m128i diff_XMM;
  __m128i temp_diff_XMM;
  __m128i temp_shift_XMM;
  __m128i temp_mask;

  diff_XMM = _mm_xor_si128(read_XMM0, ref_XMM0);
  temp_diff_XMM = _mm_xor_si128(read_XMM1, ref_XMM1);
  diff_XMM = _mm_or_si128(diff_XMM, temp_diff_XMM);

  //printf("diff_XMM: \n");
  //print128_bit_twice(diff_XMM);

  for (j = 1; j <= max_error; j++) {
    temp_mask = _mm_load_si128( (__m128i *) (MASK_SSE_BEG1 + (j - 1) *
	  SSE_BYTE_NUM));
    temp_mask = _mm_and_si128(temp_mask, mask);

    //Right shift read
    shift_XMM = shift_right_sse1(read_XMM0, j);
    temp_diff_XMM = _mm_xor_si128(shift_XMM, ref_XMM0);
    shift_XMM = shift_right_sse1(read_XMM1, j);
    temp_shift_XMM = _mm_xor_si128(shift_XMM, ref_XMM1);
    temp_diff_XMM = _mm_or_si128(temp_shift_XMM, temp_diff_XMM);
    temp_diff_XMM = _mm_and_si128(temp_diff_XMM, temp_mask);
    //		printf("Before flip: \t");
    //		print128_bit(temp_diff_XMM);
    //		flip_false_zero(temp_diff_XMM); //No flipping
    //		printf("After flip: \t");
    //		print128_bit(temp_diff_XMM);
    diff_XMM = _mm_and_si128(diff_XMM, temp_diff_XMM);

    //printf("read shift %d diff_XMM: \n", j);
    //print128_bit_twice(diff_XMM);

    //Right shift ref
    shift_XMM = shift_right_sse1(ref_XMM0, j);
    temp_diff_XMM = _mm_xor_si128(shift_XMM, read_XMM0);
    shift_XMM = shift_right_sse1(ref_XMM1, j);
    temp_shift_XMM = _mm_xor_si128(shift_XMM, read_XMM1);
    temp_diff_XMM = _mm_or_si128(temp_shift_XMM, temp_diff_XMM);
    temp_diff_XMM = _mm_and_si128(temp_diff_XMM, temp_mask);
    //		printf("Before flip: \t");
    //		print128_bit(temp_diff_XMM);
    //		flip_false_zero(temp_diff_XMM); //No flipping
    //		printf("After flip: \t");
    //		print128_bit(temp_diff_XMM);
    diff_XMM = _mm_and_si128(diff_XMM, temp_diff_XMM);

    //printf("ref shift %d diff_XMM: \n", j);
    //print128_bit_twice(diff_XMM);
  }

  total_difference = popcount1_m128i_sse(diff_XMM);

  if (total_difference > max_error)
    return 0;
  else
    return 1;
}

int bit_vec_filter_m128_sse11(uint8_t *read_vec, uint8_t *ref_vec, int length,
    int max_error) {
  const __m128i zero_mask = _mm_set1_epi8(0x00);
  const __m128i one_mask = _mm_set1_epi8(0xff);

  int total_byte = (length - 1) / BYTE_BASE_NUM11 + 1;

  int total_difference = 0;

  //Start iteration
  int i, j;
  //read data
  __m128i prev_read_XMM = _mm_set1_epi8(0x0);
  __m128i curr_read_XMM = *((__m128i *) (read_vec));
  //ref data
  __m128i prev_ref_XMM = _mm_set1_epi8(0x0);
  __m128i curr_ref_XMM = *((__m128i *) (ref_vec));

  __m128i read_XMM;
  __m128i ref_XMM;
  __m128i temp_diff_XMM;
  __m128i diff_XMM;
  __m128i mask;
  for (i = 0; i < total_byte; i += SSE_BYTE_NUM) {

    curr_read_XMM = *((__m128i *) (read_vec + i));
    curr_ref_XMM = *((__m128i *) (ref_vec + i));

    diff_XMM = _mm_xor_si128(curr_read_XMM, curr_ref_XMM);
    diff_XMM = xor11complement_sse(diff_XMM);


    if (i + SSE_BYTE_NUM >= total_byte) {
      if (length % SSE_BASE_NUM11) {
	mask = _mm_load_si128(
	    (__m128i *) (MASK_SSE_END11
	      + (length % SSE_BASE_NUM11) * SSE_BYTE_NUM));
	diff_XMM = _mm_and_si128(mask, diff_XMM);
      }
    }

    for (j = 1; j <= max_error; j++) {
      //Right shift read
      read_XMM = shift_right_sse11(prev_read_XMM, curr_read_XMM, j);

      temp_diff_XMM = _mm_xor_si128(read_XMM, curr_ref_XMM);
      temp_diff_XMM = xor11complement_sse(temp_diff_XMM);

      if (i == 0) {
	mask = _mm_load_si128(
	    (__m128i *) (MASK_SSE_BEG11 + (j - 1) * SSE_BYTE_NUM));

	temp_diff_XMM = _mm_and_si128(mask, temp_diff_XMM);
      }
      if (i + SSE_BYTE_NUM >= total_byte) {
	if (length % SSE_BASE_NUM11) {
	  mask = _mm_load_si128(
	      (__m128i *) (MASK_SSE_END11
		+ (length % SSE_BASE_NUM11) * SSE_BYTE_NUM));
	  temp_diff_XMM = _mm_and_si128(mask, temp_diff_XMM);
	}
      }

      diff_XMM = _mm_and_si128(diff_XMM, temp_diff_XMM);

      //Right shift ref
      ref_XMM = shift_right_sse11(prev_ref_XMM, curr_ref_XMM, j);

      temp_diff_XMM = _mm_xor_si128(curr_read_XMM, ref_XMM);
      temp_diff_XMM = xor11complement_sse(temp_diff_XMM);

      if (i == 0) {
	mask = _mm_load_si128(
	    (__m128i *) (MASK_SSE_BEG11 + (j - 1) * SSE_BYTE_NUM));

	temp_diff_XMM = _mm_and_si128(mask, temp_diff_XMM);
      }
      if (i + SSE_BYTE_NUM >= total_byte) {
	if (length % SSE_BASE_NUM11) {
	  mask = _mm_load_si128(
	      (__m128i *) (MASK_SSE_END11
		+ (length % SSE_BASE_NUM11) * SSE_BYTE_NUM));
	  temp_diff_XMM = _mm_and_si128(mask, temp_diff_XMM);
	}
      }

      diff_XMM = _mm_and_si128(diff_XMM, temp_diff_XMM);
    }

    total_difference += popcount11_m128i_sse(diff_XMM);

    prev_read_XMM = curr_read_XMM;
    prev_ref_XMM = curr_ref_XMM;

    if (total_difference > max_error)
      return 0;
  }

  return 1;
}

#define _MAX_LENGTH_ 320

uint8_t read_bit_t[_MAX_LENGTH_ / 4 + 16] __aligned__;
uint8_t ref_bit_t[_MAX_LENGTH_ / 4 + 16] __aligned__;

uint8_t read_vec0_t[SSE_BYTE_NUM] __aligned__;
uint8_t read_vec1_t[SSE_BYTE_NUM] __aligned__;
uint8_t ref_vec0_t[SSE_BYTE_NUM] __aligned__;
uint8_t ref_vec1_t[SSE_BYTE_NUM] __aligned__;

int bit_magnet_bin_filter_sse1(char* read, char* ref, int length, int max_error) {
  //Get ready the bits
  sse3_convert2bit1(read, read_vec0_t, read_vec1_t);
  sse3_convert2bit1(ref, ref_vec0_t, ref_vec1_t);

  //Get the mask
  __m128i mask;
  if (length >= SSE_BASE_NUM1)
    mask = _mm_set1_epi8(0xff);
  else
    mask = _mm_load_si128( (__m128i *) (MASK_SSE_END1 + (length *
	    SSE_BYTE_NUM)));

  return bit_magnet_bin_filter_m128_sse1(read_vec0_t, read_vec1_t,
      ref_vec0_t, ref_vec1_t, mask, max_error);
}

int bit_magnet_filter_sse1(char* read, char* ref, int length, int max_error) {
  //Get ready the bits
  sse3_convert2bit1(read, read_vec0_t, read_vec1_t);
  sse3_convert2bit1(ref, ref_vec0_t, ref_vec1_t);

  //Get the mask
  __m128i mask;
  if (length >= SSE_BASE_NUM1)
    mask = _mm_set1_epi8(0xff);
  else
    mask = _mm_load_si128( (__m128i *) (MASK_SSE_END1 + (length *
	    SSE_BYTE_NUM)));

  return bit_magnet_filter_m128_sse1(read_vec0_t, read_vec1_t,
      ref_vec0_t, ref_vec1_t, mask, max_error);
}

int bit_vec_filter_sse1(char* read, char* ref, int length, int max_error) {
  //Get ready the bits
  sse3_convert2bit1(read, read_vec0_t, read_vec1_t);
  sse3_convert2bit1(ref, ref_vec0_t, ref_vec1_t);

  //Get the mask
  __m128i mask;
  if (length >= SSE_BASE_NUM1)
    mask = _mm_set1_epi8(0xff);
  else
    mask = _mm_load_si128( (__m128i *) (MASK_SSE_END1 + (length *
	    SSE_BYTE_NUM)));

  return bit_vec_filter_m128_sse1(read_vec0_t, read_vec1_t,
      ref_vec0_t, ref_vec1_t, mask, max_error);
}

int bit_vec_filter_no_flipping_sse1(char* read, char* ref, int length, int max_error) {
  //Get ready the bits
  sse3_convert2bit1(read, read_vec0_t, read_vec1_t);
  sse3_convert2bit1(ref, ref_vec0_t, ref_vec1_t);

  //Get the mask
  __m128i mask;
  if (length >= SSE_BASE_NUM1)
    mask = _mm_set1_epi8(0xff);
  else
    mask = _mm_load_si128( (__m128i *) (MASK_SSE_END1 + (length *
	    SSE_BYTE_NUM)));

  return bit_vec_filter_no_flipping_m128_sse1(read_vec0_t, read_vec1_t,
      ref_vec0_t, ref_vec1_t, mask, max_error);
}

int bit_vec_filter_sse11(char* read, char* ref, int length, int max_error) {
  //Get ready the bits
  //	memcpy(read_t, read, length * sizeof(char));
  //	memcpy(ref_t, ref, length * sizeof(char));

  sse3_convert2bit11(read, length, read_bit_t);
  sse3_convert2bit11(ref, length, ref_bit_t);

  return bit_vec_filter_m128_sse11(read_bit_t, ref_bit_t, length, max_error);
}

void bit_vec_filter_no_flipping_sse_simulate1(char* read, char* ref, int length,
    int max_error, int loc_num) {
  //Get ready the bits
  sse3_convert2bit1(read, read_vec0_t, read_vec1_t);
  sse3_convert2bit1(ref, ref_vec0_t, ref_vec1_t);

  //Get the mask
  __m128i mask;
  if (length >= SSE_BASE_NUM1)
    mask = _mm_set1_epi8(0xff);
  else
    mask = _mm_load_si128( (__m128i *) (MASK_SSE_END1 + (length *
	    SSE_BYTE_NUM)));

  while (loc_num--)
    bit_vec_filter_no_flipping_m128_sse1(read_vec0_t, read_vec1_t,
	ref_vec0_t, ref_vec1_t, mask, max_error);

  return;
}

void bit_vec_filter_sse_simulate1(char* read, char* ref, int length,
    int max_error, int loc_num) {
  //Get ready the bits
  sse3_convert2bit1(read, read_vec0_t, read_vec1_t);
  sse3_convert2bit1(ref, ref_vec0_t, ref_vec1_t);

  //Get the mask
  __m128i mask;
  if (length >= SSE_BASE_NUM1)
    mask = _mm_set1_epi8(0xff);
  else
    mask = _mm_load_si128( (__m128i *) (MASK_SSE_END1 + (length *
	    SSE_BYTE_NUM)));

  while (loc_num--)
    bit_vec_filter_m128_sse1(read_vec0_t, read_vec1_t,
	ref_vec0_t, ref_vec1_t, mask, max_error);

  return;
}

void bit_vec_filter_sse_simulate11(char* read, char* ref, int length,
    int max_error, int loc_num) {
  //Get ready the bits
  //	memcpy(read_t, read, length * sizeof(char));
  //	memcpy(ref_t, ref, length * sizeof(char));

  sse3_convert2bit11(read, length, read_bit_t);
  sse3_convert2bit11(ref, length, ref_bit_t);

  while (loc_num--)
    bit_vec_filter_m128_sse11(read_bit_t, ref_bit_t, length, max_error);
}

