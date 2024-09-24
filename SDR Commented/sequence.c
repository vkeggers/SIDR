#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>

#include "sequence.h"

// Structure definition for encoding a sequence.
struct sequence_encoding {
    uint32_t capacity; // How many 64-bit integers the code can hold.
    uint64_t *code;    // Pointer to an array of 64-bit integers storing the encoded sequence.
};

// Function to create a new sequence encoding object.
// `seq_len` is the length of the sequence to be encoded.
static SEQCODE *newSEQCODE(uint32_t seq_len) {
  SEQCODE *seqcode = malloc(sizeof(SEQCODE)); // Allocate memory for the sequence encoding structure.
  assert(seqcode != 0); // Ensure memory allocation was successful.
  // Calculate how many 64-bit integers are needed to store the sequence.
  seqcode->capacity = seq_len / 32;
  if(seq_len % 32 != 0) seqcode->capacity++; // Handle sequences not perfectly divisible by 32.
  // Allocate memory for the encoded sequence.
  seqcode->code = malloc(seqcode->capacity * sizeof(uint64_t));
  assert(seqcode->code != 0); // Ensure memory allocation was successful.
  return seqcode;
}

// Function to return the first 64-bit integer from the encoded sequence.
uint64_t seed0(SEQCODE *seqcode) {
    return seqcode->code[0]; // Returns the first code piece as seed0.
}

// Function to return the bitwise NOT of the first 64-bit integer from the encoded sequence.
// This could be used as an alternative seed or for hashing purposes.
uint64_t seed1(SEQCODE *seqcode) {
    return ~seqcode->code[0]; // Returns the bitwise NOT of the first code piece as seed1.
}

// Function to free the memory allocated for a sequence encoding structure.
void freeSEQCODE(SEQCODE *seqcode) {
    free(seqcode->code); // Free the encoded sequence array.
    free(seqcode); // Free the structure itself.
}

// Function to encode a given sequence into binary representation.
// `seq` is the input character sequence, `seq_len` is its length, and `gc_count` tracks the count of 'G' and 'C' characters.
SEQCODE *encode(const char *seq, uint32_t seq_len, uint64_t *gc_count) {
    SEQCODE *seqcode = newSEQCODE(seq_len); // Initialize a new sequence encoding structure.

    for(uint32_t i = 0; i < seqcode->capacity; i++) {
        uint64_t code = 0; // Initialize the code for this segment of the sequence.
        char seq_buff[32]; // Buffer to hold a segment of the sequence for encoding.
        uint32_t buff_len; // Length of the segment to encode.

        if(i < seqcode->capacity - 1) buff_len = 32; // Full segment except for possibly the last one.
        else buff_len = seq_len % 32; // Length of the last segment if it's shorter than 32.

        // Copy the current segment of the sequence into the buffer.
        memcpy(seq_buff, &seq[i * 32], buff_len);

        // Encode each character of the segment into binary.
        for(uint32_t j = 0; j < buff_len; j++) {
            switch(seq_buff[j]) {
                  case 'C':
                  case 'c':
                      code |= ((uint64_t)0x1 << (j * 2)); // Encode 'C' or 'c'.
                      (*gc_count)++; // Increment GC count.
                      break;
                  case 'G':
                  case 'g':
                      code |= ((uint64_t)0x2 << (j * 2)); // Encode 'G' or 'g'.
                      (*gc_count)++; // Increment GC count.
                      break;
                  case 'T':
                  case 't':
                      code |= ((uint64_t)0x3 << (j * 2)); // Encode 'T' or 't'.
                      break;
            }
        }

        // Store the encoded binary segment in the sequence encoding structure.
        seqcode->code[i] = code;
    }

    return seqcode; // Return the populated sequence encoding structure.
}

// Retrieves a specific k-mer from a sequence encoding structure, SEQCODE, using zero-based indexing.
// The function calculates the k-mer based on its position (`pos`) and length (`k`) within the sequence.
// It handles sequences encoded across multiple 64-bit integers and manages edge cases where a k-mer spans these boundaries.
uint64_t get_kmer(SEQCODE *seqcode, uint32_t k, uint32_t pos) {
    uint64_t kmer = 0; // Initialize the k-mer variable.
    uint32_t div = pos / 32; // Calculate which 64-bit integer in the array contains the start of the k-mer.
    uint32_t mod = pos % 32; // Calculate the position within the 64-bit integer where the k-mer starts.
    uint64_t buffer = seqcode->code[div]; // Retrieve the 64-bit integer containing the start of the k-mer.

    // If the k-mer spans across two 64-bit integers in the array:
    if (32 - mod < k) {
        // Extract the part of the k-mer from the first integer.
        kmer |= (buffer >> (mod * 2));
        // Handle the overflow into the next 64-bit integer.
        uint64_t overflow = seqcode->code[div + 1];
        // Combine the parts of the k-mer from both integers.
        kmer |= (overflow << ((64 - mod - k) * 2) >> ((32 - k) * 2));
    }
    // If the k-mer exactly fits the remaining part of the 64-bit integer:
    else if (32 - mod == k) {
        kmer |= (buffer >> (mod * 2));
    }
    // If the k-mer is fully contained within the current 64-bit integer:
    else {
        kmer |= (buffer << ((32 - mod - k) * 2) >> ((32 - k) * 2));
    }

    return kmer; // Return the extracted k-mer.
}

// Prints the binary representation of a 64-bit integer (`n`).
// This function iterates through each bit of `n`, from the most significant bit to the least, and prints '1' or '0'.
void printBits(uint64_t n) {
    // Start with the most significant bit.
    uint64_t i = 1UL << (sizeof(n) * CHAR_BIT - 1);
    // Loop until all bits are processed.
    while (i > 0) {
        // Check if the current bit is set and print '1' or '0'.
        if (n & i)
            printf("1");
        else
            printf("0");
        // Move to the next bit.
        i >>= 1;
    }
    printf("\n\n"); // Newline for readability.
}
