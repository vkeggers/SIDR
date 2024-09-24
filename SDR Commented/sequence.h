// This directive ensures the header is only included once in a single compilation.
#ifndef sequence_h
#define sequence_h

#include <stdint.h> // For fixed-width integers.
#include <stdbool.h> // For boolean type.

// Defines a structure for sequence encoding. The exact details of the structure are not specified here.
typedef struct sequence_encoding SEQCODE;

// Function prototypes for working with sequence encoding:

// Encodes a given character sequence into a SEQCODE structure, using specified parameters for the encoding process.
extern SEQCODE *encode(const char *, uint32_t, uint64_t *);

// Inserts a specified 64-bit integer (possibly representing a k-mer or similar) into the SEQCODE structure.
extern void insertSEQCODE(SEQCODE *, uint64_t);

// Returns the initial seed (possibly the first encoded k-mer or hash value) from a SEQCODE structure.
extern uint64_t seed0(SEQCODE *);

// Returns the second seed (possibly an alternative or subsequent encoded k-mer or hash value) from a SEQCODE structure.
extern uint64_t seed1(SEQCODE *);

// Frees the memory allocated to a SEQCODE structure, cleaning up resources.
extern void freeSEQCODE(SEQCODE *);

// Retrieves a specific k-mer or encoded segment from the SEQCODE structure, based on given indices or parameters.
extern uint64_t get_kmer(SEQCODE *, uint32_t, uint32_t);

// A utility function that prints the bits of a 64-bit integer, likely for debugging or visualization purposes.
extern void printBits(uint64_t);

#endif // End of include guard for sequence_h