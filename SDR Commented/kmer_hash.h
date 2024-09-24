// Copyright 2020 Joshua J Baker. All rights reserved.
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file.

#ifndef HASHMAP_H
#define HASHMAP_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

// Forward declaration of the `kmer_hash` struct.
// This struct is defined in detail elsewhere (likely in a corresponding .c source file).
// It represents the structure of the hash table used to store and manage k-mers.
typedef struct kmer_hash HASH;

// Declaration of external functions that operate on the hash table.
// These functions must be implemented in the corresponding .c source file.

// Creates and initializes a new hash table with a specified minimum capacity
// and hash function seeds. Returns a pointer to the hash table.
extern HASH *newHASH(uint32_t, uint64_t, uint64_t);

// Clears the hash table, resetting it to its initial state.
extern void clearHASH(HASH *);

// Inserts a k-mer into the hash table.
extern void insertHASH(HASH *, uint64_t);

// Displays the contents of the hash table, typically used for debugging purposes.
extern void displayHASH(HASH *);

// Frees the hash table and releases any resources it was using.
extern void freeHASH(HASH *);

// Generates a distribution of k-mer frequencies within the hash table.
// Returns a pointer to an array containing the frequency distribution.
extern uint32_t *distHASH(HASH *, uint32_t);

// Declaration of the hash function used by the hash table.
// This function calculates a hash value based on a given k-mer and seed values.
extern uint64_t hash_func(const uint64_t *, size_t, uint64_t, uint64_t);

#endif

