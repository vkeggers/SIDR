// based on https://github.com/tidwall/hashmap.c

Sure, I'll consolidate the explanation for the entire script, providing a comprehensive overview with comments detailing each significant portion:

```c
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "kmer_hash.h"

// Structure for storing information about a single data entry in the hash table.
typedef struct data_bucket {    
    uint64_t hash:48;   // The 48-bit hash value of the k-mer. Uses a bit-field to compactly store the hash.
    uint16_t dib;       // Distance from the initial bucket (DIB) for linear probing.
    uint64_t kmer;      // The k-mer value itself, assuming it fits in 64 bits.
    uint16_t count;     // A counter for the number of occurrences of this k-mer.
} BUCKET;

// Structure for the hash table, including metadata and pointers to the data buckets.
struct kmer_hash {
    uint64_t seed0;     // First seed for the hash function.
    uint64_t seed1;     // Second seed for the hash function.
    size_t bucketsz;    // The size of each bucket in bytes.
    size_t nbuckets;    // The total number of buckets in the hash table.
    size_t count;       // The count of unique k-mers stored in the hash table.
    size_t mask;        // A bitmask used for index calculations.
    BUCKET *buckets;    // Pointer to the array of buckets.
    BUCKET *temp;       // Temporary bucket for use during insertions.
    BUCKET *entry;      // Temporary bucket for the new entry to be inserted.
};

// Function to initialize and return a new hash table.
HASH *newHASH(uint32_t min, uint64_t seed0, uint64_t seed1) {
    int cap = 256; // Start with a default capacity.
    // Increase capacity to at least 'min', doubling as needed.
    if(cap < min) {
        while (cap < min) cap *= 2;
    }
    
    // Ensure the bucket size is aligned.
    size_t bucketsz = sizeof(BUCKET);
    while (bucketsz & (sizeof(uintptr_t) - 1)) bucketsz++;
    
    // Allocate memory for the hash table structure, including space for 'temp' and 'entry' buckets.
    size_t size = sizeof(HASH) + (bucketsz * 2);
    HASH *map = malloc(size);
    assert(map != 0); // Ensure the allocation was successful.
    memset(map, 0, sizeof(HASH)); // Initialize the allocated memory.

    // Initialize the hash table properties.
    map->bucketsz = bucketsz;
    map->seed0 = seed0;
    map->seed1 = seed1;
    map->temp = (BUCKET *)((char *)map + sizeof(HASH));
    map->entry = (BUCKET *)((char *)map->temp + bucketsz);
    map->nbuckets = cap;
    map->mask = map->nbuckets - 1;
    // Allocate memory for the buckets and initialize them.
    map->buckets = malloc(map->bucketsz * map->nbuckets);
    assert(map->buckets != 0); // Ensure the allocation was successful.
    memset(map->buckets, 0, map->bucketsz * map->nbuckets);

    return map; // Return the pointer to the newly created hash table.
}

// Function to reset the hash table, clearing all stored k-mers.
void clearHASH(HASH *map) {
    map->count = 0; // Reset the count of stored k-mers.
    memset(map->buckets, 0, map->bucketsz * map->nbuckets); // Clear the buckets.
}

// Function to insert a new k-mer into the hash table.
void insertHASH(HASH *map, uint64_t kmer) {
    // Prepare the new entry.
    BUCKET *entry = map->entry;
    entry->hash = hash_func(&kmer, sizeof(uint64_t), map->seed0, map->seed1);
    entry->dib = 1; // Set the initial distance from the bucket to 1.
    entry->kmer = kmer; // Store the k-mer.
    entry->count = 1; // Initialize the count to 1.

    // Find the correct bucket for the entry.
    size_t i = entry->hash & map->mask;
    for (;;) {
        BUCKET *bucket = &map->buckets[i];
        
        // If the bucket is empty, insert the entry.
        if (bucket->dib == 0) {
            memcpy(bucket, entry, map->bucketsz);
            map->count++; // Increment the count of unique k-mers.
            break;
        }
        
        // If the bucket contains the same k-mer, increment its count.
        if (entry->hash == bucket->hash && entry->kmer == bucket->kmer){
            bucket->count++;
            break;
        }
        if (bucket->dib < entry->dib) {
            memcpy(map->temp, bucket, map->bucketsz);
            memcpy(bucket, entry, map->bucketsz);
            memcpy(entry, map->temp, map->bucketsz);
  		  }

        entry->dib++;
		    i = (i + 1) & map->mask;

	  }

}



void displayHASH(HASH *map) {
//    printf("inside displayHASH\n");

    for(size_t i = 0; i < map->nbuckets; i++) {
        BUCKET *bucket = &map->buckets[i];
        printf("[%lu, %lu]\n", bucket->kmer, bucket->count);
    }

}



void freeHASH(HASH *map) {
//    printf("inside freehash kmerhash.c\n");

    if (!map) return;
    free(map->buckets);
    free(map);

}



uint32_t *distHASH(HASH *map, uint32_t max_freq){
//    printf("inside distHASH kmer_hash.c\n");

    uint32_t *dist = malloc(max_freq * sizeof(uint32_t));
    memset(dist, 0, max_freq * sizeof(uint32_t));
    if(map->count == 0) return dist;
    for(size_t i = 0; i < map->nbuckets; i++) {
      BUCKET *bucket = &map->buckets[i];
      if (bucket->dib) {
          if (bucket->count < max_freq) dist[bucket->count]++;
          else dist[max_freq-1]++;
      }
    }

    return dist;

}



//-----------------------------------------------------------------------------
// SipHash reference C implementation
//
// Copyright (c) 2012-2016 Jean-Philippe Aumasson
// <jeanphilippe.aumasson@gmail.com>
// Copyright (c) 2012-2014 Daniel J. Bernstein <djb@cr.yp.to>
//
// To the extent possible under law, the author(s) have dedicated all copyright
// and related and neighboring rights to this software to the public domain
// worldwide. This software is distributed without any warranty.
//
// You should have received a copy of the CC0 Public Domain Dedication along
// with this software. If not, see
// <http://creativecommons.org/publicdomain/zero/1.0/>.
//
// default: SipHash-2-4
//-----------------------------------------------------------------------------
static uint64_t SIP64(const uint8_t *in, const size_t inlen, uint64_t seed0, uint64_t seed1) {
 //    printf("inside uint64_t kmer_hash.c\n");
#define U8TO64_LE(p) \
    {  (((uint64_t)((p)[0])) | ((uint64_t)((p)[1]) << 8) | \
        ((uint64_t)((p)[2]) << 16) | ((uint64_t)((p)[3]) << 24) | \
        ((uint64_t)((p)[4]) << 32) | ((uint64_t)((p)[5]) << 40) | \
        ((uint64_t)((p)[6]) << 48) | ((uint64_t)((p)[7]) << 56)) }
#define U64TO8_LE(p, v) \
    { U32TO8_LE((p), (uint32_t)((v))); \
      U32TO8_LE((p) + 4, (uint32_t)((v) >> 32)); }
#define U32TO8_LE(p, v) \
    { (p)[0] = (uint8_t)((v)); \
      (p)[1] = (uint8_t)((v) >> 8); \
      (p)[2] = (uint8_t)((v) >> 16); \
      (p)[3] = (uint8_t)((v) >> 24); }
#define ROTL(x, b) (uint64_t)(((x) << (b)) | ((x) >> (64 - (b))))
#define SIPROUND \
    { v0 += v1; v1 = ROTL(v1, 13); \
      v1 ^= v0; v0 = ROTL(v0, 32); \
      v2 += v3; v3 = ROTL(v3, 16); \
      v3 ^= v2; \
      v0 += v3; v3 = ROTL(v3, 21); \
      v3 ^= v0; \
      v2 += v1; v1 = ROTL(v1, 17); \
      v1 ^= v2; v2 = ROTL(v2, 32); }
    uint64_t k0 = U8TO64_LE((uint8_t*)&seed0);
    uint64_t k1 = U8TO64_LE((uint8_t*)&seed1);
    uint64_t v3 = UINT64_C(0x7465646279746573) ^ k1;
    uint64_t v2 = UINT64_C(0x6c7967656e657261) ^ k0;
    uint64_t v1 = UINT64_C(0x646f72616e646f6d) ^ k1;
    uint64_t v0 = UINT64_C(0x736f6d6570736575) ^ k0;
    const uint8_t *end = in + inlen - (inlen % sizeof(uint64_t));
    for (; in != end; in += 8) {
        uint64_t m = U8TO64_LE(in);
        v3 ^= m;
        SIPROUND; SIPROUND;
        v0 ^= m;
    }
    const int left = inlen & 7;
    uint64_t b = ((uint64_t)inlen) << 56;
    switch (left) {
    case 7: b |= ((uint64_t)in[6]) << 48;
    case 6: b |= ((uint64_t)in[5]) << 40;
    case 5: b |= ((uint64_t)in[4]) << 32;
    case 4: b |= ((uint64_t)in[3]) << 24;
    case 3: b |= ((uint64_t)in[2]) << 16;
    case 2: b |= ((uint64_t)in[1]) << 8;
    case 1: b |= ((uint64_t)in[0]); break;
    case 0: break;
    }
    v3 ^= b;
    SIPROUND; SIPROUND;
    v0 ^= b;
    v2 ^= 0xff;
    SIPROUND; SIPROUND; SIPROUND; SIPROUND;
    b = v0 ^ v1 ^ v2 ^ v3;
    uint64_t out = 0;
    U64TO8_LE((uint8_t*)&out, b);
    return out;
}



uint64_t hash_func(const uint64_t *item, size_t len, uint64_t seed0, uint64_t seed1){
//    printf("inside hash_func kmerhash.c\n");

    return (SIP64((uint8_t *)item, len, seed0, seed1) << 16 >> 16);

}
