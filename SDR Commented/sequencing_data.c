#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "sequencing_data.h"

// Defines a structure for holding sequencing data.
struct sequencing_data {
    char *name;              // Name of the region.
    uint32_t length;         // Length of the region in bases.
    float gc;                // GC content percentage of the region.
    float cov;               // Read coverage of the region.
    bool blast;              // Flag indicating if the region has a BLAST hit.
    bool target;             // Flag indicating if the BLAST hit matches the target classification.
    uint32_t **k_dists;      // Pointer to an array of pointers for storing k-mer coverage distributions.
};

// Function to allocate and initialize a new sequencing data record.
SEQDATA *newSEQDATA(uint32_t n_kmers) {
    SEQDATA *seqdata = malloc(sizeof(SEQDATA)); // Allocate memory for the SEQDATA structure.
    assert(seqdata != 0); // Ensure memory allocation was successful.
    // Allocate memory for the k-mer distributions array.
    seqdata->k_dists = malloc(n_kmers * sizeof(uint32_t *));
    assert(seqdata->k_dists != 0); // Ensure memory allocation for k-mer distributions was successful.
    return seqdata; // Return the pointer to the newly allocated SEQDATA structure.
}

// Function to free the memory allocated for a sequencing data record.
void freeSEQDATA(SEQDATA *seqdata, uint32_t n_kmers) {
    if (seqdata == 0) return; // If the pointer is null, return immediately.
    if (seqdata->name) {
        // Free each array of k-mer distributions.
        for (uint32_t k = 0; k < n_kmers; k++) free(seqdata->k_dists[k]);
        free(seqdata->k_dists); // Free the array of pointers to k-mer distributions.
        free(seqdata->name); // Free the memory allocated for the name.
    }
    free(seqdata); // Finally, free the SEQDATA structure itself.
}

// Function to display the contents of a sequencing data record.
void displaySEQDATA(SEQDATA *seqdata, uint32_t n_kmers, uint32_t max_freq) {
    printf("Region: %s\n", seqdata->name); // Display the region name.
    printf("Length: %d\n", seqdata->length); // Display the region length.
    printf("GC Content: %f\n", seqdata->gc); // Display the GC content.
    printf("Read Coverage: %f\n", seqdata->cov); // Display the read coverage.
    printf("BLAST Hit: %s\n", seqdata->blast ? "true" : "false"); // Display if there's a BLAST hit.
    printf("Target taxon: %s\n", seqdata->target ? "true" : "false"); // Display if the hit matches target.
    printf("Kmer dists:\n"); // Display k-mer distributions.
    for(uint32_t i = 0; i < n_kmers; i++){
        for(uint32_t j = 0; j < max_freq; j++) printf("%lu ", seqdata->k_dists[i][j]);
        printf("\n");
    }
    printf("\n");
}

// Function to return the name of the sequencing data region.
char *get_name(SEQDATA *seqdata) {
    return seqdata->name; // Return the name of the region.
}

// Function to return the length of the sequencing data region.
uint32_t get_length(SEQDATA *seqdata) {
    return seqdata->length; // Return the length of the region.
}

// Function to return the coverage of the sequencing data region.
float get_cov(SEQDATA *seqdata) {
    return seqdata->cov; // Return the read coverage.
}

// Function to return the GC content of the sequencing data region.
float get_gc(SEQDATA *seqdata) {
    return seqdata->gc; // Return the GC content.
}

// Function to return the presence of a BLAST hit for the sequencing data region.
bool get_blast(SEQDATA *seqdata) {
    return seqdata->blast; // Return the BLAST hit flag.
}

// Function to return whether the BLAST hit matches the target classification.
bool get_tax(SEQDATA *seqdata) {
    return seqdata->target; // Return the target classification match flag.
}

// Retrieves a specific value from the k-mer coverage distribution for a given sequencing data.
// `k` is the index for the specific k-mer distribution array, and `i` is the index within that array.
extern uint32_t get_kpoint(SEQDATA *seqdata, uint32_t k, uint32_t i) {
    return seqdata->k_dists[k][i]; // Returns the specific k-mer coverage count.
}

// Updates the name of a sequencing data record. This is called by the process_taxonomy function in pipeline.c.
// The new name is dynamically allocated and copied into the SEQDATA structure.
void update_name(SEQDATA *seqdata, const char *name) {
    seqdata->name = malloc(strlen(name) + 1); // Allocates memory for the new name, +1 for the null terminator.
    strcpy(seqdata->name, name); // Copies the provided name into the SEQDATA structure.
}

// Updates the length property of a sequencing data record.
// The length is a measure of the nucleotide sequence covered by the sequencing data.
void update_length(SEQDATA *seqdata, uint32_t length) {
    seqdata->length = length; // Sets the new length of the sequencing data.
}

// Updates the GC content percentage of a sequencing data record.
// The GC content is calculated as the ratio of 'G' and 'C' bases to the total number of bases, multiplied by 100.
void update_gc(SEQDATA *seqdata, uint64_t gc_count) {
    assert(seqdata->length != 0); // Ensures that the sequence length is not zero to avoid division by zero.
    float gc = ((double)gc_count / (float)seqdata->length) * 100; // Calculates the GC content percentage.
    seqdata->gc = gc; // Sets the new GC content percentage.
}

// Updates the coverage property of a sequencing data record.
// Coverage is calculated as the ratio of base pairs sequenced to the total length of the sequence.
void update_cov(SEQDATA *seqdata, uint64_t bp_count) {
    assert(seqdata->length != 0); // Ensures that the sequence length is not zero to avoid division by zero.
    float cov = (double)bp_count / (float)seqdata->length; // Calculates the coverage.
    seqdata->cov = cov; // Sets the new coverage value.
}

// Updates the flag indicating whether a sequencing data record has a BLAST hit.
// This flag can be used to indicate the presence of significant matches found by the BLAST algorithm.
void update_blast(SEQDATA *seqdata, bool has_blast) {
    seqdata->blast = has_blast; // Sets the new BLAST hit flag.
}

// Updates the flag indicating whether the BLAST hit matches a target classification.
// This can be used to indicate if the identified sequences are of particular interest, based on classification criteria.
void update_tax(SEQDATA *seqdata, bool is_target) {
    seqdata->target = is_target; // Sets the new target classification match flag.
}

// Updates the k-mer distribution for a specific k-mer size within a sequencing data record.
// `k` is the index for the k-mer distribution to update, and `dist` is the pointer to the new distribution array.
void update_kdist(SEQDATA *seqdata, uint32_t k, uint32_t *dist) {
    seqdata->k_dists[k] = dist; // Updates the k-mer distribution array for the specified k-mer size.
}
