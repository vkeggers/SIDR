// Prevents this header file from being included more than once in a single compilation.
#ifndef sequencing_data_h
#define sequencing_data_h

#include <stdbool.h> // Includes the standard boolean library for using 'bool' type.
#include <stdint.h>  // Includes the standard integer library for using fixed-width integers.

// Declares a new structure type for storing sequencing data.
typedef struct sequencing_data SEQDATA;

// Function declarations related to SEQDATA:

// Creates a new sequencing data record with specified initial capacity or parameter.
extern SEQDATA *newSEQDATA(uint32_t);

// Frees the memory allocated for a sequencing data record.
extern void freeSEQDATA(SEQDATA *, uint32_t);

// Displays information contained in a sequencing data record, potentially filtered or formatted based on the parameters.
extern void displaySEQDATA(SEQDATA *, uint32_t, uint32_t);

// Retrieves the name associated with a sequencing data record.
extern char *get_name(SEQDATA *);

// Retrieves the length (size) of the sequencing data.
extern uint32_t get_length(SEQDATA *);

// Retrieves the GC-content (percentage of guanine and cytosine) of the sequencing data.
extern float get_gc(SEQDATA *);

// Retrieves the coverage information of the sequencing data.
extern float get_cov(SEQDATA *);

// Retrieves a boolean indicating whether BLAST results are available for the sequencing data.
extern bool get_blast(SEQDATA *);

// Retrieves a boolean indicating whether taxonomic classification is available for the sequencing data.
extern bool get_tax(SEQDATA *);

// Retrieves a specific k-point (part of k-mer analysis) from the sequencing data, given indexes or parameters.
extern uint32_t get_kpoint(SEQDATA *, uint32_t, uint32_t);

// Updates the name associated with a sequencing data record.
extern void update_name(SEQDATA *, const char *);

// Updates the length (size) of the sequencing data.
extern void update_length(SEQDATA *, uint32_t);

// Updates the coverage information of the sequencing data.
extern void update_cov(SEQDATA *, uint64_t);

// Updates the GC-content (percentage of guanine and cytosine) of the sequencing data.
extern void update_gc(SEQDATA *, uint64_t);

// Updates the boolean indicating whether BLAST results are available for the sequencing data.
extern void update_blast(SEQDATA *, bool);

// Updates the boolean indicating whether taxonomic classification is available for the sequencing data.
extern void update_tax(SEQDATA *, bool);

// Updates the k-mer distribution or specific k-points in the sequencing data, given a set of values or parameters.
extern void update_kdist(SEQDATA *, uint32_t, uint32_t *);

#endif // End of include guard for sequencing_data_h
