
#ifndef param_h
#define param_h

#define PY_SSIZE_T_CLEAN
#include <Python.h>     // Include the Python API to allow embedding Python into the C program.
#include <stdint.h>     // Include the standard integer library for fixed-width integer types.

// Definition of the PARAM structure which holds configuration and file path information,
// along with Python objects necessary for the program's execution.
typedef struct {
    uint32_t n_threads;  // Number of threads for parallel processing.

    // Filepaths for input and taxonomy files.
    char *assembly;      // Path to the assembly file.
    char *alignment;     // Path to the alignment file.
    char *blast;         // Path to the BLAST query output file.
    char *delnodes;      // Path to the 'delnodes.dmp' file from NCBI Taxonomy database.
    char *merged;        // Path to the 'merged.dmp' file from NCBI Taxonomy database.
    char *nodes;         // Path to the 'nodes.dmp' file from NCBI Taxonomy database.
    char *names;         // Path to the 'names.dmp' file from NCBI Taxonomy database.

    // Taxonomy information.
    char *rank;              // Target taxonomic rank (e.g., "phylum").
    char *classification;    // Target classification to identify.

    // K-mer distribution information.
    uint32_t max_kmer_freq;  // Maximum frequency of k-mers to consider.
    uint32_t n_kmers;        // Number of k-mer sizes to analyze.
    uint32_t *kmer_list;     // List of k-mer sizes.

    PyObject *module;    // Python module for analysis, loaded at runtime.
} PARAM;

// Declarations of functions for initializing and freeing the PARAM structure.
// These functions are defined elsewhere, likely in a corresponding .c source file.
extern void initPARAM(PARAM *, int, char **);
extern void freePARAM(PARAM *);

#endif

