// This header file defines the interface for linking Python code with C code, specifically
// for running analyses and exporting sequencing data. It acts as a bridge between the two languages,
// allowing for efficient data processing and manipulation in C to be utilized from Python.

#ifndef pylink_h // Prevents the header file from being included multiple times.
#define pylink_h

#include "param.h" // Includes the definition of the PARAM structure for analysis parameters.
#include "sequencing_data.h" // Includes the definition of the SEQDATA structure for sequencing data.

// Function prototypes:

// Runs an analysis based on parameters specified in a PARAM structure.
// Returns an integer status code (e.g., success or error code).
extern int run_analysis(PARAM *);

// Exports sequencing data stored in a SEQDATA structure based on parameters specified in a PARAM structure.
// Returns an integer status code indicating success or failure of the operation.
extern int export_seqdata(PARAM *, SEQDATA *);

#endif // End of the include guard.
