// Prevents the pipeline header file from being included multiple times in a single compilation unit.
#ifndef pipeline_h
#define pipeline_h

// Includes the header file that defines the structure for parameters (PARAM), 
// which are likely configurations or settings used in the pipeline.
#include "param.h"

// Declares a function to run a processing pipeline, taking a pointer to a PARAM structure as input.
// This function presumably orchestrates a series of data processing or analysis steps based on the parameters provided.
// The return type is int, which could indicate success/failure or provide a specific exit code.
extern int run_pipeline(PARAM *);

#endif // End of the include guard for pipeline_h
