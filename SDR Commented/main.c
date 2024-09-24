//
//  main.c
//  SIDR 2.0
//
//  Created by Adam Case on 01/31/2021.
//  Copyright © 2021 Adam Case.
//  All rights reserved.
//

#include "param.h"       // Include the header for parameter handling functions.
#include "pipeline.h"    // Include the header for pipeline processing functions.
#include "pylink.h"      // Include the header for Python integration utilities.

// Entry point of the program.
int main(int argc, char *argv[]) {
    printf("starting\n"); // Notification of program start.
    Py_Initialize();       // Initializes the Python interpreter.
    printf("after initializing\n"); // Notification that Python has been initialized.

    // Initialize a PARAM structure to store parameters.
    // This likely involves parsing command-line arguments and preparing them for use.
    PARAM param;
    initPARAM(&param, argc, argv); // Calls initPARAM from param.c to initialize the PARAM structure.
    printf("after initPARAM main.ci\n");

    int ret = 0; // Variable to track the success/failure status of function calls.

    // Pipeline processing part.
    // This involves running the main data processing workflow as defined in pipeline.c.
    ret |= run_pipeline(&param); // Processes data using the pipeline, updates 'ret' with the result.
    printf("after run_pipeline, main.c \n");

    // Analysis part.
    // Further data analysis may be performed after the initial pipeline processing.
    ret |= run_analysis(&param); // Performs additional analysis, updates 'ret' with the result.
    printf("after run_analysis, main.c \n");

    // Cleanup
    // Frees any resources allocated to 'param'.
    freePARAM(&param); // Calls freePARAM from param.c to clean up.
    printf("after freePARAM main.c\n");

    return ret; // Returns the status code. '0' typically indicates success.
}

