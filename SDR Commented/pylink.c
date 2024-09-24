

// https://docs.python.org/3/extending/embedding.html

#include "pylink.h"
#include <stdio.h> 

// Runs an analysis by calling a Python function from C with the provided parameters.
int run_analysis(PARAM *param) {
    // Obtain the Python function "run_analysis" from the loaded Python module.
    PyObject *function = PyObject_GetAttrString(param->module, "run_analysis");

    // Check if the function exists and is callable.
    if (function && PyCallable_Check(function)) {

        // Prepare the arguments for the Python function call.
        PyObject *args = PyTuple_New(1); // Argument tuple with 1 element.
        PyObject *kmer_list = PyList_New(param->n_kmers); // Python list for k-mers.

        // Check if argument objects were created successfully.
        if (!args || !kmer_list) {
            printf("Cannot create analysis args\n");
            Py_Finalize(); // Clean up Python interpreter state.
            return -1; // Indicate failure.
        }

        // Populate the Python list with k-mers from the parameters.
        for(int n = 0; n < param->n_kmers; n++) {
            PyList_SetItem(kmer_list, n, Py_BuildValue("i", param->kmer_list[n]));
        }

        // Set the Python list as the argument for the Python function.
        PyTuple_SetItem(args, 0, kmer_list);

        // Call the Python function with the prepared arguments.
        PyObject *ret = PyObject_CallObject(function, args);
        Py_DECREF(args); // Decrease reference count for args to prevent memory leak.

        // Check the return value to determine if the call was successful.
        if (ret != NULL){
            Py_Finalize(); // Clean up Python interpreter state.
            return 0; // Indicate success.
        }
        else {
            PyErr_Print(); // Print Python error if any.
            printf("Analysis call failed\n");
            Py_Finalize(); // Clean up Python interpreter state.
            return -1; // Indicate failure.
        }
    }
    else {
        if (PyErr_Occurred()) PyErr_Print(); // Print Python error if any.
        printf("Cannot find analysis function\n");
        Py_Finalize(); // Clean up Python interpreter state.
        return -1; // Indicate failure.
    }
}

// Exports sequencing data to a Python environment for further processing.
int export_seqdata(PARAM *param, SEQDATA *seqdata) {
    // Obtain the Python function "import_seqdata" from the loaded Python module.
    PyObject *function = PyObject_GetAttrString(param->module, "import_seqdata");

    // Check if the function exists and is callable.
    if (function && PyCallable_Check(function)) {
        PyObject *args = PyTuple_New(7); // Prepare a tuple for 7 arguments.
        PyObject *dist_list = PyList_New(param->n_kmers); // List for k-mer distributions.

        // Check if argument objects were created successfully.
        if (!args || !dist_list) {
            printf("Cannot create export args\n");
            return -1; // Indicate failure.
        }

        // Populate the k-mer distributions list for each k-mer.
        for(uint32_t n = 0; n < param->n_kmers; n++) {
            PyObject *dist = PyList_New(0); // List for individual k-mer distribution.
            for(int i = 0; i < param->max_kmer_freq; i++){
                // Append k-mer frequency to the distribution list.
                PyList_Append(dist, Py_BuildValue("i", get_kpoint(seqdata, n, i)));
            }
            // Set the distribution list in the k-mer distributions list.
            PyList_SetItem(dist_list, n, dist);
        }

        // Open a file for appending the sequence data.
        FILE *fp;
        fp = fopen("/jlf/peadams/SIDR3_runs/DF5018_Osc/next_New/C_output_test.txt", "a");

        // Write various sequence data attributes to the file.
        fprintf(fp, "%s\t", get_name(seqdata));
        fprintf(fp, "%i\t", get_length(seqdata));
        fprintf(fp, "%f\t", get_gc(seqdata));
        fprintf(fp, "%f\t", get_cov(seqdata));
        fprintf(fp, "%i\t", get_blast(seqdata));
        fprintf(fp, "%i\t", get_tax(seqdata));
// Paula's print statements for the kmer_distribution. need to figure out exact formating for this.
	fprintf(fp,"["); // Start of kmer_distribution JSON array.
	for(uint32_t n = 0; n < param->n_kmers; n++) { // Iterate over each k-mer.
	    fprintf(fp, "["); // Start of this k-mer's frequency array.
            for(int i = 0; i < param->max_kmer_freq; i++){ // Iterate over frequencies of this k-mer.
		if((i+1) < param->max_kmer_freq) {
                    fprintf(fp,"%i, ", get_kpoint(seqdata, n, i)); // Print frequency, followed by comma for non-last elements.
                }
		else {
		   fprintf(fp, "%i", get_kpoint(seqdata, n, i)); // Print last frequency without a trailing comma.
                }
            }
	    if((n+1) < param->n_kmers) {
                fprintf(fp, "], "); // End of this k-mer's frequency array, with comma for non-last k-mers.
            }
            else {
		fprintf(fp, "]"); // End of the last k-mer's frequency array.
            }
        }
	fprintf(fp, "]\n"); // End of the kmer_distribution JSON array.
	
        fclose(fp); // Close the file pointer.

// The commented-out Python integration code is presumably for adding k-mer data to a Python list and creating a Python tuple with sequence data attributes.
// This could be part of a larger application where C code interacts with Python, possibly using Python's C API.

        PyTuple_SetItem(args, 0, Py_BuildValue("s", get_name(seqdata))); // Set the sequence name in a Python tuple.
        PyTuple_SetItem(args, 1, Py_BuildValue("I", get_length(seqdata))); // Set the sequence length in the tuple.
        PyTuple_SetItem(args, 2, Py_BuildValue("f", get_gc(seqdata))); // Set the GC content in the tuple.
        PyTuple_SetItem(args, 3, Py_BuildValue("f", get_cov(seqdata))); // Set the coverage in the tuple.
        PyTuple_SetItem(args, 4, Py_BuildValue("i", get_blast(seqdata))); // Set the BLAST hit presence in the tuple.
        PyTuple_SetItem(args, 5, Py_BuildValue("i", get_tax(seqdata))); // Set the taxonomic identification in the tuple.
        PyTuple_SetItem(args, 6, dist_list); // Set the k-mer distribution list in the tuple.

        PyObject *ret = PyObject_CallObject(function, args); // Call a Python function with the tuple as an argument.
	Py_DECREF(args); // Decrease the reference count for the Python tuple to prevent memory leaks.
        freeSEQDATA(seqdata, param->n_kmers); // Free the SEQDATA structure after use.

        if (ret != NULL){
            Py_DECREF(ret); // Decrease the reference count for the Python function return value.
//            printf("exporting data \n"); // Indicate successful data export.
            return 0; // Return success.
        }

        else {
            PyErr_Print(); // Print Python error if any.
            printf("Export call failed\n"); // Indicate a failure in the export call.
            return -1; // Return failure.
        }
    }

    else {
        if (PyErr_Occurred()) PyErr_Print(); // Check and print any Python errors if the function wasn't found.
        printf("Cannot find export function\n"); // Indicate the export function wasn't found.
        return -1; // Return failure.
    }
}