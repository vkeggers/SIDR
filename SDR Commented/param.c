
#include <argp.h>
#include "param.h"

// Function to construct a full file path from a directory and a filename.
// It ensures the directory path ends with a slash before appending the filename.
static char *getFilepathFromDir(const char *dir, const char *filename) {
    int dlen = strlen(dir); // Length of the directory string.
    int flen = strlen(filename); // Length of the filename string.
    // Allocate enough memory for the full path plus a slash and null terminator.
    char *filepath = malloc(dlen + 1 + flen + 1);
    assert(filepath != 0); // Ensure the memory allocation was successful.
    // Construct the filepath. If the directory string does not end with a slash, add one.
    if(dir[dlen-1] == '/') strcat(strcpy(filepath, dir), filename);
    else strcat(strcat(strcpy(filepath, dir), "/"), filename);
    return filepath; // Return the dynamically allocated full file path.
}

// Commented-out function intended to parse k-mer related arguments.
// static void parse_kmers(struct argp_state *state, const char *arg, const char *info) {}

// Function to parse command-line options and arguments.
static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    PARAM *param = state->input; // Access the PARAM structure passed as input.

    // Switch based on the key (option) encountered.
    switch (key) {
        case 'f':
            param->assembly = arg; // Assembly file path.
            break;
        case 'a':
            param->alignment = arg; // Alignment file path.
            break;
        case 'b':
            param->blast = arg; // BLAST query output file path.
            break;
        case 't':
            // For NCBI Taxdump files, construct paths using the directory provided.
            param->delnodes = getFilepathFromDir(arg, "delnodes.dmp");
            param->merged = getFilepathFromDir(arg, "merged.dmp");
            param->nodes = getFilepathFromDir(arg, "nodes.dmp");
            param->names = getFilepathFromDir(arg, "names.dmp");
            break;
        case 'r':
            param->rank = arg; // Target taxonomic rank.
            break;
        case 'c':
            param->classification = arg; // Target classification.
            break;
        case 'n':
            param->n_threads = (uint32_t)strtol(arg, 0, 0); // Number of threads.
            break;
        case ARGP_KEY_ARGS:
            // Handle non-option arguments here if necessary.
            break;
        case ARGP_KEY_END:
            // Validate required options/arguments are provided.
            // argp_error is used to show an error and exit if something is missing.
            if (!param->assembly || !param->alignment || !param->blast ||
                !param->delnodes || !param->merged || !param->nodes || !param->names || !param->classification) {
                argp_error(state, "Missing required files or directories.");
            }
            break;
        default:
            return ARGP_ERR_UNKNOWN; // Return error for unknown options.
    }
    return 0; // Return success.
}

// Initializes parameters, reads in arguments, and calls parse_opt to handle them.
void initPARAM(PARAM *param, int argc, char **argv) {
    // Set default values for various parameters.
    param->rank = "phylum";
    param->max_kmer_freq = 256;
    param->n_kmers = 11;
    // Define a default list of k-mer sizes to analyze.
    uint32_t kmer_list[11] = {15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
    param->kmer_list = malloc(sizeof(kmer_list)); // Allocate memory for k-mer list.
    memcpy(param->kmer_list, kmer_list, sizeof(kmer_list)); // Copy the default list into the allocated space.
    param->n_threads = 4; // Default number of threads.

    // Embed Python initialization, setting the sys.path to include a specific directory.
    PyRun_SimpleString("import sys\nsys.path.insert(0, \"/path/to/directory/\")\n");
    PyObject *file = PyUnicode_DecodeFSDefault("analysis");
    param->module = PyImport_Import(file); // Import the Python module for analysis.
    assert(param->module != 0 && "Failed to load module\n"); // Ensure the module was loaded successfully.
    Py_DECREF(file); // Decrease the reference count for the file object to prevent memory leaks.

    //


    struct argp_option options[] = {

            {"fa", 'f', "FILE", 0, "Assembly file (fasta/fa/fna)", 0},
            {"aln", 'a', "FILE", 0, "Alignment file (bam/sam/cram)", 0},
            {"blast", 'b', "FILE", 0, "NCBI Blast+ query file", 0},
            {"tax", 't', "DIRECTORY", 0, "NCBI Taxdump directory", 0},
            {"rank", 'r', "STRING", 0, "Target taxonomic rank, default = 'phylum'", 0},
            {"class", 'c', "STRING", 0, "Target classification to identify", 0},

            {"kmers", 'k', "[a,c-m,p ...]", 0, "Kmers to include in analysis, default = [15-25]", 1},
            {"threads", 'n', "N", 0, "Number of processing threads, default = 4", 1},

            {0}

    };

    struct argp argp = {options, parse_opt, "ARGS...", "Sequence Identification"};
    argp_parse(&argp, argc, argv, 0, 0, param);

}

void freePARAM(PARAM *param) {
//    printf("inside freeparam param.c\n");

    free(param->delnodes);
    free(param->merged);
    free(param->nodes);
    free(param->names);
    free(param->kmer_list);

}
