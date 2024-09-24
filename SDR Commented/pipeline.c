#include <pthread.h>

#include "htslib/thread_pool.h" // For using HTSlib's thread pool functionality.
#include "htslib/faidx.h" // For FASTA index handling in HTSlib.
#include "htslib/sam.h" // For SAM/BAM/CRAM operations in HTSlib.
#include "pipeline.h" // Custom header for pipeline operations.
#include "taxonomy.h" // Custom header for taxonomy-related operations.
#include "sequence.h" // Custom header for sequence encoding and processing.
#include "kmer_hash.h" // Custom header for k-mer hashing functionality.
#include "pylink.h" // Custom header for linking Python functionality, if applicable.

// Forward declarations of static functions that represent steps or tasks in the pipeline.
static void *pipe_in2tax(void *arg);
static void *process_taxonomy(void *arg);
static void *pipe_tax2seq(void *arg);
static void *process_sequence(void *arg);

// Structure representing options and shared data for the pipeline.
typedef struct {
    hts_tpool *p; // Pointer to a thread pool for managing concurrent tasks.
    hts_tpool_process *q1; // Queue for the first stage of the pipeline (e.g., taxonomy processing).
    hts_tpool_process *q2; // Queue for the second stage of the pipeline (e.g., sequence processing).
    PARAM *param; // Custom structure containing parameters for the pipeline.
    uint32_t n_seq; // Number of sequences to process.
    char **id_list; // List of sequence identifiers.
    TAX **blast_table; // Table containing BLAST results for sequences.
    TAX *check_table; // Table for tracking which sequences have been processed or checked.
    pthread_mutex_t *lock; // Mutex for synchronizing access to shared resources.
} pipe_opt;

// Structure for a job in the pipeline, encapsulating task-specific data.
typedef struct {
    pipe_opt *o; // Options and shared data for the pipeline.
    SEQDATA *seqdata; // Data structure for holding sequence information.
    uint32_t id; // Identifier for the current job or sequence.
} pipe_job;

// Entry point for the pipeline's first stage, initializing sequence data and dispatching taxonomy processing tasks.
// This function is responsible for creating jobs for processing taxonomy information and submitting them to the pipeline.
static void *pipe_in2tax(void *arg) {
    pipe_opt *o = (pipe_opt *)arg;

    // Iterate over all sequences, creating and dispatching jobs for taxonomy processing.
    for (uint32_t i = 0; i < o->n_seq; i++) {
        pipe_job *job = malloc(sizeof(pipe_job));
        assert(job != 0);
        job->o = o;
        job->seqdata = newSEQDATA(o->param->n_kmers); // Initialize sequence data for the job.
        job->id = i;

        // Dispatch the job for taxonomy processing. On failure, clean up and exit the thread.
        if (hts_tpool_dispatch(o->p, o->q1, process_taxonomy, job) != 0) {
            free(job);
            pthread_exit((void *)1);
        }
    }

    pthread_exit(0); // Exit the thread once all jobs have been dispatched.
}

// Processes taxonomy information for a given job, updating sequence data based on taxonomy results.
// This function performs taxonomy-related updates to sequence data and manages synchronization for shared resources.
static void *process_taxonomy(void *arg) {
    pipe_job *job = (pipe_job *)arg;

    // Update sequence name and free the corresponding identifier string.
    update_name(job->seqdata, job->o->id_list[job->id]);
    free(job->o->id_list[job->id]);

    // If there are BLAST results for the current sequence, process them.
    if(job->o->blast_table[job->id] != 0) {
        update_blast(job->seqdata, true); // Indicate that BLAST results are present.

        // Iterate over BLAST hits, performing taxonomy checks and updates.
        for(uint32_t i = 0; i < numHits(job->o->blast_table[job->id]); i++) {
            uint32_t curr_TaxID = getTaxID(job->o->blast_table[job->id], i); // Retrieve the current TaxID from BLAST results.

          // first check to see if current ID was already identified
          pthread_mutex_lock(job->o->lock);
          if(checkTAX(job->o->check_table, curr_TaxID)) { //calls checkTAX in taxonomy.c
              pthread_mutex_unlock(job->o->lock);
              update_tax(job->seqdata, true);
              break;
          }
// Unlocking the mutex before potentially long-running operations to avoid blocking other threads.
pthread_mutex_unlock(job->o->lock);

// If current ID is unidentified, then parse taxdump to verify and classify the ID.
uint32_t check = curr_TaxID;
// Continue to the next iteration if the current TaxID is marked as deleted.
if(check_delnodes(job->o->param->delnodes, check)) continue;
// Check and update the TaxID if it has been merged into another ID.
check_merged(job->o->param->merged, &check);
// Update the check ID based on the nodes database, adjusting for hierarchical changes.
check_nodes(job->o->param->nodes, &check, check, job->o->param->rank);
// If the TaxID matches a specified classification, proceed to update the check table and tax flag.
if(check_names(job->o->param->names, check, job->o->param->classification)) {
    // Locking the mutex to safely update shared resources.
    pthread_mutex_lock(job->o->lock);
    // Insert the current TaxID into the check table.
    insertTAX(job->o->check_table, curr_TaxID);
    // Unlock the mutex after the update is complete.
    pthread_mutex_unlock(job->o->lock);
    // Update the tax flag for the sequence data to true, indicating a match.
    update_tax(job->seqdata, true);
}

// Free the allocated memory for the blast table entry corresponding to the job ID.
freeTAX(job->o->blast_table[job->id]);

return (void *)job; // Return the job as void pointer, signaling completion.

// This function serves as a bridge in the pipeline, handling the transition from taxonomic analysis to sequence processing.
static void *pipe_tax2seq(void *arg) {
    pipe_opt *o = (pipe_opt *)arg; // Cast the argument back to the expected type.
    hts_tpool_result *result;

    // Continuously process results from the previous queue as long as they are available.
    while ((result = hts_tpool_next_result_wait(o->q1))) {
        pipe_job *job = (pipe_job *)hts_tpool_result_data(result);
        hts_tpool_delete_result(result, 0); // Free the result after retrieving the job.

        // Dispatch the job to the next stage in the pipeline; exit thread on failure.
        if (hts_tpool_dispatch(job->o->p, job->o->q2, process_sequence, job) != 0) pthread_exit((void *)1);

        // Check if the input queue is empty and cleanup if true.
        if (hts_tpool_process_empty(o->q1)) {
            free(o->id_list); // Free the ID list used for tracking job IDs.
            free(o->blast_table); // Free the blast table containing results for each ID.
            freeTAX(o->check_table); // Free the check table used for taxonomic verification.
            break; // Exit the loop as processing is complete.
        }
    }

    pthread_exit(0); // Exit the thread cleanly after processing all items.

}

// Processes each sequence, computing GC content, k-mer distribution, and read coverage.
static void *process_sequence(void *arg) {
    pipe_job *job = (pipe_job *)arg; // Cast the argument back to the expected type.

    // ========== GC CONTENT ========== //
    // Load the fasta index for the assembly associated with the job.
    faidx_t *fa_idx = fai_load(job->o->param->assembly);
    assert(fa_idx != 0); // Ensure the fasta index was successfully loaded.
    uint32_t seq_len; // Variable to hold the length of the fetched sequence.
    // Fetch the sequence based on the name retrieved from seqdata.
    char *seq = fai_fetch(fa_idx, get_name(job->seqdata), &seq_len);
    // Update the sequence length in the seqdata structure.
    update_length(job->seqdata, seq_len);

    // Calculate GC content and encode the sequence.
    uint64_t gc_count = 0; // Variable to hold the count of GC bases.
    SEQCODE *seqcode = encode(seq, seq_len, &gc_count); // Encode the sequence, counting GC bases.
    update_gc(job->seqdata, gc_count); // Update the GC content in seqdata.

    free(seq); // Free the fetched sequence string.
    fai_destroy(fa_idx); // Free the fasta index.

// ================ KMER DISTRIBUTION ================ //

// Initializes a new hash table for k-mer distribution analysis using initial seeds from the SEQCODE.
HASH *map = newHASH(seq_len, seed0(seqcode), seed1(seqcode));
// Iterates over the list of k-mer lengths defined in the job parameters.
for(uint32_t k = 0; k < job->o->param->n_kmers; k++) {
    uint32_t kmer_len = job->o->param->kmer_list[k];
    // For each k-mer length, iterates over the sequence to insert k-mers into the hash table.
    for(int i = 0; i <= (signed)(seq_len - kmer_len); i++) insertHASH(map, get_kmer(seqcode, kmer_len, i));
    // Updates the k-mer distribution in SEQDATA structure for the current k-mer length.
    update_kdist(job->seqdata, k, distHASH(map, job->o->param->max_kmer_freq));
    // Clears the hash table for the next iteration.
    clearHASH(map);
}
// Frees the hash table and the SEQCODE structure after completing k-mer distribution analysis.
freeHASH(map);
freeSEQCODE(seqcode);

// ================ READ COVERAGE ================ //

// Opens a SAM/BAM file for reading alignments.
samFile *sam = sam_open(job->o->param->alignment, "r");
assert(sam != 0); // Ensures the file is successfully opened.
// Loads the index file for the opened SAM/BAM file.
hts_idx_t *hts_idx = sam_index_load(sam, job->o->param->alignment);
assert(hts_idx != 0); // Ensures the index is successfully loaded.
// Reads the header from the SAM/BAM file.
bam_hdr_t *bam_hdr = sam_hdr_read(sam);
assert(bam_hdr != 0); // Ensures the header is successfully read.
// Ensures the sequence length matches the target length specified in the BAM header.
assert(seq_len == bam_hdr->target_len[job->id]);
// Initializes an iterator for reading over specified region or regions.
hts_itr_t *itr = sam_itr_querys(hts_idx, bam_hdr, get_name(job->seqdata));
assert(itr != 0); // Ensures the iterator is successfully initialized.
// Initializes a structure for storing alignment results.
bam1_t *itr_data = bam_init1();
assert(itr_data != 0); // Ensures the structure is successfully initialized.

// Counts the base pairs aligned within the sequence length to calculate coverage.
uint64_t bp_count = 0;
while (sam_itr_next(sam, itr, itr_data) >= 0) {
    bam1_core_t *c = &(itr_data->core);
    // Adds the aligned sequence length to the count, adjusting for sequences spilling over the end.
    if((c->pos + c->l_qseq) < seq_len) bp_count += c->l_qseq;
    else bp_count += (seq_len - c->pos);
}
// Updates the SEQDATA structure with the calculated read coverage.
update_cov(job->seqdata, bp_count);

// Cleans up allocated structures and closes the SAM/BAM file.
bam_destroy1(itr_data);
hts_itr_destroy(itr);
bam_hdr_destroy(bam_hdr);
hts_idx_destroy(hts_idx);
sam_close(sam);

// Exports the processed SEQDATA.
export_seqdata(job->o->param, job->seqdata);

// Frees the job structure after processing is complete.
free(job);


// Function that initializes and manages the pipeline process.
int run_pipeline(PARAM *param) {
    // Load quantity and names of sequences from a FASTA index file.
    faidx_t *fa_idx = fai_load(param->assembly);
    assert(fa_idx != 0); // Ensures the index is successfully loaded.
    uint32_t n_seq = faidx_nseq(fa_idx); // Gets the number of sequences.
    // Allocates memory for storing sequence names.
    char **id_list = malloc(n_seq * sizeof(char *));
    // Copies sequence names into the allocated array.
    for(uint32_t i = 0; i < n_seq; i++){
        const char *seq_name = faidx_iseq(fa_idx, i);
        id_list[i] = malloc(strlen(seq_name) + 1);
        strcpy(id_list[i], seq_name);
    }
    // Frees the FASTA index structure after use.
    fai_destroy(fa_idx);

    // Parses BLAST results and initializes taxonomy processing.
TAX **blast_table = parse_blast(param->blast, id_list, n_seq);
TAX *check_table = newTAX(16); // Initializes a new table with a capacity for 16 items.

// Calls functions defined in taxonomy.c to load BLAST hits into the blast_table and initializes a new check_table.
TAX **blast_table = parse_blast(param->blast, id_list, n_seq);
TAX *check_table = newTAX(16); // A fresh taxonomy table for checks is initialized again.
pthread_mutex_t lock;
pthread_mutex_init(&lock, NULL); // Initializes a mutex lock for thread synchronization.

// Initialize thread pool for concurrent processing.
hts_tpool *p = hts_tpool_init(param->n_threads); // Initializes a thread pool with a specified number of threads.
hts_tpool_process *q1 = hts_tpool_process_init(p, n_seq*2, 0);    // Initializes a process queue for taxonomy with double the number of sequences as its depth.
hts_tpool_process *q2 = hts_tpool_process_init(p, n_seq, 1);      // Initializes a process queue for sequences, with the sequence count as its depth and ordered processing enabled.
pipe_opt o = {p, q1, q2, param, n_seq, id_list, blast_table, check_table, &lock}; // Structures options for pipeline processing.

// Launch data source and sink threads.
// Two threads are created: one to pipe input to taxonomy processing, and another to pipe taxonomy to sequence processing.
pthread_t tidIto1, tid1to2; // Thread identifiers.
pthread_create(&tidIto1, NULL, pipe_in2tax, &o); // Creates a thread to process input to taxonomy.
pthread_create(&tid1to2, NULL, pipe_tax2seq, &o); // Creates a thread to process taxonomy to sequence.

// Wait for the pipeline threads to complete.
void *retv; // Pointer to capture return value from threads.
int ret = 0; // Variable to accumulate any error codes.
pthread_join(tidIto1, &retv); ret |= (retv != NULL); // Joins the first thread and checks for errors.
pthread_join(tid1to2, &retv); ret |= (retv != NULL); // Joins the second thread and checks for errors.

// Cleanup after thread pool processing.
pthread_mutex_destroy(&lock); // Destroys the mutex lock.
hts_tpool_process_destroy(q1); // Destroys the taxonomy processing queue.
hts_tpool_process_flush(q2); // Flushes the sequence processing queue to ensure all tasks are completed.
hts_tpool_process_destroy(q2); // Destroys the sequence processing queue.
hts_tpool_destroy(p); // Destroys the thread pool.

return ret; // Returns the accumulated error code.
}