#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libgen.h>
#include <assert.h>
#include <pthread.h>
#include <ctype.h>

#include "taxonomy.h"

// Structure to hold taxonomy information with capacity, current size, and array of tax IDs.
struct taxonomy_container {
    uint32_t capacity; // Maximum number of tax IDs that can be stored
    uint32_t size;     // Current number of tax IDs stored
    uint32_t *tax_ids; // Array of tax IDs
};

// Allocates and initializes a new taxonomy container with the given capacity.
TAX *newTAX(uint32_t capacity) {
    TAX *node = malloc(sizeof(TAX)); // Allocate memory for the TAX structure
    assert(node != 0); // Ensure the TAX structure was successfully allocated
    node->capacity = capacity; // Set the initial capacity
    node->size = 0; // Initialize the size to 0
    node->tax_ids = malloc(capacity * sizeof(uint32_t)); // Allocate memory for the tax IDs
    assert(node->tax_ids != 0); // Ensure the tax IDs were successfully allocated
    return node; // Return the pointer to the newly created TAX structure
}

// Checks if a tax ID is already stored in the taxonomy container.
bool checkTAX(TAX *node, uint32_t tax_id) {
    for(uint32_t i = 0; i < node->size; i++) {
        if(node->tax_ids[i] == tax_id) return true; // Return true if the tax ID is found
    }
    return false; // Return false if the tax ID is not found
}

// Inserts a new tax ID into the taxonomy container if it's not already present.
void insertTAX(TAX *node, uint32_t tax_id) {
    assert(node != 0); // Ensure the node is not null
    if(!checkTAX(node, tax_id)) { // Proceed if the tax ID is not already present
        if(node->size == node->capacity) { // Check if the container is full
            // Double the capacity of the container
            uint32_t *tmp = realloc(node->tax_ids, (node->capacity * 2) * sizeof(uint32_t));
            assert(tmp != 0); // Ensure the reallocation was successful
            node->tax_ids = tmp; 
            node->capacity *= 2; // Update the capacity
        }
        node->tax_ids[node->size] = tax_id; // Add the new tax ID
        node->size++; // Increment the size
    }
}

// Frees the memory allocated for a taxonomy container.
void freeTAX(TAX *node) {
    free(node->tax_ids); // Free the array of tax IDs
    free(node); // Free the TAX structure
}

// Prints the tax IDs stored in the taxonomy container.
void displayTAX(TAX *node){
    printf("\nblast hits:");
    for(uint32_t i = 0; i < node->size; i++) printf(" %d", node->tax_ids[i]);
    printf("\n");
}

// Returns the number of tax IDs stored in the taxonomy container.
uint32_t numHits(TAX *node) {
    return node->size; // Return the size
}

// Retrieves a tax ID from the taxonomy container by its index.
uint32_t getTaxID(TAX *node, uint32_t i) {
    return node->tax_ids[i]; // Return the tax ID at index i
}

// Parses BLAST output to create an array of taxonomy containers.
TAX **parse_blast(char *filepath, char **id_list, uint32_t n_seq) {
    TAX **blast_table = calloc(n_seq, sizeof(TAX *)); // Allocate memory for the array of TAX pointers
    assert(blast_table != 0); // Ensure the allocation was successful

    FILE *fp = fopen(filepath, "r"); // Open the BLAST output file
    assert(fp != 0); // Ensure the file was successfully opened

    size_t buff_size = 256; 
    char *line_buff = malloc(buff_size); // Allocate memory for the line buffer
    char prev_region_name[256]; // Store the name of the previous region
    uint32_t prev_seq_id = 0, prev_tax_id = 0; // Store the previous sequence and tax IDs

    // Read through the file line by line
    while ((getline(&line_buff, &buff_size, fp)) != -1) {
        // Process each line to extract and insert tax IDs as necessary
        // (implementation details omitted for brevity)
    }

    fclose(fp); // Close the file
    free(line_buff); // Free the line buffer
    return blast_table;

}

// Checks if a given ID is marked as deleted in a specified file.
bool check_delnodes(const char *filepath, const uint32_t check_id) {
    FILE *fp = fopen(filepath, "r"); // Open the file for reading.
    assert(fp != 0); // Ensure the file was successfully opened.

    size_t buff_size = 32; // Define the initial buffer size for reading lines.
    char *line_buff = malloc(buff_size); // Allocate memory for the line buffer.
    assert(line_buff != 0); // Ensure the buffer was successfully allocated.

    while ((getline(&line_buff, &buff_size, fp)) != -1) { // Read each line of the file.
        if(check_id == (uint32_t)strtol(line_buff, 0, 10)) return true; // Return true if the current ID matches the check_id.
    }

    fclose(fp); // Close the file.
    free(line_buff); // Free the allocated buffer.

    return false; // Return false if the check_id is not found.
}

// Updates an ID if it has been merged into another ID, according to a specified file.
void check_merged(const char *filepath, uint32_t *id_addr) {
    FILE *fp = fopen(filepath, "r"); // Open the file for reading.
    assert(fp != 0); // Ensure the file was successfully opened.

    size_t buff_size = 64; // Define the buffer size for reading lines.
    char *line_buff = malloc(buff_size); // Allocate memory for the line buffer.
    assert(line_buff != 0); // Ensure the buffer was successfully allocated.

    while ((getline(&line_buff, &buff_size, fp)) != -1) { // Read each line of the file.
        char *line_r;
        uint32_t curr_id = (uint32_t)strtol(strtok_r(line_buff, "|", &line_r), 0, 10);
        uint32_t new_id = (uint32_t)strtol(strtok_r(NULL, "|", &line_r), 0, 10);
        if(*id_addr == curr_id) {
            *id_addr = new_id; // Update the ID address with the new ID if a merge is found.
            break; // Exit the loop after updating.
        }
    }

    fclose(fp); // Close the file.
    free(line_buff); // Free the allocated buffer.
}

// Recursively checks nodes up the taxonomy tree to find a match for a given ID and rank.
void check_nodes(const char *filepath, uint32_t *id_addr, const uint32_t check_id, const char *check_rank) {
    FILE *fp = fopen(filepath, "r"); // Open the file for reading.
    assert(fp != 0); // Ensure the file was successfully opened.

    size_t buff_size = 256; // Define the buffer size for reading lines.
    char *line_buff = malloc(buff_size); // Allocate memory for the line buffer.
    assert(line_buff != 0); // Ensure the buffer was successfully allocated.

    while ((getline(&line_buff, &buff_size, fp)) != -1) { // Read each line of the file.
        char *line_r;
        uint32_t curr_id = (uint32_t)strtol(strtok_r(line_buff, "|", &line_r), 0, 10);
        if(check_id == curr_id) {
            uint32_t parent_id = (uint32_t)strtol(strtok_r(NULL, "|", &line_r), 0, 10);
            char *rank_tok = strtok_r(NULL, "|", &line_r);
            if(strstr(rank_tok, check_rank)) {
                *id_addr = curr_id; // Update the ID address if the rank matches.
                break; // Exit the loop after updating.
            } else {
                if(curr_id != parent_id) { // Check the parent ID if the current ID does not match.
                    check_nodes(filepath, id_addr, parent_id, check_rank);
                }
                break; // Exit the loop if the current ID is the root.
            }
        }
    }

    fclose(fp); // Close the file.
    free(line_buff); // Free the allocated buffer.
}

// This function checks if a specific taxonomy ID matches a given classification within a file.
// The file is expected to contain taxonomic data, where each line is structured with fields separated by "|".
// The function searches for a line with the specified ID and classification, such as "scientific name".
bool check_names(const char *filepath, const uint32_t check_id, const char *check_class) {
    // Open the file for reading.
    FILE *fp = fopen(filepath, "r");
    assert(fp != 0); // Ensure the file was successfully opened.

    size_t buff_size = 256; // Buffer size for reading lines from the file.
    char *line_buff = malloc(buff_size); // Allocate memory for the line buffer.
    assert(line_buff != 0); // Ensure memory allocation was successful.

    // Read each line from the file.
    while ((getline(&line_buff, &buff_size, fp)) != -1) {
        char *line_r; // Pointer for use with strtok_r function.
        // Parse the current line to extract the taxonomy ID.
        uint32_t curr_id = (uint32_t)strtol(strtok_r(line_buff, "|", &line_r), 0, 10);
        // Parse the next field to get the class, and convert it to lowercase.
        char *class = strlower(strtok_r(NULL, "|", &line_r));
        // Check if the current line's ID matches the specified ID.
        if(check_id == curr_id) {
            // Check if the current line is marked with "scientific name".
            if(strstr(line_r, "scientific name")) {
                // If the class matches the specified classification, return true.
                if(strstr(class, check_class)) return true;
                else break; // Exit the loop if there's a class mismatch.
            }
            else continue; // Continue to the next line if it's not a "scientific name" entry.
        }
    }

    // Clean up: close the file and free the allocated buffer.
    fclose(fp);
    free(line_buff);

    return false; // Return false if no match was found.
}

// This function converts all characters in a string to lowercase.
// It iterates over each character in the string, applying the tolower function.
char *strlower(char *str) {
    uint32_t len = strlen(str); // Get the length of the string.
    for(uint32_t i = 0; i < len; i++){
        str[i] = tolower(str[i]); // Convert each character to lowercase.
    }
    return str; // Return the modified string.
}