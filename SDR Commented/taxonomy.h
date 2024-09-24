// Prevents multiple inclusions of this header file
#ifndef taxonomy_h
#define taxonomy_h

// Includes the header file for sequencing data processing
#include "sequencing_data.h"

// Defines a type alias for a taxonomy container structure
typedef struct taxonomy_container TAX;

// Function declarations for managing TAX objects:
extern TAX *newTAX(uint32_t); // Creates a new TAX object
extern bool checkTAX(TAX *, uint32_t); // Checks for a specific condition in a TAX object
extern void insertTAX(TAX *, uint32_t); // Inserts a value into a TAX object
extern void freeTAX(TAX *); // Frees the memory allocated for a TAX object
extern void displayTAX(TAX *); // Displays the contents of a TAX object
extern uint32_t numHits(TAX *); // Returns the number of hits in a TAX object
extern uint32_t getTaxID(TAX *, uint32_t); // Retrieves a taxonomic ID from a TAX object

// Parses BLAST output to create TAX objects
extern TAX **parse_blast(char *, char **, uint32_t);

// Functions for checking the integrity and information of taxonomy data:
extern bool check_delnodes(const char *, const uint32_t); // Checks for deleted nodes in taxonomy data
extern void check_merged(const char *, uint32_t *); // Checks for merged nodes in taxonomy data
extern void check_nodes(const char *, uint32_t *, const uint32_t, const char *); // Checks nodes in taxonomy data
extern bool check_names(const char *, const uint32_t, const char *); // Checks names in taxonomy data

// Utility function to convert a string to lowercase
extern char *strlower(char *);

#endif // End of include guard
