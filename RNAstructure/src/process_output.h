#ifndef _HMM_OUTPUT_
#define _HMM_OUTPUT_

#define MAX_N_NUCS (5000)

// Manages the observations (outputs) for calculating probabilities, training, etc.
extern int length1;
extern int length2;

// Reads the next symbol from output file (or whatever) and returns the linear value of symbol.
// Accessing and indexing can be in any way so parameters in the get_symbol can be changed from one hmm to another.
void load_output(char* seq1_name, char* seq2_name);
int get_symbol(int nuc1, int nuc2);

// This function can be used for mapping particular symbol into linear symbol that
// is used while defining hmm model. Implementation of this function is completely output 
// format specific, for alignment model, it should implement a sequence interface and 
// read nucleotides from sequences.
int map_symbol_to_linear();

double get_trans_emit_prob(int first_state, int second_state, int n1, int n2);

void load_raw_output(char* _seq1_nucs, char* _seq2_nucs);

char get_nuc(int seq_index, int nuc_index);

char generate_random_nuc();

#endif
