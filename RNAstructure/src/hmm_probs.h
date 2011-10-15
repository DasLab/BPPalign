#ifndef _HMM_PROBS_
#define _HMM_PROBS_

#define N_STATES (3)
//#define N_OPS (3)

// Transition and emission probabilities are to be defined in hmm
// description file, it sohuld be noted that these two can 
// define a hmm by themselves, I do not need anything else.
// However a better description is needed to generalize for-back algortihm 
// because states and output might not be synchronized as in the case of 
// pairwise hmm for alignment, in that case gaps are not seen at the output.
// This makes automation harder because I need a mapping from particular 
// outputs (that is, the outputs that are seen) into the linear outputs that 
// define hmm in the definition file, the emission probabilities are given for
// linear output space, that is, the outputs are mapped onto the integers.

// Following function loads an hmm model from hmm definition file.
// Currently I am not using a hmm definition file, for future implementations.
bool load_hmm_model();

double get_trans_prob(int from, int to);
double get_log_trans_prob(int from, int to);

// get_emit_prob gets prob. of emission of a symbol 
// out of state in the parameter. symbol parameter must be 
// mapped into linear symbol in hmm definition file.
double get_emit_prob(int state, int symbol);
double get_log_emit_prob(int state, int symbol);

#endif

