#ifndef _HMM_ARRAYS_
#define _HMM_ARRAYS_

struct t_annot_aln_info;

enum{ALN_NACK, ALN_ACK, ALN_ACTUAL};

/*
Software topology is like this:
hmm_array.cpp calculates all arrays needed, DOES NOT DO ANY PROCESSING, it just makes arrays ready.
For processing these (for example using annotated alignments, like in RFAM) other programs can link to these
obj files (or a library for pairwise alignment.) Note that there is global variable global_aln_info which can be 
included by external linkage to other objs.

In that sense hmm_arrays.cpp should not be calculating any external data since all this can be done externally using 
global_aln_info. Instead hmm_arrays should be clearly well running concentrating on correctness and efficiency of 
pairwise alignment hmm. However ther is still the question of what services hmm_array.cpp should be offering, 
for instance is it necessary to include alignment envelope calculating code into hmm_array? NO. This should also be
an external service, the core hmm calculations (i.e. hmm_arrays.cpp) should be kept simple.
*/

/*
Following are core pairwise alignment hmm services.
*/

struct t_hmm_array
{
	int n_length1;
	int n_length2;

	double*** probs;

	double alignment_prob;
};

struct t_aln_info
{
	double op_prob;
	t_hmm_array* fore_hmm_array;
	t_hmm_array* back_hmm_array;

	double** aln_probs;
	double* ins1_probs;
	double* ins2_probs;
	double* ind_aln1_probs;
	double* ind_aln2_probs;

	// Alignment envelope parameters.
	//bool** aln_env;
	char** corr_aln_env; // This alignment envelope is the one that services should be using, service should not be changing it.
	char** dump_aln_env; // This is the alignment envelope to illustrate things on so that it is changeable without any problem
	int* smallest_M; // This is the smallest M parameter that yields M constraint which will cover all probabilistic aln. env. constraint.
	int n_aln_env;

	int length1;
	int length2;
};

// All following are core pairwise hmm operations which are not changeable from hmm_arrays.cpp
void init_global_aln_info();

void allocate_array(t_hmm_array* hmm_array);
void free_array(t_hmm_array* hmm_array);

void free_global_aln_info();

bool init_forward_array(t_hmm_array*, short**);
bool init_backward_array(t_hmm_array*);
bool init_backward_array2(t_hmm_array*, short**);

void calculate_forward_probs(t_hmm_array* fore_hmm_array, short**);
void calculate_backward_probs(t_hmm_array* back_hmm_array);
void calculate_backward_probs2(t_hmm_array* back_hmm_array, short**);

void calculate_aln_prob_array(t_hmm_array* fore_hmm_array, t_hmm_array* back_hmm_array);

void calculate_ins_prob_array(t_hmm_array* fore_hmm_array, t_hmm_array* back_hmm_array);
void calculate_ins_prob_array2(t_hmm_array* fore_hmm_array, t_hmm_array* back_hmm_array);

void verify_forward_backward();

void get_aln_permissions(short** forcealign, 
			bool& forbid_STATE_ALN, 
			bool& forbid_STATE_INS1, 
			bool& forbid_STATE_INS2, 
			int i, 
			int k);
//////////////////////////////////////////////////////////////////////////////////////////////////

#endif
