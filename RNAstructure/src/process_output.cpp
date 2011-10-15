#include "process_output.h"
#include <time.h>
#include "structure.h"
#include <string.h>
#include "alignment_hmm_model.h"
#include "hmm_probs.h"
#include "hmm_arrays.h"
#include "xlog_math.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>

using namespace std;

/*
static char seq1_nucs[MAX_N_NUCS];
static char seq2_nucs[MAX_N_NUCS];
*/

static char* seq1_nucs;
static char* seq2_nucs;

int length1;
int length2;

extern t_aln_info global_aln_info;

extern double emit_probs[N_OUTPUTS][N_STATES];
extern double trans_probs[N_STATES][N_STATES];

// NOTE THAT OUTPUT SEQUENCES ARE INDEXED STARTING FROM 0 HOWEVER HMM OUTPUTS ARE INDEXED STARTING FROM 1 IN THE CODE.
// Assuming that outputs are in .seq format.
void load_output(char* seq1_name, char* seq2_name)
{
	//t_sequence* seq1 = new t_sequence(NULL, seq1_name);
	//t_sequence* seq2 = new t_sequence(NULL, seq2_name);
	structure seq1;
	structure seq2;

	openseq(&seq1, seq1_name);
	openseq(&seq2, seq2_name);

	seq1_nucs = (char*)malloc(sizeof(char) * (seq1.numofbases + 3));
	seq2_nucs = (char*)malloc(sizeof(char) * (seq2.numofbases + 3));

	printf("Loading output with %d  and %d nucleotides...\n", seq1.numofbases, seq2.numofbases);

	srand((unsigned int) time(0));
	
	// Take care of unknown nucleotides by randomly selecting them.
	//for(int cnt = 0; cnt < strlen(seq1.nucs); cnt++)
	for(int cnt = 1; cnt < seq1.numofbases; cnt++)
	{
		seq1.nucs[cnt] = toupper(seq1.nucs[cnt]);		
		if(seq1.nucs[cnt] != 'A' &&
		seq1.nucs[cnt] != 'C' &&
		seq1.nucs[cnt] != 'G' &&
                seq1.nucs[cnt] != 'T'&&
		seq1.nucs[cnt] != 'U')
		{
			seq1.nucs[cnt] = generate_random_nuc();
		}
	}

	//for(int cnt = 0; cnt < strlen(seq2.nucs); cnt++)
	for(int cnt = 1; cnt < seq2.numofbases; cnt++)
	{
		seq2.nucs[cnt] = toupper(seq2.nucs[cnt]);
		if(seq2.nucs[cnt] != 'A' &&
		seq2.nucs[cnt] != 'C' &&
		seq2.nucs[cnt] != 'G' &&
                seq2.nucs[cnt] != 'T' &&
		seq2.nucs[cnt] != 'U')
		{
			seq2.nucs[cnt] = generate_random_nuc();
		}		
	}

	strcpy(seq1_nucs, &seq1.nucs[1]);
	strcpy(seq2_nucs, &seq2.nucs[1]);

	printf("individual nucs:\n%s\n%s\n", seq1_nucs, seq2_nucs);

	FILE* prob_log_file = fopen("prob_log.txt", "a");
	//fprintf(prob_log_file, "\nAligning\n%s\n%s\n\n", seq1_nucs, seq2_nucs);
	fclose(prob_log_file);

	length1 = seq1.numofbases;
	length2 = seq2.numofbases;

	global_aln_info.length1 = length1;
	global_aln_info.length2 = length2;
}

// Assuming that outputs are read externally and loaded to here.
void load_raw_output(char* _seq1_nucs, char* _seq2_nucs)
{
/*	
        printf("Loading output with %d  and %d nucleotides...\n", strlen(_seq1_nucs+1), strlen(_seq2_nucs+1));
	printf("seq1_nucs: %s\n\n", _seq1_nucs + 1);
	printf("seq2 nucs: %s\n", _seq2_nucs + 1);
*/
	
	seq1_nucs = (char*)malloc(sizeof(char) * (strlen(_seq1_nucs + 1) + 3));
	seq2_nucs = (char*)malloc(sizeof(char) * (strlen(_seq2_nucs + 1) + 3));

        // Following shift in output is for taking care of changing indices between hmm caclulations and structure class.
        strcpy(seq1_nucs, _seq1_nucs + 1);
        strcpy(seq2_nucs, _seq2_nucs + 1);

        length1 = (int) strlen(seq1_nucs);
        length2 = (int) strlen(seq2_nucs);

        // Take care of unknown nucleotides by randomly selecting them.
        //for(int cnt = 0; cnt < strlen(seq1.nucs); cnt++)
        for(int cnt = 0; cnt < length1; cnt++)
        {
                seq1_nucs[cnt] = toupper(seq1_nucs[cnt]);
                if(seq1_nucs[cnt] != 'A' &&
                seq1_nucs[cnt] != 'C' &&
                seq1_nucs[cnt] != 'G' &&
                seq1_nucs[cnt] != 'T' &&
                seq1_nucs[cnt] != 'U')
                {
                        seq1_nucs[cnt] = generate_random_nuc();
                }
//                printf("%c", seq1_nucs[cnt]);
        }

//        printf("\n");

        //for(int cnt = 0; cnt < strlen(seq2.nucs); cnt++)
        for(int cnt = 0; cnt < length2; cnt++)
        {
                seq2_nucs[cnt] = toupper(seq2_nucs[cnt]);
                if(seq2_nucs[cnt] != 'A' &&
                seq2_nucs[cnt] != 'C' &&
                seq2_nucs[cnt] != 'G' &&
                seq2_nucs[cnt] != 'T' &&
                seq2_nucs[cnt] != 'U')
                {
                        seq2_nucs[cnt] = generate_random_nuc();
                }
//                printf("%c", seq2_nucs[cnt]);
        }

//        printf("\n");

/*
	strcpy(seq1_nucs, _seq1_nucs);
	strcpy(seq2_nucs, _seq2_nucs);

	length1 = strlen(seq1_nucs);
	length2 = strlen(seq2_nucs);
*/
	global_aln_info.length1 = length1;
	global_aln_info.length2 = length2;

	//FILE* prob_log_file = fopen("prob_log.txt", "a");
	//fprintf(prob_log_file, "\nAligning\n%s\n%s\n\n", seq1_nucs, seq2_nucs);
	//printf("\nAligning\n%s\n%s\n\n", seq1_nucs, seq2_nucs);
	//cin.get();
	//fclose(prob_log_file);
}

// Return linear symbol value.
// I am indexing nucs from 1-N1 in probability calculations, in sequence interface they are indexed
// from 0-(N1-1), to compensate for that, subtract 1 from nuc1 and nuc2 if they are not -1 (meaning a gap).
int get_symbol(int nuc1, int nuc2)
{
	//printf("nuc1 = %d, nuc2 = %d, length1 = %d, length2 = %d\n", nuc1, nuc2, length1, length2);

	if(nuc1 == length1 + 1 && nuc2 == length2 + 1)
	{
		return(25); // 26th symbol is ending symbol, emitted only by STATE_ALN state.
	}

	if(nuc1 == 0 && nuc2 == 0)
	{
		return(26); // 27th symbol is starting symbol, emitted only by STATE_ALN state.
	}

	// If this is not ending or starting symbol, treat the ends as gaps.
	if(nuc1 == length1 + 1 || nuc1 == 0)
	{
		nuc1 = -1;
	}
	
	if(nuc2 == length2 + 1 || nuc2 == 0)
	{
		nuc2 = -1;
	}

	if(nuc1 != -1 && nuc2 != -1)
	{
		nuc1--;
		nuc2--;
		if(nuc1 < 0 || nuc2 < 0)
		{
			exit(0);
		}
		//printf("%c\n%c", seq1_nucs[nuc1],seq2_nucs[nuc2] );
	}	
	else if(nuc1 == -1 && nuc2 != -1)
	{
		nuc2--;
		//printf("-\n%c", seq2_nucs[nuc2] );
	}
	else if(nuc2 == -1 && nuc1 != -1)
	{
		nuc1--;
		//printf("%c\n-", seq1_nucs[nuc1] );
	}
	else
	{
		//printf("-\n-", seq1_nucs[nuc1] );
	}
	
	int digit1, digit2;
	
	if(nuc1 == -1)
	{
		digit1 = 4;
	}
	else	
	{
		switch(seq1_nucs[nuc1])
		{
			case 'A':
				digit1 = 0;
				break;
			case 'C':
				digit1 = 1;
				break;
			case 'G':
				digit1 = 2;
				break;
			case 'U':
				digit1 = 3;
				break;
                        case 'T':
                                digit1 = 3;
                                break;
			default:
				printf("Problematic nucleotide: %c at %s(%d)\n", seq1_nucs[nuc1], __FILE__, __LINE__);
				exit(0);
				break;
		}
	}

	if(nuc2 == -1)
	{
		digit2 = 4;
	}
	else	
	{
		switch(seq2_nucs[nuc2])
		{
			case 'A':
				digit2 = 0;
				break;
			case 'C':
				digit2 = 1;
				break;
			case 'G':
				digit2 = 2;
				break;
			case 'U':
				digit2 = 3;
				break;
                        case 'T':
                                digit2 = 3;
                                break;
			default:
				printf("Problematic nucleotide at %d: %c at %s(%d)\n", nuc2, seq2_nucs[nuc2], __FILE__, __LINE__);
				exit(0);
				break;
		}
	}

	int index = digit1 * 5 + digit2;

	//printf("\n%d\n", index);

	return(index);
}

// n1 and n2 are assumed to be the next OBSERVED outputs after first_state's outputs.
// transition from first_state to second_state and emission based on second_state of symbols starting fron n1 and n2.
double get_trans_emit_prob(int first_state, int second_state, int n1, int n2)
{	
	// If this is ending state, return 1/3 transition probability.
	// Since I am doing an emission in ending state while calculating backward variable,
	// I have to assign correct trans_emit probabilities for that because otherwise
	// emission of ending symbols is calculated as 0 probability, which is wrong. 
	/*
	if(n1 == length1+1 && n2 == length2+1 && second_state == STATE_ALN)
	{
		//return(xlog((double)1/N_STATES));
		return(xlog((double)0));
	}
	*/
	
	// For starting state, I am assuming the same transition probabilities are applied, however
	// starting state does not emit any symbols so I do not need to handle the case of 0,0 (aligned state)
	// however I am assuming that the probability of starting actual emitting states is the transition 
	// probability into these states from an aligned state, I will not take this probability as 
	// uniform.

	// Must check whether this is the ending state for backward variable and starting state for forward variable.
	// transition probability.
	double trans_prob = get_log_trans_prob(first_state, second_state);
	//printf("transition_prob = %f, emitting nucs @ indices (%d, %d)\n", trans_prob, n1, n2);

	//printf("n1 = %d, n2 = %d\n", n1, n2);

	// emission of second_state based on n1 and n2 and emitted symbol depends on second state.
	int linear_sym;

	if(n1 == 0 && n2 == 0)
	{
		printf("Exiting, 0 indices problem..\n");
		exit(0);
	}

	switch(second_state)
	{
		// Emit n1 and a gap.
		case(STATE_INS1):
			linear_sym = get_symbol(n1, -1);	
		break;

		case(STATE_INS2):
			linear_sym = get_symbol(-1, n2);
		break;

		case(STATE_ALN):
			linear_sym = get_symbol(n1, n2);
		break;
	}

	double emit_prob = get_log_emit_prob(second_state, linear_sym);

	return(xlog_mul(trans_prob, emit_prob));
}


char get_nuc(int seq_index, int nuc_index)
{
	if(seq_index == 1)
	{
		return(seq1_nucs[nuc_index]);
	}
	else
	{
		return(seq2_nucs[nuc_index]);
	}
}

// Following function generates random nucleotide in place of an unknown nucleotide.
char generate_random_nuc()
{
	int rand_nuc = rand() % 4;

	switch(rand_nuc)
	{
	case 0:
		return('A');
		break;
	case 1:
		return('C');
		break;
	case 2:
		return('G');
		break;
	case 3:
		return('U');
		break;
	default:
		printf("Invalid random nuc @ %s(%d).\n", __FILE__, __LINE__);
		exit(0);
		break;
	};
}
