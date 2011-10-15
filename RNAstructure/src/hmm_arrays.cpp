#include <stdio.h>
#include <stdlib.h>
#include "hmm_arrays.h"
#include "hmm_probs.h"
#include "alignment_hmm_model.h"
#include "xlog_math.h"
#include "process_output.h"
#include <iostream>
#include <math.h>

using namespace std;

t_aln_info global_aln_info;

extern double emit_probs[N_OUTPUTS][N_STATES];
extern double trans_probs[N_STATES][N_STATES];

extern int length1;
extern int length2;

void init_global_aln_info()
{
	global_aln_info.aln_probs = NULL;
	global_aln_info.ins1_probs = NULL;
	global_aln_info.ins2_probs = NULL;
	global_aln_info.ind_aln1_probs = NULL;
	global_aln_info.ind_aln2_probs = NULL;
	global_aln_info.corr_aln_env = NULL;
	global_aln_info.dump_aln_env = NULL;
	global_aln_info.smallest_M = NULL;

	global_aln_info.fore_hmm_array = NULL;
	global_aln_info.back_hmm_array = NULL;
}

void allocate_array(t_hmm_array* hmm_array)
{
	int n_nucs1 = length1;
	int n_nucs2 = length2;
	hmm_array->n_length1 = n_nucs1;
	hmm_array->n_length2 = n_nucs2;

	double*** _hmm_array = (double***)malloc( sizeof(double**) * (n_nucs1+2) );

	if(_hmm_array == NULL)
	{
		printf("malloc failed at %s (%d)\n", __FILE__, __LINE__);
	}

	for(int cnt = 0; cnt < n_nucs1+2; cnt++)
	{
		_hmm_array[cnt] = (double**)malloc( sizeof(double*) * (n_nucs2+2) );

		if(_hmm_array[cnt] == NULL)
		{
			printf("malloc failed at %s (%d)\n", __FILE__, __LINE__);
		}

		// Now allocate N_STATES many states for each 
		for(int cnt2 = 0; cnt2 < n_nucs2+2; cnt2++)
		{
			_hmm_array[cnt][cnt2] = (double*)malloc( sizeof(double) * N_STATES );
		
			for(int i_init = 0; i_init < N_STATES; i_init++)
			{
				_hmm_array[cnt][cnt2][i_init] = xlog(0.0f);
			}
		}
	}

	hmm_array->probs = _hmm_array;
}

void free_array(t_hmm_array* hmm_array)
{
	//double*** _hmm_array = (double***)malloc( sizeof(double**) * (n_nucs1+2) );

	int n_nucs1 = hmm_array->n_length1;
	int n_nucs2 = hmm_array->n_length2;

	if(hmm_array == NULL)
	{
		printf("malloc failed at %s (%d)\n", __FILE__, __LINE__);
	}

	for(int cnt = 0; cnt < n_nucs1+2; cnt++)
	{
		//_hmm_array[cnt] = (double**)malloc( sizeof(double*) * (n_nucs2+2) );

		// Now allocate N_STATES many states for each 
		for(int cnt2 = 0; cnt2 < n_nucs2+2; cnt2++)
		{
			//_hmm_array[cnt][cnt2] = (double*)malloc( sizeof(double) * N_STATES );
			free(hmm_array->probs[cnt][cnt2]);
		}

		free(hmm_array->probs[cnt]);
	}

	free(hmm_array->probs);	

	//hmm_array->n_length1 = -1;
	//hmm_array->n_length2 = -1;
}

// Parameter to following function should be a forward hmm array.
bool init_forward_array(t_hmm_array* hmm_array, short** forcealign)
{
	global_aln_info.fore_hmm_array = hmm_array;

	// Set 0,0 alignment as 1 prob.
	hmm_array->probs[0][0][STATE_ALN] = xlog(1.0);
	hmm_array->probs[0][0][STATE_INS1] = xlog(0);
	hmm_array->probs[0][0][STATE_INS2] = xlog(0);

	//printf("init. (%d, %d) = %f\n", 0, 0, hmm_array->probs[0][0][STATE_ALN] );

	for(int cnt = 1; cnt <= hmm_array->n_length1; cnt++)
	{
		// Get prob. state for emission, assuming that this is 1 for now.
		//double emit_prob = xlog(1.0/24);		

		if(cnt > 1)
		{
			// Transition probability between forward boundary state.
			//double trans_prob = get_log_trans_prob(STATE_INS1, STATE_INS1);

			double trans_emit_prob = get_trans_emit_prob(STATE_INS1, STATE_INS1, cnt, 0);

			hmm_array->probs[cnt][0][STATE_ALN] = xlog(0);
			//hmm_array->probs[cnt][0][STATE_INS1] = xlog_mul(hmm_array->probs[cnt - 1][0][STATE_INS1], xlog_mul(trans_prob, emit_prob));
			hmm_array->probs[cnt][0][STATE_INS1] = xlog_mul(hmm_array->probs[cnt - 1][0][STATE_INS1], trans_emit_prob);
			hmm_array->probs[cnt][0][STATE_INS2] = xlog(0);
		}
		else // if cnt is 1 initing special case, we are doing a switch from alignment state (state of 0,0) into ins1 state.
		{
			//double trans_prob = get_log_trans_prob(STATE_ALN, STATE_INS1);
			double trans_emit_prob = get_trans_emit_prob(STATE_ALN, STATE_INS1, cnt, 0);

			hmm_array->probs[cnt][0][STATE_ALN] = xlog(0);
			//hmm_array->probs[cnt][0][STATE_INS1] = xlog_mul(0 , xlog_mul(trans_prob, emit_prob));
			hmm_array->probs[cnt][0][STATE_INS1] = xlog_mul(0 , trans_emit_prob);
			hmm_array->probs[cnt][0][STATE_INS2] = xlog(0);
		}

		//printf("init. (%d, %d) = %f\n", cnt, 0, hmm_array->probs[cnt][0][STATE_INS1] );
	}

	for(int cnt = 1; cnt <= hmm_array->n_length2; cnt++)
	//for(cnt = 1; cnt <= hmm_array->n_length2; cnt++)
	{
		// Get prob. state for emission, assuming that this is 1 for now.
		//double emit_prob = xlog(1.0/24);

		if(cnt > 1)
		{
			// Transition probability between forward boundary state.
			//double trans_prob = get_log_trans_prob(STATE_INS2, STATE_INS2);
			double trans_emit_prob = get_trans_emit_prob(STATE_INS2, STATE_INS2, 0, cnt);

			hmm_array->probs[0][cnt][STATE_ALN] = xlog(0);
			hmm_array->probs[0][cnt][STATE_INS1] = xlog(0);
			//hmm_array->probs[0][cnt][STATE_INS2] = xlog_mul(hmm_array->probs[0][cnt - 1][STATE_INS2], xlog_mul(trans_prob, emit_prob));
			hmm_array->probs[0][cnt][STATE_INS2] = xlog_mul(hmm_array->probs[0][cnt - 1][STATE_INS2], trans_emit_prob);
		}
		else // if cnt is 1 initing special case, we are doing a switch from alignment state (state of 0,0) into ins1 state.
		{
			//double trans_prob = get_log_trans_prob(STATE_ALN, STATE_INS2);
			double trans_emit_prob = get_trans_emit_prob(STATE_ALN, STATE_INS2, 0, cnt);

			hmm_array->probs[0][cnt][STATE_ALN] = xlog(0);
			hmm_array->probs[0][cnt][STATE_INS1] = xlog(0);
			//hmm_array->probs[0][cnt][STATE_INS2] = xlog_mul(0, xlog_mul(trans_prob, emit_prob));
			hmm_array->probs[0][cnt][STATE_INS2] = xlog_mul(0, trans_emit_prob);
		}	

		//printf("init. (%d, %d) = %f\n", 0, cnt, hmm_array->probs[0][cnt][STATE_INS2] );
	}

	return(true);
}

void calculate_forward_probs(t_hmm_array* hmm_array, short** forcealign)
{
	//printf("Calculating forward variable...\n");

	// Order of calculation is important. Do not recalculate boundaries, use them.
	for(int cnt1 = 1; cnt1 <= hmm_array->n_length1; cnt1++)
	{
		for(int cnt2 = 1; cnt2 <= hmm_array->n_length2; cnt2++)
		{	
			bool forbid_STATE_INS1 = false;
			bool forbid_STATE_INS2 = false;
			bool forbid_STATE_ALN = false;

			get_aln_permissions(forcealign,
				forbid_STATE_ALN, 
				forbid_STATE_INS1, 
				forbid_STATE_INS2, 
				cnt1,
				cnt2);
			
			// This loop is for iterating over possible states in this alignment pair.
			for(int current_state = 0; current_state < N_STATES; current_state++)
			{
				//printf("------------------------------------\n");

				hmm_array->probs[cnt1][cnt2][current_state] = xlog(0);

				// This loop is for marginalizing over previous state.
				for(int prev_state = 0; prev_state < N_STATES; prev_state++)
				{
					//double trans_prob = get_log_trans_prob(prev_state, current_state);

					// This state is emitting starting with cnt1 and cnt2.
					//get_trans_emit_prob(prev_state, current_state, cnt1, cnt2);

					if(!forbid_STATE_ALN &&
						current_state == STATE_ALN)
					{
						//double emit_prob = xlog(1.0/24); // Take emission prob. as 1 for now.						

						//printf("prev_state = %s\n", state_names[prev_state]);
						//printf("%s -> %s\n", state_names[prev_state], state_names[current_state]);
						double trans_emit_prob = get_trans_emit_prob(prev_state, current_state, cnt1, cnt2);
						//printf("Recursing on (%d, %d) with prob. %f.\n", cnt1-1, cnt2-1, hmm_array->probs[cnt1-1][cnt2-1][prev_state]);

						hmm_array->probs[cnt1][cnt2][current_state] = xlog_sum(hmm_array->probs[cnt1][cnt2][current_state], xlog_mul(hmm_array->probs[cnt1-1][cnt2-1][prev_state], trans_emit_prob));
					}
					else if(!forbid_STATE_INS1 && 
						current_state == STATE_INS1)
					{
						//double emit_prob = xlog(1.0/24); // Take emission prob. as 1 for now.
						//get_symbol(cnt1, -1);

						//printf("prev_state = %s\n", state_names[prev_state]);
						//printf("%s -> %s\n", state_names[prev_state], state_names[current_state]);
						double trans_emit_prob = get_trans_emit_prob(prev_state, current_state, cnt1, cnt2);
						//printf("Recursing on (%d, %d) with prob. %f.\n", cnt1-1, cnt2, hmm_array->probs[cnt1-1][cnt2][prev_state]);

						//hmm_array->probs[cnt1][cnt2][current_state] = xlog_sum(hmm_array->probs[cnt1][cnt2][current_state], xlog_mul(hmm_array->probs[cnt1-1][cnt2][prev_state], xlog_mul(trans_prob, emit_prob)));
						hmm_array->probs[cnt1][cnt2][current_state] = xlog_sum(hmm_array->probs[cnt1][cnt2][current_state], xlog_mul(hmm_array->probs[cnt1-1][cnt2][prev_state], trans_emit_prob));
					}
					else if(!forbid_STATE_INS2)
					{
						//double emit_prob = xlog(1.0/24); // Take emission prob. as 1 for now.
						//get_symbol(-1, cnt2);

						//printf("prev_state = %s\n", state_names[prev_state]);
						//printf("%s -> %s\n", state_names[prev_state], state_names[current_state]);
						double trans_emit_prob = get_trans_emit_prob(prev_state, current_state, cnt1, cnt2);
						//printf("Recursing on (%d, %d) with prob. %f.\n", cnt1, cnt2-1, hmm_array->probs[cnt1][cnt2-1][prev_state]);

						hmm_array->probs[cnt1][cnt2][current_state] = xlog_sum(hmm_array->probs[cnt1][cnt2][current_state], xlog_mul(hmm_array->probs[cnt1][cnt2-1][prev_state], trans_emit_prob));
					}
				}
				
				//printf("%d, %d, current state: %s, prob = %f (%f)\n", cnt1, cnt2, state_names[current_state], hmm_array->probs[cnt1][cnt2][current_state], xexp(hmm_array->probs[cnt1][cnt2][current_state]) );
			}
		}
	}

	// Calculate output probabilities.
	double output_prob = LOG_OF_ZERO;

	for(int cnt = 0; cnt < N_STATES; cnt++)
	{
		double last_state_prob = get_log_trans_prob(cnt, STATE_ALN); // Prob. of finishing alignment with this state.
		output_prob = xlog_sum(output_prob, xlog_mul(hmm_array->probs[hmm_array->n_length1][hmm_array->n_length2][cnt], last_state_prob));
	}

	//printf("Output probability of alignment from n+1sts is %f\n", hmm_array->probs[hmm_array->n_length1 + 1][hmm_array->n_length2 + 1][STATE_ALN]);
	//printf("Output probability of alignment from n+1sts is %f\n", hmm_array->probs[hmm_array->n_length1 + 1][hmm_array->n_length2 + 1][STATE_INS1]);
	//printf("Output probability of alignment from n+1sts is %f\n", hmm_array->probs[hmm_array->n_length1 + 1][hmm_array->n_length2 + 1][STATE_INS2]);

	global_aln_info.op_prob = output_prob; // Save output prob.

	//cin.get();
}

// Allocate and calculate probability array for alignment pairs.
void calculate_aln_prob_array(t_hmm_array* fore_hmm_array, t_hmm_array* back_hmm_array)
{
	double** aln_pair_probs = (double**)malloc((fore_hmm_array->n_length1+1) * sizeof(double*));

	global_aln_info.aln_probs = aln_pair_probs;

	//FILE* prob_log_file = fopen("prob_log.txt", "a");
	//FILE* posterior_prob_file = fopen("prob_matrix.txt", "w");

	//fprintf(prob_log_file, "Posterior probs:\n");
	//fprintf(posterior_prob_file, "Posterior probs:\n");

	//for(int cnt1 = 1; cnt1 < fore_hmm_array->n_length1+1; cnt1++)
	for(int cnt1 = 0; cnt1 < fore_hmm_array->n_length1+1; cnt1++)
	{
		*(aln_pair_probs + cnt1) = (double*)malloc((fore_hmm_array->n_length2+1) * sizeof(double));

		//for(int cnt2 = 1; cnt2 < fore_hmm_array->n_length2+1; cnt2++)
		for(int cnt2 = 0; cnt2 < fore_hmm_array->n_length2+1; cnt2++)
		{
			// Calculate alignment pair probabilities by multiplying 
			// corresponding forward and backward variable array entries.
			aln_pair_probs[cnt1][cnt2] = xlog_mul(fore_hmm_array->probs[cnt1][cnt2][STATE_ALN], back_hmm_array->probs[cnt1][cnt2][STATE_ALN]);

			//printf("prob of aligning %d, %d is %f\n", cnt1, cnt2, aln_pair_probs[cnt1][cnt2]);

			//fprintf(prob_log_file, "prob of aligning %d, %d is %f\n", cnt1, cnt2, aln_pair_probs[cnt1][cnt2]);

			//fprintf(posterior_prob_file, "%-6.3f ", aln_pair_probs[cnt1][cnt2]);
		}

		//fprintf(posterior_prob_file, "\n");
	}

	//fclose(prob_log_file);

	//fclose(posterior_prob_file);
}

void verify_forward_backward()
{
	printf("\nCalculating alignment and insertion probabilities of 1st sequence..\n");

	for(int cnt1 = 1; cnt1 < global_aln_info.fore_hmm_array->n_length1+1; cnt1++)
	{
		double ins_prob = global_aln_info.ins1_probs[cnt1];

		double aln_prob = LOG_OF_ZERO;

		for(int cnt2 = 1; cnt2 < global_aln_info.fore_hmm_array->n_length2+1; cnt2++)
		{		
			aln_prob = xlog_sum(aln_prob, global_aln_info.aln_probs[cnt1][cnt2]);
		}

		printf("%d. base in seq1, prob. of aln + insertion = %f + %f = %.20f, prob. of o/p = %.20f\n", cnt1, aln_prob, ins_prob, xlog_sum(aln_prob, ins_prob), global_aln_info.op_prob);

		if(xlog_sum(aln_prob, ins_prob) != global_aln_info.op_prob)
		{
			printf("\tWrong prob. calculation at %d. base in 1st sequence. \n", cnt1);
		}
	}

	printf("\nCalculating alignment and insertion probabilities of 2nd sequence..\n");

	for(int cnt2 = 1; cnt2 < global_aln_info.fore_hmm_array->n_length2+1; cnt2++)
	{
		double ins_prob = global_aln_info.ins2_probs[cnt2];

		double aln_prob = LOG_OF_ZERO;

		for(int cnt1 = 1; cnt1 < global_aln_info.fore_hmm_array->n_length1+1; cnt1++)
		{		
			aln_prob = xlog_sum(aln_prob, global_aln_info.aln_probs[cnt1][cnt2]);
		}

		printf("%d. base in seq2, prob. of aln + insertion = %f + %f = %.20f, prob. of o/p = %.20f\n", cnt2, aln_prob, ins_prob, xlog_sum(aln_prob, ins_prob), global_aln_info.op_prob);

		if(xlog_sum(aln_prob, ins_prob) != global_aln_info.op_prob)
		{
			printf("\tWrong prob. calculation at %d. base in 2nd sequence. \n", cnt2);
		}	
	}
}

void free_global_aln_info()
{

	if(global_aln_info.aln_probs != NULL)
	{
		for(int cnt1 = 0; cnt1 < global_aln_info.fore_hmm_array->n_length1+1; cnt1++)
		{
			free( global_aln_info.aln_probs[cnt1] );
		}

		free(global_aln_info.aln_probs);
	}


	if(global_aln_info.corr_aln_env != NULL)
	{
		for(int cnt = 1; cnt < global_aln_info.fore_hmm_array->n_length1+1; cnt++)
		{
			free( global_aln_info.corr_aln_env[cnt] );
			free( global_aln_info.dump_aln_env[cnt] );
		}	

		free(global_aln_info.corr_aln_env);	
		free(global_aln_info.dump_aln_env);	
	}

	if(global_aln_info.smallest_M != NULL)
	{
		free( global_aln_info.smallest_M );
	}

/*
	global_aln_info.aln_probs = NULL;
	global_aln_info.ins1_probs = NULL;
	global_aln_info.ins2_probs = NULL;
	global_aln_info.ind_aln1_probs = NULL;
	global_aln_info.ind_aln2_probs = NULL;
	global_aln_info.corr_aln_env = NULL;
	global_aln_info.smallest_M = NULL;

	global_aln_info.fore_hmm_array = NULL;
	global_aln_info.back_hmm_array = NULL;
*/
}

void get_aln_permissions(short** forcealign,
                        bool& forbid_STATE_ALN,
                        bool& forbid_STATE_INS1,
                        bool& forbid_STATE_INS2,
                        int i,
                        int k)
{
	if(forcealign == NULL)
	{
		forbid_STATE_ALN = false;
		forbid_STATE_INS1 = false;
		forbid_STATE_INS2 = false;			
	}
	else
	{
		short* seq1_aln_const = forcealign[0]; 
		short* seq2_aln_const = forcealign[1]; 


		if(seq1_aln_const[i] != 0)
		{
			if(seq1_aln_const[i] == k)
			{
				forbid_STATE_ALN = false;
				forbid_STATE_INS1 = true;
				forbid_STATE_INS2 = true;
			}
			else
			{
				if(seq1_aln_const[i] < k)
				{
					forbid_STATE_ALN = true;
					forbid_STATE_INS1 = true;
					forbid_STATE_INS2 = false;
				}
				else
				{
					forbid_STATE_ALN = true;
					forbid_STATE_INS1 = true;
					forbid_STATE_INS2 = true;
				}
			}
		}
		else if(seq2_aln_const[k] != 0)
		{
                        if(seq2_aln_const[k] == i)
                        {
                                forbid_STATE_ALN = false;
                                forbid_STATE_INS1 = true;
                                forbid_STATE_INS2 = true;
                        }
                        else
                        {

				if(seq2_aln_const[k] < i)
				{
					forbid_STATE_ALN = true;
					forbid_STATE_INS1 = false;
					forbid_STATE_INS2 = true;
				}
				else
				{
					forbid_STATE_ALN = true;
					forbid_STATE_INS1 = true;
					forbid_STATE_INS2 = true;
				}
			}
		}
		else 
		{
			forbid_STATE_ALN = false;
			forbid_STATE_INS1 = false;
			forbid_STATE_INS2 = false;				
		}
	} 
}

// Backward variable calculations:
// Parameter to following function should be a backward hmm array.
bool init_backward_array2(t_hmm_array* hmm_array, short** forcealign)
{
	global_aln_info.back_hmm_array = hmm_array;

	// Set 0,0 alignment as 1 prob.
	int n1 = hmm_array->n_length1;
	int n2 = hmm_array->n_length2;
	hmm_array->probs[n1+1][n2+1][STATE_ALN] = xlog(1.0);
	hmm_array->probs[n1+1][n2+1][STATE_INS1] = xlog(0);
	hmm_array->probs[n1+1][n2+1][STATE_INS2] = xlog(0);

	for(int cnt = n1; cnt >= 1; cnt--)
	{
		// Get prob. state for emission, assuming that this is 1 for now.
		//double emit_prob = xlog(1.0/24);		

		if(cnt < n1)
		{
			// Transition probability between forward boundary state.
			//double trans_prob = get_log_trans_prob(STATE_INS1, STATE_INS1);
			// Emitted bases are cnt+1, n2+1 for ins1.
			double trans_emit_prob = get_trans_emit_prob(STATE_INS1, STATE_INS1, cnt+1, n2+1);

			hmm_array->probs[cnt][n2+1][STATE_ALN] = xlog(0);
			//hmm_array->probs[cnt][n2+1][STATE_INS1] = xlog_mul(hmm_array->probs[cnt+1][n2+1][STATE_INS1], xlog_mul(trans_prob, emit_prob));
			hmm_array->probs[cnt][n2+1][STATE_INS1] = xlog_mul(hmm_array->probs[cnt+1][n2+1][STATE_INS1], trans_emit_prob);
			hmm_array->probs[cnt][n2+1][STATE_INS2] = xlog(0);
		}
		else // if cnt is 1 initing special case, we are doing a switch from ins1 state (state of 0,0) into alignment state.
		{
			//double trans_prob = get_log_trans_prob(STATE_ALN, STATE_INS1);
			double trans_emit_prob = get_trans_emit_prob(STATE_INS1, STATE_ALN, cnt+1, n2+1);

			hmm_array->probs[cnt][n2+1][STATE_ALN] = xlog(0);
			//hmm_array->probs[cnt][n2+1][STATE_INS1] = xlog_mul(0, xlog_mul(trans_prob, emit_prob));
			hmm_array->probs[cnt][n2+1][STATE_INS1] = xlog_mul(0, trans_emit_prob);
			hmm_array->probs[cnt][n2+1][STATE_INS2] = xlog(0);
		}

		//printf("init. (%d, %d) = %f\n", cnt, n2+1, hmm_array->probs[cnt][n2+1][STATE_INS1] );
	}

	for(int cnt = n2; cnt >= 1; cnt--)
	//for(cnt = n2; cnt >= 1; cnt--)
	{
		// Get prob. state for emission, assuming that this is 1 for now.
		//double emit_prob = xlog(1.0/24);

		if(cnt < n2)
		{
			// Transition probability between forward boundary state.
			//double trans_prob = get_log_trans_prob(STATE_INS2, STATE_INS2);

			double trans_emit_prob = get_trans_emit_prob(STATE_INS2, STATE_INS2, n1+1, cnt+1);

			hmm_array->probs[n1+1][cnt][STATE_ALN] = xlog(0);
			hmm_array->probs[n1+1][cnt][STATE_INS1] = xlog(0);
			//hmm_array->probs[n1+1][cnt][STATE_INS2] = xlog_mul(hmm_array->probs[n1+1][cnt+1][STATE_INS2], xlog_mul(trans_prob, emit_prob));
			hmm_array->probs[n1+1][cnt][STATE_INS2] = xlog_mul(hmm_array->probs[n1+1][cnt+1][STATE_INS2], trans_emit_prob);
		}
		else // if cnt is 1 initing special case, we are doing a switch from alignment state (state of 0,0) into ins1 state.
		{
			//double trans_prob = get_log_trans_prob(STATE_ALN, STATE_INS2);
			double trans_emit_prob = get_trans_emit_prob(STATE_INS2, STATE_ALN, n1+1, cnt+1);

			hmm_array->probs[n1+1][cnt][STATE_ALN] = xlog(0);
			hmm_array->probs[n1+1][cnt][STATE_INS1] = xlog(0);
			//hmm_array->probs[n1+1][cnt][STATE_INS2] = xlog_mul(0, xlog_mul(trans_prob, emit_prob));
			hmm_array->probs[n1+1][cnt][STATE_INS2] = xlog_mul(0, trans_emit_prob);
		}

		//printf("init. (%d, %d) = %f\n", n1+1, cnt, hmm_array->probs[n1+1][cnt][STATE_INS2]);
	}

	return(true);
}

// Note that I am not considering emission of current state, considering emission of next state
// based on current set of observed outputs.
void calculate_backward_probs2(t_hmm_array* hmm_array, short** forcealign)
{
	for(int cnt1 = hmm_array->n_length1; cnt1 >= 0; cnt1--)
	{
		for(int cnt2 = hmm_array->n_length2; cnt2 >= 0; cnt2--)
		{
                        bool forbid_STATE_INS1 = false;
                        bool forbid_STATE_INS2 = false;
                        bool forbid_STATE_ALN = false;

                        get_aln_permissions(forcealign,
                                forbid_STATE_ALN,
                                forbid_STATE_INS1,
                                forbid_STATE_INS2,
                                cnt1,
                                cnt2);

			// This loop is for iterating over possible states in this alignment pair.
			for(int current_state = 0; current_state < N_STATES; current_state++)
			{
				//printf("------------------------------------\n");
				
				bool forbid_cur_state = false;

				if(current_state == STATE_ALN && forbid_STATE_ALN)
				{
					forbid_cur_state = true;
				}
				else if(current_state == STATE_INS1 && forbid_STATE_INS1)
                                {
                                        forbid_cur_state = true;
                                }
                                else if(current_state == STATE_INS2 && forbid_STATE_INS2)
                                {
                                        forbid_cur_state = true;
                                }


				hmm_array->probs[cnt1][cnt2][current_state] = xlog(0);

				// This loop is for marginalizing over next state.
				for(int next_state = 0; !forbid_cur_state && next_state < N_STATES; next_state++)
				{
					//double trans_prob = get_log_trans_prob(current_state, next_state);

					// Transition to next state and emission of next pair of symbols.
					if(next_state == STATE_ALN)
					{
						// Next state is emitting starting with cnt1 + 1 and cnt2 + 1.
						//printf("next_state = %s (prob = %f)\n", state_names[next_state], hmm_array->probs[cnt1+1][cnt2+1][next_state]);
						//printf("%s -> %s\n", state_names[current_state], state_names[next_state]);
						double trans_emit_prob = get_trans_emit_prob(current_state, next_state, cnt1+1, cnt2+1);

						//double emit_prob = xlog(1.0/24); // Take emission prob. as 1 for now.

						//printf("Recursing on (%d, %d) with prob. %f.\n", cnt1+1, cnt2+1, hmm_array->probs[cnt1+1][cnt2+1][next_state]);
						//hmm_array->probs[cnt1][cnt2][current_state] = xlog_sum(hmm_array->probs[cnt1][cnt2][current_state], xlog_mul(hmm_array->probs[cnt1+1][cnt2+1][next_state], xlog_mul(trans_prob, emit_prob)));
						hmm_array->probs[cnt1][cnt2][current_state] = xlog_sum(hmm_array->probs[cnt1][cnt2][current_state], xlog_mul(hmm_array->probs[cnt1+1][cnt2+1][next_state], trans_emit_prob));
					}
					else if(next_state == STATE_INS1)
					{
						// Next state is emitting starting with cnt1 + 1 and cnt2.
						//printf("next_state = %s (prob = %f)\n", state_names[`next_state], hmm_array->probs[cnt1+1][cnt2][next_state]);
						//printf("%s -> %s\n", state_names[current_state], state_names[next_state]);
						double trans_emit_prob = get_trans_emit_prob(current_state, next_state, cnt1+1, cnt2);

						//double emit_prob = xlog(1.0/24); // Take emission prob. as 1 for now.
						
						//printf("Recursing on (%d, %d) with prob. %f.\n", cnt1+1, cnt2, hmm_array->probs[cnt1+1][cnt2][next_state]);
						//hmm_array->probs[cnt1][cnt2][current_state] = xlog_sum(hmm_array->probs[cnt1][cnt2][current_state], xlog_mul(hmm_array->probs[cnt1+1][cnt2][next_state], xlog_mul(trans_prob, emit_prob)));
						hmm_array->probs[cnt1][cnt2][current_state] = xlog_sum(hmm_array->probs[cnt1][cnt2][current_state], xlog_mul(hmm_array->probs[cnt1+1][cnt2][next_state], trans_emit_prob));
					}
					else
					{
						// Next state is emitting starting with cnt1 and cnt2 + 1.
						//printf("next_state = %s (prob = %f)\n\n", state_names[next_state], hmm_array->probs[cnt1][cnt2+1][next_state]);
						//printf("%s -> %s\n", state_names[current_state], state_names[next_state]);
						double trans_emit_prob = get_trans_emit_prob(current_state, next_state, cnt1, cnt2+1);

						//double emit_prob = xlog(1.0/24); // Take emission prob. as 1 for now.
						
						//printf("Recursing on (%d, %d) with prob. %f.\n", cnt1, cnt2+1, hmm_array->probs[cnt1+1][cnt2+1][next_state]);
						//hmm_array->probs[cnt1][cnt2][current_state] = xlog_sum(hmm_array->probs[cnt1][cnt2][current_state], xlog_mul(hmm_array->probs[cnt1][cnt2+1][next_state], xlog_mul(trans_prob, emit_prob)));
						hmm_array->probs[cnt1][cnt2][current_state] = xlog_sum(hmm_array->probs[cnt1][cnt2][current_state], xlog_mul(hmm_array->probs[cnt1][cnt2+1][next_state], trans_emit_prob));
					}	
				}

				//printf("%d, %d, current state: %s, prob = %f (%f)\n", cnt1, cnt2, state_names[current_state], hmm_array->probs[cnt1][cnt2][current_state], xexp(hmm_array->probs[cnt1][cnt2][current_state]) );
			}			
		}
	}

	//printf("Backward prob. of sequences is %f\n", hmm_array->probs[0][0][STATE_ALN]);
	//printf("Backward prob. of sequences is %f\n", hmm_array->probs[0][0][STATE_INS1]);
	//printf("Backward prob. of sequences is %f\n", hmm_array->probs[0][0][STATE_INS2]);
	//printf("Backward prob. of sequences is %f\n", xlog_sum(hmm_array->probs[0][0][STATE_INS2], xlog_sum(hmm_array->probs[0][0][STATE_INS1], hmm_array->probs[0][0][STATE_ALN])) );
	//cin.get();
/*
	// Also have to initialize 0th entries of backward array for calculating insertion probabilities correctly.
	for(int cnt = n1; cnt >= 0; cnt)
	{
		if(cnt > 0)
		{
			double trans_emit_prob = get_trans_emit_prob(STATE_INS2, STATE_INS2, 1, cnt+1);

			hmm_array->probs[0][cnt][STATE_ALN] = xlog(0);
			hmm_array->probs[0][cnt][STATE_INS1] = xlog(0);
			hmm_array->probs[0][cnt][STATE_INS2] = xlog_mul(hmm_array->probs[1][cnt+1][STATE_INS2], trans_emit_prob);
		}
		else
		{
			// Skip this case for now.
		}
	}

	for(int cnt = n2; cnt >= 0; cnt)
	{
		if(cnt > 0)
		{
			double trans_emit_prob = get_trans_emit_prob(STATE_INS1, STATE_ALN, cnt+1, 1);

			hmm_array->probs[cnt][0][STATE_ALN] = xlog(0);
			hmm_array->probs[cnt][n2+1][STATE_INS1] = xlog_mul(hmm_array->probs[cnt+1][1][STATE_INS1], trans_emit_prob);
			hmm_array->probs[cnt][0][STATE_INS2] = xlog(0);
		}
		else
		{
			// Skip this case for now.
		}
	}
*/
}

void calculate_ins_prob_array2(t_hmm_array* fore_hmm_array, t_hmm_array* back_hmm_array)
{
	double* ins1_probs = (double*)malloc((fore_hmm_array->n_length1+1) * sizeof(double));
	double* ins2_probs = (double*)malloc((fore_hmm_array->n_length2+1) * sizeof(double));
	double* ind_aln1_probs = (double*)malloc((fore_hmm_array->n_length1+1) * sizeof(double));
	double* ind_aln2_probs = (double*)malloc((fore_hmm_array->n_length2+1) * sizeof(double));

	global_aln_info.ins1_probs = ins1_probs;
	global_aln_info.ins2_probs = ins2_probs;
	global_aln_info.ind_aln1_probs = ind_aln1_probs;
	global_aln_info.ind_aln2_probs = ind_aln2_probs;

	// Calculate insertion probabilities of first sequence.
	for(int cnt1 = 1; cnt1 < fore_hmm_array->n_length1+1; cnt1++)
	{
		ins1_probs[cnt1] = LOG_OF_ZERO;		

		for(int cnt2 = 0; cnt2 < fore_hmm_array->n_length2+1; cnt2++)
		{
			//double cur_ins1_prob = xlog_mul(fore_hmm_array->probs[cnt1][cnt2][STATE_INS1], back_hmm_array->probs[cnt1][cnt2+1][STATE_INS1]);
			double cur_ins1_prob = xlog_mul(fore_hmm_array->probs[cnt1][cnt2][STATE_INS1], back_hmm_array->probs[cnt1][cnt2][STATE_INS1]);

			ins1_probs[cnt1] = xlog_sum(ins1_probs[cnt1], cur_ins1_prob);
		}
	}

	// Calculate insertion probabilities of second sequence.
	for(int cnt2 = 1; cnt2 < fore_hmm_array->n_length2+1; cnt2++)
	{
		ins2_probs[cnt2] = LOG_OF_ZERO;		

		for(int cnt1 = 0; cnt1 < fore_hmm_array->n_length1+1; cnt1++)
		{
			//double cur_ins2_prob = xlog_mul(fore_hmm_array->probs[cnt1][cnt2][STATE_INS2], back_hmm_array->probs[cnt1+1][cnt2][STATE_INS2]);
			double cur_ins2_prob = xlog_mul(fore_hmm_array->probs[cnt1][cnt2][STATE_INS2], back_hmm_array->probs[cnt1][cnt2][STATE_INS2]);

			ins2_probs[cnt2] = xlog_sum(ins2_probs[cnt2], cur_ins2_prob);
		}
	}

	// Calculate individual alignment probabilities of bases in the alignment.
	//for(cnt1 = 1; cnt1 < fore_hmm_array->n_length1+1; cnt1++)
	for(int cnt1 = 1; cnt1 < fore_hmm_array->n_length1+1; cnt1++)
	{
		ind_aln1_probs[cnt1] = LOG_OF_ZERO;

		for(int cnt2 = 1; cnt2 < fore_hmm_array->n_length2+1; cnt2++)
		{
			double cur_aln1_prob = xlog_mul(fore_hmm_array->probs[cnt1][cnt2][STATE_ALN], back_hmm_array->probs[cnt1][cnt2][STATE_ALN]);

			ind_aln1_probs[cnt1] = xlog_sum(ind_aln1_probs[cnt1], cur_aln1_prob);
		}
	}

	for(int cnt2 = 1; cnt2 < fore_hmm_array->n_length2+1; cnt2++)
	//for(cnt2 = 1; cnt2 < fore_hmm_array->n_length2+1; cnt2++)
	{
		ind_aln2_probs[cnt2] = LOG_OF_ZERO;

		for(int cnt1 = 1; cnt1 < fore_hmm_array->n_length1+1; cnt1++)
		{
			double cur_aln2_prob = xlog_mul(fore_hmm_array->probs[cnt1][cnt2][STATE_ALN], back_hmm_array->probs[cnt1][cnt2][STATE_ALN]);

			ind_aln2_probs[cnt2] = xlog_sum(ind_aln2_probs[cnt2], cur_aln2_prob);
		}
	}
}

