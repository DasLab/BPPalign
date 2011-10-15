#include <cstdlib>
#include <stdio.h>
#include "hmm_arrays.h"
#include "misc_decoding_services.h"
#include "hmm_probs.h"
#include "alignment_hmm_model.h"
#include "xlog_math.h"
#include "process_output.h"
//#include "annot_aln_services.h"
#include <iostream>
#include <math.h>

//#include <conio.h>

// Do the posterior decoding based on taking max. prob. path on the next class of steps
// It should be noted that I cannot apply a fully posterior decoding to this problem,
// because output indexing is dependent on states.
// My approach will not be a full posterior decoding but it is worth the try for seeing
// the results and comparing it with actual alignment.
extern t_aln_info global_aln_info;

extern double emit_probs[N_OUTPUTS][N_STATES];
extern double trans_probs[N_STATES][N_STATES];

char* ml_aln_line1;
char* ml_aln_line2;

// Output an alignment file.
void decode_posterior()
{
	// Following are for storing alignment data.
	char post_aln_line1[500];
	char post_aln_line2[500];

	post_aln_line1[0] = 0;
	post_aln_line2[0] = 0;

	// Start from (0,0) and find next most probable alignment pair, including 
	// insertions.

	int cnt1 = 1;
	int cnt2 = 1;

	while(cnt1 != global_aln_info.length1+1 && cnt2 != global_aln_info.length2+1)
	{
		// Find the most probable state based on current output, write to alignment strings.
		// ALN, INS1 or INS2 is most probable.
		int opt_state = STATE_ALN;
		double opt_prob = global_aln_info.aln_probs[cnt1][cnt2];

		double cur_ins1_prob = xlog_mul(global_aln_info.fore_hmm_array->probs[cnt1][cnt2][STATE_INS1], global_aln_info.back_hmm_array->probs[cnt1][cnt2+1][STATE_INS1]);
		double cur_ins2_prob = xlog_mul(global_aln_info.fore_hmm_array->probs[cnt1][cnt2][STATE_INS2], global_aln_info.back_hmm_array->probs[cnt1+1][cnt2][STATE_INS2]);

		if(cur_ins1_prob > opt_prob)
		{
			opt_state = STATE_INS1;
			opt_prob = cur_ins1_prob;
		}

		if(cur_ins2_prob > opt_prob)
		{
			opt_state = STATE_INS2;
			opt_prob = cur_ins2_prob;
		}

		// I have optimal prob. and optimal state for this case. Now update alignment strings.
		switch(opt_state)
		{
		case(STATE_ALN):
			sprintf(post_aln_line1, "%s%c", post_aln_line1, get_nuc(1, cnt1-1));
			sprintf(post_aln_line2, "%s%c", post_aln_line2, get_nuc(2, cnt2-1));
			cnt1++;
			cnt2++;
			break;

		case (STATE_INS1):
			sprintf(post_aln_line1, "%s%c", post_aln_line1, get_nuc(1, cnt1-1));
			sprintf(post_aln_line2, "%s-", post_aln_line2);
			cnt1++;
			break;

		case(STATE_INS2):
			sprintf(post_aln_line1, "%s-", post_aln_line1);
			sprintf(post_aln_line2, "%s%c", post_aln_line2, get_nuc(2, cnt2-1));
			cnt2++;
			break;
		}

		printf("Current strings (%d, %d):\n%s\n%s\n", cnt1, cnt2, post_aln_line1, post_aln_line2);
	}

	// Add last parts to alignment strings.
	if(cnt1 != global_aln_info.length1 + 1)
	{
		while(cnt1 < global_aln_info.length1 + 1)
		{
			sprintf(post_aln_line1, "%s%c", post_aln_line1, get_nuc(1, cnt1-1));
			sprintf(post_aln_line2, "%s-", post_aln_line2);
			cnt1++;
		}
	}

	if(cnt2 != global_aln_info.length2 + 1)
	{
		while(cnt2 < global_aln_info.length2 + 1)
		{
			sprintf(post_aln_line1, "%s-", post_aln_line1);
			sprintf(post_aln_line2, "%s%c", post_aln_line2, get_nuc(2, cnt2-1));
			cnt2++;
		}
	}

	printf("Posterior decoded alignment:\n%s\n%s\n", post_aln_line1, post_aln_line2);

	FILE* aln_file = fopen("post_aln.txt", "w");
	fprintf(aln_file, "%s\n%s\n", post_aln_line1, post_aln_line2);
	fclose(aln_file);
}

void get_windowed_ins_probs(int win_length)
{
	// Allocate two arrays for ins. prob over each window of length win_length over each sequence.
	double* winned_ins_probs1 = (double*)malloc(sizeof(double) * global_aln_info.length1 - win_length + 1);
	double* winned_ins_probs2 = (double*)malloc(sizeof(double) * global_aln_info.length2 - win_length + 1);

	char win_ins_prob_file_name[100];
	sprintf(win_ins_prob_file_name, "win_ins_probs1_%d.txt", win_length);
	FILE* win_ins_prob_file = fopen(win_ins_prob_file_name, "w");

	for(int cnt1 = 1; cnt1 < global_aln_info.length1 - win_length + 1; cnt1++)
	{
		// cnt1 is the index of starting of current window.
		winned_ins_probs1[cnt1] = LOG_OF_ZERO;

		// cnt2 is the index in second sequence where window is betweenm i.e., 
		// the window is assumed to be inserted between cnt2 and cnt2 + 1.
		// however I do not put the constraint that cnt1-1 is not aligned to cnt2 and cnt2 is not aligned to 
		// cnt1 + win_length + 1.
		for(int cnt2 = 1; cnt2 < global_aln_info.length2 + 1; cnt2++)
		{
			// cur_win_ins_prob includes probability of all enumerations of alignments where 
			// currently set window is inserted in first sequences.
			double cur_win_ins_prob	= xlog_mul(global_aln_info.fore_hmm_array->probs[cnt1][cnt2][STATE_INS1], global_aln_info.back_hmm_array->probs[cnt1 + win_length][cnt2 + 1][STATE_INS1]);			

			// Calculate the emission probability of inserted sequence by going over all inserted sequence.
			// emit the nucleotides in the inserted sequence is sequence 1.
			double int_emit_prob = xlog(1);

			for(int emit_cnt = cnt1 + 1; emit_cnt <= cnt1 + win_length; emit_cnt++)
			{
				double trans_emit_prob;

				if(emit_cnt == cnt1 + win_length)
				{
					trans_emit_prob = get_trans_emit_prob(STATE_INS1, STATE_INS1, emit_cnt, cnt2);
				}
				else
				{
					trans_emit_prob = get_trans_emit_prob(STATE_INS1, STATE_INS1, emit_cnt, cnt2);
				}

				int_emit_prob = xlog_mul(int_emit_prob, trans_emit_prob);
			}

			cur_win_ins_prob = xlog_mul(cur_win_ins_prob, int_emit_prob);

			winned_ins_probs1[cnt1] = xlog_sum(winned_ins_probs1[cnt1], cur_win_ins_prob);			
		}	

		fprintf(win_ins_prob_file, "%f ", xlog_div(winned_ins_probs1[cnt1], global_aln_info.op_prob));
	}

	fprintf(win_ins_prob_file, "\n");

	fclose(win_ins_prob_file);

/*
	for(int cnt2 = 1; cnt2 < global_aln_info.length2 - win_length + 1; cnt2++)
	{
		winned_ins_probs2[cnt2] = LOG_OF_ZERO;

		for(int cnt1 = 1; cnt1 < global_aln_info.length1 + 1; cnt1++)
		{
			double cur_win_ins_prob	= global_aln_info.fore_hmm_array->probs[cnt1][cnt2][STATE_ALN] * global_aln_info.back_hmm_array->probs[cnt1 + 1][cnt2 + win_length + 1][ALN];
			winned_ins_probs2[cnt2] = xlog_sum(winned_ins_probs2[cnt2], cur_win_ins_prob);
		}
	}
*/
}

void get_windowed_ins_probs2(int win_length)
{
	// Allocate two arrays for ins. prob over each window of length win_length over each sequence.
	double* winned_ins_probs1 = (double*)malloc(sizeof(double) * global_aln_info.length1 - win_length + 1);
	double* winned_ins_probs2 = (double*)malloc(sizeof(double) * global_aln_info.length2 - win_length + 1);

	char win_ins_prob_file_name[100];
	sprintf(win_ins_prob_file_name, "win_ins_probs1_%d.txt", win_length);
	FILE* win_ins_prob_file = fopen(win_ins_prob_file_name, "w");

	for(int cnt1 = 1; cnt1 < global_aln_info.length1 - win_length + 1; cnt1++)
	{
		// cnt1 is the index of starting of current window.
		winned_ins_probs1[cnt1] = LOG_OF_ZERO;

		// cnt2 is the index in second sequence where window is betweenm i.e., 
		// the window is assumed to be inserted between cnt2 and cnt2 + 1.
		// however I do not put the constraint that cnt1-1 is not aligned to cnt2 and cnt2 is not aligned to 
		// cnt1 + win_length + 1.
		for(int cnt2 = 1; cnt2 < global_aln_info.length2 + 1; cnt2++)
		{
			// cur_win_ins_prob includes probability of all enumerations of alignments where 
			// currently set window is inserted in first sequences.
			double cur_win_ins_prob	= xlog_mul(global_aln_info.fore_hmm_array->probs[cnt1][cnt2][STATE_INS1], global_aln_info.back_hmm_array->probs[cnt1 + win_length][cnt2][STATE_INS1]);			

			// Calculate the emission probability of inserted sequence by going over all inserted sequence.
			// emit the nucleotides in the inserted sequence is sequence 1.
			double int_emit_prob = xlog(1);

			for(int emit_cnt = cnt1 + 1; emit_cnt <= cnt1 + win_length; emit_cnt++)
			{
				double trans_emit_prob;

				if(emit_cnt == cnt1 + win_length)
				{
					trans_emit_prob = get_trans_emit_prob(STATE_INS1, STATE_INS1, emit_cnt, cnt2);
				}
				else
				{
					trans_emit_prob = get_trans_emit_prob(STATE_INS1, STATE_INS1, emit_cnt, cnt2);
				}

				int_emit_prob = xlog_mul(int_emit_prob, trans_emit_prob);
			}

			cur_win_ins_prob = xlog_mul(cur_win_ins_prob, int_emit_prob);

			winned_ins_probs1[cnt1] = xlog_sum(winned_ins_probs1[cnt1], cur_win_ins_prob);			
		}	

		fprintf(win_ins_prob_file, "%f ", xlog_div(winned_ins_probs1[cnt1], global_aln_info.op_prob));
	}

	fprintf(win_ins_prob_file, "\n");

	fclose(win_ins_prob_file);
}

// Apply a viterbi like algorithm to forward array to determine ML alignment to estimate alignment similarity
// of 
void decode_ML(short** forcealign)
{
	//printf("Setting model parameters as ML_ parameters.\n");
	set_hmm_parameters(ML_emit_probs, ML_trans_probs);

	//printf("\n\nML decoding alignment...\n");

	double*** max_score_array;
	char*** ML_path_array;
	int cnt1, cnt2;

	max_score_array = (double***)malloc(sizeof(double**) * (global_aln_info.length1 + 3));
	// Allocate and init. max scores matrix.
	for(cnt1 = 0; cnt1 <= global_aln_info.length1 + 1; cnt1++)
	{
		max_score_array[cnt1] = (double**)malloc(sizeof(double*) * (global_aln_info.length2 + 3));

		for(cnt2 = 0; cnt2 <= global_aln_info.length2 + 1; cnt2++)
		{
			max_score_array[cnt1][cnt2] = (double*)malloc(sizeof(double) * N_STATES);

			// Initialize starting state max score as 1. 
			// Note following is done over and over again which is ok from implementation point of view 
			// to simplify init.ing code.
			max_score_array[0][0][STATE_ALN] = xlog(1);
			max_score_array[0][0][STATE_INS1] = xlog(0);
			max_score_array[0][0][STATE_INS2] = xlog(0);

			// init max scores.
			if(cnt1 == 0 && cnt2 > 0)
			{
				max_score_array[cnt1][cnt2][STATE_ALN] = xlog(0);
				max_score_array[cnt1][cnt2][STATE_INS1] = xlog(0); 
				max_score_array[cnt1][cnt2][STATE_INS2] = xlog(0); 
				if(cnt2 > 1)
				{
					max_score_array[cnt1][cnt2][STATE_INS2] = xlog_mul(max_score_array[cnt1][cnt2 - 1][STATE_INS2], get_trans_emit_prob(STATE_INS2, STATE_INS2, cnt1, cnt2));
				}
				else
				{
					max_score_array[cnt1][cnt2][STATE_INS2] = xlog_mul(max_score_array[cnt1][cnt2 - 1][STATE_ALN], get_trans_emit_prob(STATE_ALN, STATE_INS2, cnt1, cnt2));					
				}
			}
			else if(cnt1 > 0 && cnt2 == 0)
			{
				max_score_array[cnt1][cnt2][STATE_ALN] = xlog(0);
				max_score_array[cnt1][cnt2][STATE_INS1] = xlog(0); 
				max_score_array[cnt1][cnt2][STATE_INS2] = xlog(0); 
				if(cnt1 > 1)
				{
					max_score_array[cnt1][cnt2][STATE_INS1] = xlog_mul(max_score_array[cnt1 - 1][cnt2][STATE_INS1], get_trans_emit_prob(STATE_INS1, STATE_INS1, cnt1, cnt2));
				}
				else
				{
					max_score_array[cnt1][cnt2][STATE_INS1] = xlog_mul(max_score_array[cnt1 - 1][cnt2][STATE_ALN], get_trans_emit_prob(STATE_ALN, STATE_INS1, cnt1, cnt2));
				}
			}
			else
			{
				max_score_array[cnt1][cnt2][STATE_ALN] = xlog(0);
				max_score_array[cnt1][cnt2][STATE_INS1] = xlog(0); 
				max_score_array[cnt1][cnt2][STATE_INS2] = xlog(0); 
			}
		}
	}

	ML_path_array = (char***)malloc(sizeof(char**) * (global_aln_info.length1 + 3));
	// Allocate and init. max scoring path matrix for traceback.
	for(cnt1 = 0; cnt1 < global_aln_info.length1 + 2; cnt1++)
	{
		ML_path_array[cnt1] = (char**)malloc(sizeof(char*) * (global_aln_info.length2 + 3));
		for(cnt2 = 0; cnt2 < global_aln_info.length2 + 2; cnt2++)
		{
			ML_path_array[cnt1][cnt2] = (char*)malloc(sizeof(char) * N_STATES);

			ML_path_array[0][0][STATE_ALN] = -1; // Beginning of everything, origin of HMM.
			ML_path_array[0][0][STATE_INS1] = -2;
			ML_path_array[0][0][STATE_INS2] = -2;

			if(cnt1 == 0 && cnt2 > 0)
			{
				ML_path_array[cnt1][cnt2][STATE_ALN] = -2;
				ML_path_array[cnt1][cnt2][STATE_INS1] = -2; 
				ML_path_array[cnt1][cnt2][STATE_INS2] = -2; 

				if(cnt2 > 1)
				{
					ML_path_array[cnt1][cnt2][STATE_INS2] = STATE_INS2;
				}
				else
				{
					ML_path_array[cnt1][cnt2][STATE_INS2] = STATE_ALN;
				}
			}
			else if(cnt1 > 0 || cnt2 == 0)
			{
				ML_path_array[cnt1][cnt2][STATE_ALN] = -2;
				ML_path_array[cnt1][cnt2][STATE_INS1] = -2; 
				ML_path_array[cnt1][cnt2][STATE_INS2] = -2; 

				if(cnt1 > 1)
				{
					ML_path_array[cnt1][cnt2][STATE_INS1] = STATE_INS1;
				}
				else
				{
					ML_path_array[cnt1][cnt2][STATE_INS1] = STATE_ALN;
				}
			}
			else
			{
				ML_path_array[cnt1][cnt2][STATE_ALN] = -2;
				ML_path_array[cnt1][cnt2][STATE_INS1] = -2; 
				ML_path_array[cnt1][cnt2][STATE_INS2] = -2; 
			}
		}
	}

	// Calculate max scores matrix and set max scoring path all along.
	for(cnt1 = 1; cnt1 <= global_aln_info.length1 + 1; cnt1++)
	{
		for(cnt2 = 1; cnt2 <= global_aln_info.length2 + 1; cnt2++)
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

			//printf("(%d, %d)\r", cnt1, cnt2);
			// Go over cur states.
			for(int cur_state = 0; cur_state < N_STATES; cur_state++)
			{
				// Choose max path through next states using transition and emission of next symbol.

				double trans_emit_prob = 0;

				double max_score = xlog(0);
				char max_scoring_state = 0;

				for(int prev_state = 0; prev_state < N_STATES; prev_state++)
				{
					if(!forbid_STATE_ALN &&
						cur_state == STATE_ALN)
					{
						//printf("trans_emit = %.3f\n", trans_emit_prob);
						trans_emit_prob = get_trans_emit_prob(prev_state, cur_state, cnt1, cnt2);

						// Set max, that is min for probabilities.
						if(xlog_mul(max_score_array[cnt1 - 1][cnt2 - 1][prev_state], trans_emit_prob) > max_score)
						{
							//max_score_array[cnt1][cnt2] = xlog_mul(max_score_array[cnt1][cnt2], trans_emit_prob);
							max_score = xlog_mul(max_score_array[cnt1 - 1][cnt2 - 1][prev_state], trans_emit_prob);

							max_scoring_state = prev_state;
						}
					}
					else if(!forbid_STATE_INS1 &&
						cur_state == STATE_INS1)
					{
						//printf("trans_emit = %.3f\n", trans_emit_prob);

						trans_emit_prob = get_trans_emit_prob(prev_state, cur_state, cnt1, cnt2);

						// Set max, that is min for probabilities.
						if(xlog_mul(max_score_array[cnt1 - 1][cnt2][prev_state], trans_emit_prob) > max_score)
						{
							//max_score_array[cnt1][cnt2] = xlog_mul(max_score_array[cnt1][cnt2], trans_emit_prob);
							max_score = xlog_mul(max_score_array[cnt1 - 1][cnt2][prev_state], trans_emit_prob);

							max_scoring_state = prev_state;
						}
					}
					else if(!forbid_STATE_INS2 &&
						cur_state == STATE_INS2)
					{

						trans_emit_prob = get_trans_emit_prob(prev_state, cur_state, cnt1, cnt2);
						//printf("trans_emit = %.3f\n", trans_emit_prob);

						// Set max, that is min for probabilities.
						if(xlog_mul(max_score_array[cnt1][cnt2 - 1][prev_state], trans_emit_prob) > max_score)
						{
							//max_score_array[cnt1][cnt2] = xlog_mul(max_score_array[cnt1][cnt2], trans_emit_prob);
							max_score = xlog_mul(max_score_array[cnt1][cnt2 - 1][prev_state], trans_emit_prob);

							max_scoring_state = prev_state;
						}
					}
				} // next_state

				// Copy max score.
				max_score_array[cnt1][cnt2][cur_state] = max_score;
				//printf("(%d, %d): STATE = %d, max. score = %.3f\n", cnt1, cnt2, cur_state, max_score);

				// Set max scoring path for current state.
				ML_path_array[cnt1][cnt2][cur_state] = max_scoring_state;

			} // current_state
		} // cnt2
	} // cnt1

	// Calculate score for ending sequence and set ML path for ending sequence.
	double max_score = xlog(0);
	char max_scoring_state = 0;
	for(int cnt = 0; cnt < N_STATES; cnt++)
	{
		double log_trans_prob = get_log_trans_prob(cnt, STATE_ALN);
	
		if(xlog_mul(max_score_array[global_aln_info.length1][global_aln_info.length2][cnt], log_trans_prob) > max_score)
		{
			 max_score = xlog_mul(max_score_array[global_aln_info.length1][global_aln_info.length2][cnt], log_trans_prob);

			 max_scoring_state = cnt;
		}
	}

	//printf("Max ML prob = %f\n", max_score);

	ML_path_array[global_aln_info.length1 + 1][global_aln_info.length2 + 1][STATE_ALN] = max_scoring_state;

	// Do traceback over max scoring path to get ML alignment.
	//printf("Tracing back ML path of seq length %d, %d\n", global_aln_info.length1, global_aln_info.length2);
	cnt1 = global_aln_info.length1 + 1;
	cnt2 = global_aln_info.length2 + 1;
	char cur_state = STATE_ALN;

	int ins1_length = 0; // This is for determining length of alignment strings.
	int ins2_length = 0; // This is for determining length of alignment strings.

	// cur_state is the current state that traceback is in.
	while(cur_state != -1)
	{
		switch(cur_state)
		{
		case STATE_ALN:
			cur_state = ML_path_array[cnt1][cnt2][STATE_ALN];
			cnt1--;
			cnt2--;
			break;

		case STATE_INS1:
			ins1_length++;
			cur_state = ML_path_array[cnt1][cnt2][STATE_INS1];
			cnt1--;

			break;

		case STATE_INS2:
			ins2_length++;
			cur_state = ML_path_array[cnt1][cnt2][STATE_INS2];
			cnt2--;			
			break;
		}

		//printf("%s (%d, %d)\n", state_names[cur_state], cnt1 , cnt2);
	}

	int len_aln_line = (global_aln_info.length1  + ins1_length +  global_aln_info.length2 + ins2_length) / 2;
	//printf("Length of alignment strings is %d with %d ins1 and %d ins2.\n", len_aln_line, ins1_length, ins2_length);

	// Dump ML alignment in appropriate format.

	// Allocate alignment strings.
	ml_aln_line1 = (char*)malloc(sizeof(char) * (len_aln_line + 2));
	ml_aln_line2 = (char*)malloc(sizeof(char) * (len_aln_line + 2));
	ml_aln_line1[len_aln_line] = 0;
	ml_aln_line2[len_aln_line] = 0;
	int aln_line_index = len_aln_line - 1;

	// Do traceback again and copy strings.
	//printf("RE-Tracing back for copying alignment strings for seq. of length %d %d.\n", global_aln_info.length1, global_aln_info.length2);
	cnt1 = global_aln_info.length1 + 1;
	cnt2 = global_aln_info.length2 + 1;
	cur_state = STATE_ALN;

	// Switch to emitting states to correctly log alignment.
	switch(cur_state)
	{
	case STATE_ALN:
		cur_state = ML_path_array[cnt1][cnt2][STATE_ALN];
		cnt1--;
		cnt2--;
		break;

	case STATE_INS1:
		cur_state = ML_path_array[cnt1][cnt2][STATE_INS1];
		cnt1--;
		break;

	case STATE_INS2:
		cur_state = ML_path_array[cnt1][cnt2][STATE_INS2];
		cnt2--;		
		break;
	}

	// cur_state is the current state that traceback is in.
	// Copy one character at the end of each of alignment lines.
	while(cur_state != -1 && aln_line_index >= 0)
	{
		switch(cur_state)
		{
		case STATE_ALN:
			if(aln_line_index < 0)
			{
				printf("PROBLEM !!!, %d\n", cur_state);
			}

			ml_aln_line1[aln_line_index] = get_nuc(1, cnt1 - 1);
			ml_aln_line2[aln_line_index] = get_nuc(2, cnt2 - 1);

			//printf("%s (%d, %d), %d(%c, %c)\n", state_names[cur_state], cnt1 , cnt2, aln_line_index, get_nuc(1, cnt1 - 1), get_nuc(2, cnt2 - 1));

			cur_state = ML_path_array[cnt1][cnt2][STATE_ALN];
			cnt1--;
			cnt2--;

			break;

		case STATE_INS1:
			if(aln_line_index < 0)
			{
				printf("PROBLEM !!!, %d\n", cur_state);
			}

			ml_aln_line1[aln_line_index] = get_nuc(1, cnt1 - 1);
			ml_aln_line2[aln_line_index] = '.';

			//printf("%s (%d, %d), %d(%c, %c)\n", state_names[cur_state], cnt1 , cnt2, aln_line_index, get_nuc(1, cnt1 - 1), '.');

			cur_state = ML_path_array[cnt1][cnt2][STATE_INS1];
			cnt1--;

			break;

		case STATE_INS2:
			if(aln_line_index < 0)
			{
				printf("PROBLEM !!!, %d\n", cur_state);
			}

			ml_aln_line1[aln_line_index] = '.';
			ml_aln_line2[aln_line_index] = get_nuc(2, cnt2 - 1);

			//printf("%s (%d, %d), %d(%c, %c)\n", state_names[cur_state], cnt1 , cnt2, aln_line_index, '.', get_nuc(2, cnt2 - 1));

			cur_state = ML_path_array[cnt1][cnt2][STATE_INS2];
			cnt2--;		

			break;
		}

		//printf("Index = %d, cur_state = %d\n", aln_line_index, cur_state);

		aln_line_index--;
	}

	// Check for error.
	if(aln_line_index != -1 || cnt1 != 0 || cnt2 != 0)
	{
		printf("Indexing traceback problem at %s(%d)\n", __FILE__, __LINE__);

		exit(0);
	}

	//printf("%s\n%s\n", ml_aln_line1, ml_aln_line2);;

	// Dump alignment strings.
	//FILE* aln_file = fopen("ml_alignment.txt", "w");
	//fprintf(aln_file, "%s\n%s\n", ml_aln_line1, ml_aln_line2);
	//fclose(aln_file);

	// Free array memory.
	free_ML_arrays(max_score_array, ML_path_array);
}

void free_ML_arrays(double*** score_array, char*** ML_path_array)
{
	// Allocate and init. max scores matrix.
	for(int cnt1 = 0; cnt1 <= global_aln_info.length1 + 1; cnt1++)
	{
		//max_score_array[cnt1] = (double**)malloc(sizeof(double*) * (global_aln_info.length2 + 1));
		for(int cnt2 = 0; cnt2 <= global_aln_info.length2 + 1; cnt2++)
		{
			//max_score_array[cnt1][cnt2] = (double*)malloc(sizeof(double) * N_STATES);

			free(score_array[cnt1][cnt2]);
		}

		free(score_array[cnt1]);
	}

	free(score_array);

	// Allocate and init. max scoring path matrix for traceback.
	for(int cnt1 = 0; cnt1 < global_aln_info.length1 + 2; cnt1++)
	{
		//ML_path_array[cnt1] = (char**)malloc(sizeof(char*) * (global_aln_info.length2 + 3));
		for(int cnt2 = 0; cnt2 < global_aln_info.length2 + 2; cnt2++)
		{
			//ML_path_array[cnt1][cnt2] = (char*)malloc(sizeof(char) * N_STATES);

			free(ML_path_array[cnt1][cnt2]);
		}

		free(ML_path_array[cnt1]);
	}

	free(ML_path_array);


}
