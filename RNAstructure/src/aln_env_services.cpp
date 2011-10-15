#include "aln_env_services.h"
#include <stdio.h>
#include "hmm_arrays.h"
#include "misc_decoding_services.h"
#include "hmm_probs.h"
#include "alignment_hmm_model.h"
#include "xlog_math.h"
#include "process_output.h"
#include "annot_aln_services.h"
#include <iostream>
#include <math.h>

extern t_aln_info global_aln_info;

// Use alignment probabilities to calculate an alignment envelope to input to dynalign.
// This function also calculates a value for M whenever the alignment envelope need not be calculated
// that is an M parameter rather than an alignment envelope is to be used, if an alignment
// envelope is used, M can be set to max value posibble which is length of second sequence.
void get_alignment_envelopes(double threshold_prob, bool** aln_env, int* M)
{
	if(aln_env == NULL)
	{
		aln_env = (bool**)malloc((global_aln_info.fore_hmm_array->n_length1+2) * sizeof(bool*));

		for(int cnt = 0; cnt < global_aln_info.fore_hmm_array->n_length1+2; cnt++)
		{
			aln_env[cnt] = (bool*)malloc((global_aln_info.fore_hmm_array->n_length2+1) * sizeof(bool));
		}		
	}

	double log_threshold = xlog(threshold_prob);

	int last1 = 1;
	int last2 = 1;

	for(int cnt1 = 1; cnt1 < global_aln_info.fore_hmm_array->n_length1+1; cnt1++)
	{
		for(int cnt2 = 1; cnt2 < global_aln_info.fore_hmm_array->n_length2+1; cnt2++)
		{
			//printf("(%d, %d): %f, %f\n", cnt1, cnt2, xlog_div(global_aln_info.aln_probs[cnt1][cnt2], global_aln_info.op_prob), log_threshold);
			double ins1_prob = xlog_mul(global_aln_info.fore_hmm_array->probs[cnt1][cnt2][STATE_INS1], global_aln_info.back_hmm_array->probs[cnt1][cnt2][STATE_INS1]);
			double ins2_prob = xlog_mul(global_aln_info.fore_hmm_array->probs[cnt1][cnt2][STATE_INS2], global_aln_info.back_hmm_array->probs[cnt1][cnt2][STATE_INS2]);
			double three_plane_sum = xlog_sum(global_aln_info.aln_probs[cnt1][cnt2], xlog_sum(ins1_prob, ins2_prob));

			//if( xlog_div(global_aln_info.aln_probs[cnt1][cnt2], global_aln_info.op_prob) < log_threshold )
			if( xlog_div(three_plane_sum, global_aln_info.op_prob) < log_threshold )	
			{
				aln_env[cnt1][cnt2] = false;
			}
			else
			{
				aln_env[cnt1][cnt2] = true;
			}
		}
	}

	// Calculate pruned alignment envelope and set it to global_aln_info's alignment envelope, 
	// calculate also the size of alignment envelope.
	#define _PRUNE_ALN_
	#ifdef _PRUNE_ALN_
		prune_aln_env(aln_env);
	#else
		copy_aln_env(aln_env);
	#endif

	FILE* unpruned_aln_file = fopen("unpruned_aln.txt", "w");

	// Dump unpruned alignment envelope.
	for(int cnt1 = 1; cnt1 < global_aln_info.fore_hmm_array->n_length1+1; cnt1++)
	{
		for(int cnt2 = 1; cnt2 < global_aln_info.fore_hmm_array->n_length2+1; cnt2++)
		{
			fprintf(unpruned_aln_file, "%d ", aln_env[cnt1][cnt2]);
		}

		fprintf(unpruned_aln_file, "\n");
	}

	fclose(unpruned_aln_file);

	// Free unpruned aln. env. since it is of no use any more.
	if(aln_env != NULL)
	{
		for(int cnt = 1; cnt < global_aln_info.fore_hmm_array->n_length1+1; cnt++)
		{
			free(aln_env[cnt] );
		}	

		free(aln_env);	
	}

	// Calculate smallest M constraint that includes probabilistic alignment constraint.
	global_aln_info.smallest_M = M;

	*M = 1;
	calculate_M_from_aln_env(M);

	//printf("M from prob. aln. env. is %d\n", *M);
	////////////////////////////////////////////////////////////////////////////
}

void prune_aln_env(bool** aln_env)
{
	if(global_aln_info.corr_aln_env == NULL)
	{
		global_aln_info.corr_aln_env = (char**)malloc((global_aln_info.fore_hmm_array->n_length1+1) * sizeof(char*));
		global_aln_info.dump_aln_env = (char**)malloc((global_aln_info.fore_hmm_array->n_length1+1) * sizeof(char*));

		if(global_aln_info.corr_aln_env == NULL || global_aln_info.dump_aln_env == NULL)
		{
			printf("NULL ptr @ %s(%d)\n", __FILE__, __LINE__);
			exit(0);
		}

		for(int cnt = 0; cnt < global_aln_info.fore_hmm_array->n_length1+1; cnt++)
		{
			global_aln_info.corr_aln_env[cnt] = (char*)malloc((global_aln_info.fore_hmm_array->n_length2+1) * sizeof(char));
			global_aln_info.dump_aln_env[cnt] = (char*)malloc((global_aln_info.fore_hmm_array->n_length2+1) * sizeof(char));

			if(global_aln_info.corr_aln_env[cnt] == NULL || global_aln_info.dump_aln_env[cnt] == NULL)
			{
				printf("NULL ptr @ %s(%d)\n", __FILE__, __LINE__);
				exit(0);
			}
		}
	}

	// Size of alignment constraint.
	int n_aln_env = 0;

	// Find the first 1 in the first row of aln_env or column.
	for(int cnt1 = 1; cnt1 < global_aln_info.fore_hmm_array->n_length1+1; cnt1++)
	{
		for(int cnt2 = 1; cnt2 < global_aln_info.fore_hmm_array->n_length2+1; cnt2++)
		{
			if(check_connection(cnt1, cnt2, aln_env, global_aln_info.corr_aln_env))
			{
				global_aln_info.corr_aln_env[cnt1][cnt2] = ALN_ACK;
				global_aln_info.dump_aln_env[cnt1][cnt2] = ALN_ACK;
			}
			else
			{
				global_aln_info.corr_aln_env[cnt1][cnt2] = ALN_NACK;
				global_aln_info.dump_aln_env[cnt1][cnt2] = ALN_NACK;
			}

		}
	}

	for(int cnt1 = global_aln_info.fore_hmm_array->n_length1; cnt1 >= 1; cnt1--)
	{
		for(int cnt2 = global_aln_info.fore_hmm_array->n_length2; cnt2 >= 1; cnt2--)
		{
			if(check_forward_connection(cnt1, cnt2, global_aln_info.corr_aln_env, global_aln_info.corr_aln_env))
			{
				global_aln_info.corr_aln_env[cnt1][cnt2] = ALN_ACK;
				global_aln_info.dump_aln_env[cnt1][cnt2] = ALN_ACK;
				n_aln_env++;
			}
			else
			{
				global_aln_info.corr_aln_env[cnt1][cnt2] = ALN_NACK;
				global_aln_info.dump_aln_env[cnt1][cnt2] = ALN_NACK;
			}
		}
	}

	// Calculate size of probabilistic alignment envelope.
	//printf("Size of aln. env. is %d\n", n_aln_env);
	global_aln_info.n_aln_env = n_aln_env;
}

// Copies and envelope to correct envelope of golobal info, this is an alternative to pruning.
void copy_aln_env(bool** aln_env)
{
	if(global_aln_info.corr_aln_env == NULL)
	{
		global_aln_info.corr_aln_env = (char**)malloc((global_aln_info.fore_hmm_array->n_length1+1) * sizeof(char*));
		global_aln_info.dump_aln_env = (char**)malloc((global_aln_info.fore_hmm_array->n_length1+1) * sizeof(char*));

		for(int cnt = 0; cnt < global_aln_info.fore_hmm_array->n_length1+1; cnt++)
		{
			global_aln_info.corr_aln_env[cnt] = (char*)malloc((global_aln_info.fore_hmm_array->n_length2+1) * sizeof(char));
			global_aln_info.dump_aln_env[cnt] = (char*)malloc((global_aln_info.fore_hmm_array->n_length2+1) * sizeof(char));
		}
	}

	int n_aln_env = 0;
	for(int cnt1 = global_aln_info.fore_hmm_array->n_length1; cnt1 >= 1; cnt1--)
	{
		for(int cnt2 = global_aln_info.fore_hmm_array->n_length2; cnt2 >= 1; cnt2--)
		{
			if(aln_env[cnt1][cnt2])
			{
				n_aln_env++;
				global_aln_info.corr_aln_env[cnt1][cnt2] = ALN_ACK;
				global_aln_info.dump_aln_env[cnt1][cnt2] = ALN_ACK;
			}
			else
			{
				global_aln_info.corr_aln_env[cnt1][cnt2] = ALN_NACK;
				global_aln_info.dump_aln_env[cnt1][cnt2] = ALN_NACK;
			}
		}
	}

	//printf("Size of aln. env. is %d\n", n_aln_env);
	global_aln_info.n_aln_env = n_aln_env;	
}

// Checks the connection between two points in aln_env, checks immediate neighbors.
bool check_connection(int cnt1, int cnt2, bool** aln_env, char** corr_aln_env)
{
	if(aln_env[cnt1][cnt2])
	{
		// If at the first row or first column, just return true.
		if(cnt1 == 1 || cnt2 == 1)
		{
			return(true);
		}
		else
		{
			if(corr_aln_env[cnt1-1][cnt2] == ALN_ACK || corr_aln_env[cnt1][cnt2-1] == ALN_ACK || corr_aln_env[cnt1-1][cnt2-1] == ALN_ACK)
			{
				return(true);
			}
			else
			{
				return(false);
			}
		}
	}
	else
	{
		return(false);
	}
}


// Checks the connection between two points in aln_env, checks immediate neighbors.
bool check_forward_connection(int cnt1, int cnt2, char** pruned_aln_env, char** corr_aln_env)
{
	if(pruned_aln_env[cnt1][cnt2] == ALN_ACK)
	{
		// If at the first row or first column, just return true.
		if(cnt1 == global_aln_info.fore_hmm_array->n_length1 || cnt2 == global_aln_info.fore_hmm_array->n_length2)
		{
			return(true);
		}
		else
		{
			if(corr_aln_env[cnt1+1][cnt2] == ALN_ACK || corr_aln_env[cnt1][cnt2+1] == ALN_ACK || corr_aln_env[cnt1+1][cnt2+1] == ALN_ACK)
			{
				return(true);
			}
			else
			{
				return(false);
			}
		}
	}
	else
	{
		return(false);
	}
}

void dump_M_alignment_plane(int M)
{
	char file_name[100];
	sprintf(file_name, "%d_alignment_plane.txt", M);
	FILE* M_alignment_plane = fopen(file_name, "w");

	for(int cnt1 = 1; cnt1 < global_aln_info.fore_hmm_array->n_length1+1; cnt1++)
	{
		int min = cnt1 * global_aln_info.fore_hmm_array->n_length2 / global_aln_info.fore_hmm_array->n_length1 - M;
		int max = cnt1 * global_aln_info.fore_hmm_array->n_length2 / global_aln_info.fore_hmm_array->n_length1 + M;

		if(min <= 0)
		{
			min = 1;
		}

		if(max > global_aln_info.fore_hmm_array->n_length2)
		{
			max = global_aln_info.fore_hmm_array->n_length2;
		}

		for(int cnt2 = 1; cnt2 < global_aln_info.fore_hmm_array->n_length2+1; cnt2++)
		{
			if(cnt2 < min || cnt2 > max)
			{
				fprintf(M_alignment_plane, " 0");
			}
			else
			{
				fprintf(M_alignment_plane, " 1");
			}
		}

		fprintf(M_alignment_plane, "\n");
	}

	fclose(M_alignment_plane);
}

void calculate_M_from_aln_env(int* M)
{
	int last1, last2;
	for(int cnt1 = 1; cnt1 < global_aln_info.fore_hmm_array->n_length1+1; cnt1++)
	{
		for(int cnt2 = 1; cnt2 < global_aln_info.fore_hmm_array->n_length2+1; cnt2++)
		{
			// Must do a length scaling for determining M parameter.
			if(global_aln_info.corr_aln_env[cnt1][cnt2] == ALN_ACK && *M < abs((cnt1 * global_aln_info.fore_hmm_array->n_length2) / global_aln_info.fore_hmm_array->n_length1 - cnt2))
			{
				*M = abs((cnt1 * global_aln_info.fore_hmm_array->n_length2) / global_aln_info.fore_hmm_array->n_length1 - cnt2);					
			
				last1 = cnt1;
				last2 = cnt2;

				//printf("M = %d, abs = %d, cnt1 = %d, cnt2 = %d, l1 = %d, l2 = %d, \n", *M, abs((cnt1 * global_aln_info.fore_hmm_array->n_length2) / global_aln_info.fore_hmm_array->n_length1 - cnt2), cnt1, cnt2, global_aln_info.fore_hmm_array->n_length1, global_aln_info.fore_hmm_array->n_length2 );
			}
		}
	}
}

char** get_aln_env(char* seq1_name, char* seq2_name)
{
	load_output(seq1_name, seq2_name);

	// forward calculations.
	t_hmm_array hmm_fore_array;

	allocate_array(&hmm_fore_array);

	init_forward_array(&hmm_fore_array);

	calculate_forward_probs(&hmm_fore_array);

	// backward calculations.
	t_hmm_array hmm_back_array;

	allocate_array(&hmm_back_array);

	//init_backward_array2(&hmm_back_array);
	init_backward_array(&hmm_back_array);

	calculate_backward_probs(&hmm_back_array);

	calculate_aln_prob_array(&hmm_fore_array, &hmm_back_array);

	calculate_ins_prob_array(&hmm_fore_array, &hmm_back_array);

	int M;
	double P_thresh = -10.0;
	get_alignment_envelopes(exp(P_thresh), NULL, &M);

	return(global_aln_info.corr_aln_env);
}

// Checks if k is in the range of i using constraint k.
bool check_M_range(int M, int i, int k, int N1, int N2)
{
	int lowlimit_i = i * N2 / N1 - M;
	int highlimit_i = i * N2 / N1 + M;
	if(k <= highlimit_i && k >= lowlimit_i)
	{
		return(true);
	}
	else
	{
		//printf("%d is not in range of %d, (%d, %d) (l1, l2) = (%d, %d)\n", k, i, lowlimit_i, highlimit_i, N1, N2);
		return(false);
	}
}

int get_n_M_const(int M, int N1, int N2)
{
	int size = 0;

	for(int cnt = 1; cnt < N1 + 1; cnt++)
	{
		int lowlimit_i = cnt * N2 / N1 - M;
		int highlimit_i = cnt * N2 / N1 + M;

		if(lowlimit_i <= 0)
		{
			lowlimit_i = 1;
		}

		if(highlimit_i > N2)
		{
			highlimit_i = N2;
		}		
		
		size += (highlimit_i - lowlimit_i + 1);
	}

	return(size);
}

// Following function checks for the connectivity of a whole alignment envelope,
// Basically what it needs is it has make sure that the alignment envelope can give
// at least one valid alignment between two sequences.
bool check_connectivity(char** aln_env)
{
	// One algorithm is iteratively going over rows and columns of alignment envelope until 
	// there is no change that is directed by: Mark all the cells as connected that are on same row or column of
	// an aligned cell, do this over rows and columns iteratively until there is no change, if the
	// connectivity is reached at last column or row, this alignment enveope is connected.

	// Allocate connectivity matrix same size of aln. env.
	bool** conn_matrix = (bool**)malloc(sizeof(bool*) * (global_aln_info.fore_hmm_array->n_length1 + 3));

	for(int cnt1 = 0; cnt1 < global_aln_info.fore_hmm_array->n_length1 + 2; cnt1++)
	{
		conn_matrix[cnt1] = (bool*)malloc(sizeof(bool) * (global_aln_info.fore_hmm_array->n_length2 + 3));

		for(int cnt2 = 0; cnt2 < global_aln_info.fore_hmm_array->n_length2 + 2; cnt2++)
		{
			conn_matrix[cnt1][cnt2] = false; // Set all entries as false at first.
		}
	}

	bool start_connection = false;
	bool end_connection = false;

	bool changed = true;
	while(changed)
	{
		changed = false;
		// Go over rows and mark the cells at current row which have connectivity to aligned cells in aln. env as connected.
		for(int cnt1 = 1; cnt1 < global_aln_info.fore_hmm_array->n_length1 + 1; cnt1++)
		{
			// Going over rows.
			for(int cnt2 = 1; cnt2 < global_aln_info.fore_hmm_array->n_length2 + 1; cnt2++)
			{
				if(aln_env[cnt1][cnt2] == ALN_ACK && !conn_matrix[cnt1][cnt2])
				{
					changed = true;
					conn_matrix[cnt1][cnt2] = true;

					// If iterations reached the border, connectivity is established.
					if(cnt1 == 1 || cnt2 == 1)
					{
						start_connection = true;
					}

					// If iterations reached the border, connectivity is established.
					if(cnt1 == global_aln_info.fore_hmm_array->n_length1 || cnt2 == global_aln_info.fore_hmm_array->n_length2)
					{
						end_connection = true;
					}
				}
			}
		}

		// Go over columns and do the same thing for columns.
		for(int cnt2 = 1; cnt2 < global_aln_info.fore_hmm_array->n_length2 + 1; cnt2++)
		{
			// Going over columns.
			for(int cnt1 = 1; cnt1 < global_aln_info.fore_hmm_array->n_length1 + 1; cnt1++)
			{
				if(aln_env[cnt1][cnt2] == ALN_ACK && !conn_matrix[cnt1][cnt2])
				{
					changed = true;
					conn_matrix[cnt1][cnt2] = true;

					// If iterations reached the border, start connectivity is established.
					if(cnt1 == 1 || cnt2 == 1)
					{
						start_connection = true;
					}

					// If iterations reached the border, end connectivity is established.
					if(cnt1 == global_aln_info.fore_hmm_array->n_length1 || cnt2 == global_aln_info.fore_hmm_array->n_length2)
					{
						end_connection = false;
					}
				}
			}
		}

		if(start_connection && end_connection)
		{
			return(true);
		}

		// Meanwhile check if there is any change made, set changed accordingly.

		// If any cells are changed, then a new iteration begins over rows and columns.
	}

	return(false);
}



