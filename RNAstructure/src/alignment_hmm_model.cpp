#include "alignment_hmm_model.h"
#include <stdio.h>
#include <stdlib.h>

// Following are emission and transition probabilities that current ML/forward-backward calculations are using.
// Those should be loaded correctly before doing these calculations.
double emit_probs[N_OUTPUTS][N_STATES];
double trans_probs[N_STATES][N_STATES];

double fam_hmm_pars[N_BINZ * (N_STATES + N_OUTPUTS) * N_STATES];
double fam_thresholds[N_BINZ];

void set_hmm_parameters(double new_emit_probs[N_OUTPUTS][N_STATES], double new_trans_probs[N_STATES][N_STATES])
{
	// Copy transition matrix.
	for(int cnt1 = 0; cnt1 < N_STATES; cnt1++)
	{
		for(int cnt2 = 0; cnt2 < N_STATES; cnt2++)
		{
			trans_probs[cnt1][cnt2] = new_trans_probs[cnt1][cnt2];
		}
	}

	// Copy emission probabilities.
	for(int cnt1 = 0; cnt1 < N_OUTPUTS; cnt1++)
	{
		for(int cnt2 = 0; cnt2 < N_STATES; cnt2++)
		{
			emit_probs[cnt1][cnt2] = new_emit_probs[cnt1][cnt2];
		}
	}
}

void load_fam_hmm_parameters(char* fam_par_file_name)
{
	FILE* fam_par_file = fopen(fam_par_file_name, "r");

	// Load all parameters from file.
	for(int cnt = 0; cnt < N_BINZ * (N_STATES + N_OUTPUTS) * N_STATES; cnt++)
	{
		//double cur_val = 0.0;
		char cur_num_str[20];
		fscanf(fam_par_file, "%s", cur_num_str);
		fam_hmm_pars[cnt] = atof(cur_num_str);
		//printf("%d: %.3f\n", cnt, fam_hmm_pars[cnt]);
	}

	//printf("\n\n\nReading thresholds!!!\n");
	// Read thresholds.
	for(int cnt = 0; cnt < N_BINZ; cnt++)
	{
		char cur_num_str[20];
		fscanf(fam_par_file, "%s", cur_num_str);
		fam_thresholds[cnt] = atof(cur_num_str);	
		//printf("%d: %f\n", cnt, fam_thresholds[cnt]);
	}

	fclose(fam_par_file);
}

// Choosese and loads hmm parameters from family parameters.
void choose_hmm_parameters(double similarity)
{
	//int fam_par_set_index = (int)(similarity * (double)N_BINZ);
	int fam_par_set_index = get_bin_index(similarity, N_BINZ);

	// Load emission probabilities.
	// Each parameter set is (N_STATES + N_OUTPUTS) * N_STATES doubles long.
	int start_linear_index = (N_STATES + N_OUTPUTS) * N_STATES * get_bin_index(similarity, N_BINZ);
	double* par_ptr = fam_hmm_pars + start_linear_index;

	for(int cnt1 = 0; cnt1 < N_OUTPUTS; cnt1++)
	{
		for(int cnt2 = 0; cnt2 < N_STATES; cnt2++)
		{
			emit_probs[cnt1][cnt2] = *(par_ptr + cnt1 * N_STATES + cnt2);
		}
	}

	start_linear_index = (N_STATES + N_OUTPUTS) * N_STATES * fam_par_set_index + N_STATES * N_OUTPUTS;
	par_ptr = fam_hmm_pars + start_linear_index;

	// Load trans probabilities.
	for(int cnt1 = 0; cnt1 < N_STATES; cnt1++)
	{
		for(int cnt2 = 0; cnt2 < N_STATES; cnt2++)
		{
			trans_probs[cnt1][cnt2] = *(par_ptr + cnt1 * N_STATES + cnt2);
		}
	}		

	// Dump emission probabilities.
/*
	for(int cnt1 = 0; cnt1 < N_OUTPUTS; cnt1++)
	{
		for(int cnt2 = 0; cnt2 < N_STATES; cnt2++)
		{
			printf("%.3f ", emit_probs[cnt1][cnt2]);
		}

		printf("\n");
	}

	// Dump transition probabilities.
	printf("\n");
	for(int cnt1 = 0; cnt1 < N_STATES; cnt1++)
	{
		for(int cnt2 = 0; cnt2 < N_STATES; cnt2++)
		{
			printf("%.3f ", trans_probs[cnt1][cnt2]);
		}

		printf("\n");
	}
*/
}

// Get index of bin of parameters for a sequence alignment.
int get_bin_index(double similarity, int n_bins)
{
	if(similarity == 1.0)
	{
		return(n_bins - 1);
	}
	else
	{
		return((int)(n_bins * similarity));
	}
}

double get_fam_threshold(double similarity)
{
	int bin_index = get_bin_index(similarity, N_BINZ);

	return(fam_thresholds[bin_index]);
}
