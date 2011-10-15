#include "align_services.h"
#include "misc_decoding_services.h"
#include <stdio.h>
#include "hmm_arrays.h"
#include "xlog_math.h"
#include <string.h>
#include "process_output.h"
#include "structure.h"
#include "alignment_hmm_model.h"
#include "aln_env_utils.h"
#include <stdlib.h>

extern t_aln_info global_aln_info;

extern char* ml_aln_line1;
extern char* ml_aln_line2;

int calculate_aln_probs_env(structure *ct1, structure *ct2, double** align_probs, bool** aln_env, double p_thresh)
{
    // Load these sequence strings to hmm output manipulator.
    load_raw_output(ct1->nucs, ct2->nucs);

    // Calculate pair alignment probabilities.
    // forward calculations.
    t_hmm_array hmm_fore_array;

    allocate_array(&hmm_fore_array);

    init_forward_array(&hmm_fore_array, NULL);

    calculate_forward_probs(&hmm_fore_array, NULL);

    // backward calculations.
    t_hmm_array hmm_back_array;

    allocate_array(&hmm_back_array);

    init_backward_array2(&hmm_back_array, NULL);

    calculate_backward_probs2(&hmm_back_array, NULL);
    calculate_aln_prob_array(&hmm_fore_array, &hmm_back_array);

	
	// Copy alignment probabilities to argument.
	for(int cnt1 = 1; cnt1 <= hmm_back_array.n_length1; cnt1++)
	{
		for(int cnt2 = 1; cnt2 <= hmm_fore_array.n_length2; cnt2++)
		{
			double fore_aln_prob = hmm_fore_array.probs[cnt1][cnt2][STATE_ALN];
			double back_aln_prob = hmm_back_array.probs[cnt1][cnt2][STATE_ALN];

			// If align_probs is not NULL, calculate align probs.
			if(align_probs != NULL)
			{
				align_probs[cnt1][cnt2] = xlog_div(xlog_mul(fore_aln_prob, back_aln_prob), global_aln_info.op_prob);
			}

			if(aln_env != NULL)
			{
				// Determine alignment envelope.
				if(xlog_div(xlog_mul(fore_aln_prob, back_aln_prob), global_aln_info.op_prob) >= p_thresh)
				{
					aln_env[cnt1][cnt2] = true;
				}
				else
				{
					aln_env[cnt1][cnt2] = false;
				}
			}
		}
	}

	return(1); // return 1 on success.
}

int calculate_coinc_probs_env(structure *ct1, structure *ct2, double** coinc_probs, bool** aln_env, char* data_path, short** forcealign)
{
	//printf("In calculate_coinc_probs_env\n");

	// Init all knowing structure.
	init_global_aln_info();

	// Initialize family of parameters.
	char hmm_par_file_path[1000];

	// If argument to data path is NULL, ask for it from environment variables.
	if(data_path == NULL)
	{
		// Care must be taken since data_path is on stack.
		data_path = getenv("DATAPATH");
	}

	// Do environment variables contain data path?
        if(data_path != NULL) 
	{
		sprintf(hmm_par_file_path, "%s/fam_hmm_pars.dat", data_path);
	}
        else // If environment variables do not contain a data path, assume that parameter file is in current working directory.
	{
		sprintf(hmm_par_file_path, "fam_hmm_pars.dat");
	}

	FILE* test_par_file = fopen(hmm_par_file_path, "r");
	if(test_par_file == NULL) // Check if assumed file path exists and contains hmm parameter file.
	{
		printf("Could not locate HMM parameter file %s at %s(%d)\n", hmm_par_file_path, __FILE__, __LINE__);
		exit(0);
	}
	else
	{
		fclose(test_par_file);
	}

	//load_fam_hmm_parameters("fam_hmm_pars.txt");
	load_fam_hmm_parameters(hmm_par_file_path);

	// Load output.
	// NOTE THAT ct1->nucs IS A STRING WHICH STARTS FROM INDEX 1. THIS MUST BE COMPENSATED IN OUTPUT LOADING!!!
    	load_raw_output(ct1->nucs, ct2->nucs);

	//printf("ML computation!\n");

	// Do ML decoding on alignment.
	decode_ML(forcealign);

	//printf("ML computation ok\n");

	//exit(0);

	// Calculate similarity of alignment from ML alignment.
	double ml_similarity = get_aln_similarity(ml_aln_line1, ml_aln_line2);

	// Load correct parameters from family of parameters for forward-backward calculation.
	choose_hmm_parameters(ml_similarity);

	// forward calculations.
	t_hmm_array hmm_fore_array;

	allocate_array(&hmm_fore_array);

	init_forward_array(&hmm_fore_array, forcealign);

	calculate_forward_probs(&hmm_fore_array, forcealign);

	// backward calculations.
	t_hmm_array hmm_back_array;

	allocate_array(&hmm_back_array);

	init_backward_array2(&hmm_back_array, forcealign);

	// Calculate backward probabilities.
	calculate_backward_probs2(&hmm_back_array, forcealign);

	// Calculate alignment probabilities.
	calculate_aln_prob_array(&hmm_fore_array, &hmm_back_array);

	// Get probability of thresh.
	double P_thresh = get_fam_threshold(ml_similarity);
	//printf("sim = %.3f -> P_thresh = %.3f\n", ml_similarity, P_thresh);


	//FILE* f_aln_env = fopen("aln_env.txt", "w"); 
	// Copy alignment probabilities to argument.
	for(int cnt1 = 1; cnt1 <= hmm_back_array.n_length1; cnt1++)
	{
		for(int cnt2 = 1; cnt2 <= hmm_fore_array.n_length2; cnt2++)
		{
			double fore_aln_prob = hmm_fore_array.probs[cnt1][cnt2][STATE_ALN];
			double back_aln_prob = hmm_back_array.probs[cnt1][cnt2][STATE_ALN];

			double ins1_prob = xlog_mul(global_aln_info.fore_hmm_array->probs[cnt1][cnt2][STATE_INS1], global_aln_info.back_hmm_array->probs[cnt1][cnt2][STATE_INS1]);
			double ins2_prob = xlog_mul(global_aln_info.fore_hmm_array->probs[cnt1][cnt2][STATE_INS2], global_aln_info.back_hmm_array->probs[cnt1][cnt2][STATE_INS2]);
			double three_plane_sum = xlog_sum(global_aln_info.aln_probs[cnt1][cnt2], xlog_sum(ins1_prob, ins2_prob));

			// If align_probs is not NULL, calculate align probs.
			if(coinc_probs != NULL)
			{
				coinc_probs[cnt1][cnt2] = xlog_div(three_plane_sum, global_aln_info.op_prob);
			}

			if(aln_env != NULL)
			{
				// Determine alignment envelope.
				if(xlog_div(three_plane_sum, global_aln_info.op_prob) >= P_thresh)
				{
					aln_env[cnt1][cnt2] = true;
				//	fprintf(f_aln_env, "1 ");
				}
				else
				{
					aln_env[cnt1][cnt2] = false;
					//fprintf(f_aln_env, "0 ");
				}
			}
		}
		//fprintf(f_aln_env, "\n");
	}
	//fclose(f_aln_env);

	// Free!!!
	free_global_aln_info();
	free_array(&hmm_fore_array);
	free_array(&hmm_back_array);

	return(1); // return 1 on success.
}
