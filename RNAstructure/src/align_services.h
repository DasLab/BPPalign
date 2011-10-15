#ifndef _DAVE_SERVICES_
#define _DAVE_SERVICES_

#include "structure.h"

int calculate_aln_probs_env(structure *ct1, structure *ct2, double** align_probs, bool** aln_env, double p_thresh);

int calculate_coinc_probs_env(structure *ct1, structure *ct2, double** coinc_probs, bool** aln_env, char* data_path, short** forcealign);

#endif

