#include "hmm_probs.h"
#include "alignment_hmm_model.h"
#include "xlog_math.h"
#include <stdio.h>

extern double emit_probs[N_OUTPUTS][N_STATES];
extern double trans_probs[N_STATES][N_STATES];

double get_trans_prob(int from, int to)
{
	return(trans_probs[from][to]);
}

double get_log_trans_prob(int from, int to)
{
	if(trans_probs[from][to] < 0)
	{
		printf("trans_prob < 0: %d -> %d\n", from, to);
	}

	return(xlog(trans_probs[from][to]));
}

double get_log_emit_prob(int state, int symbol)
{
	if(emit_probs[symbol][state] < 0)
	{
		printf("emit_prob < 0: %d emitting %d\n", state, symbol);
	}

	return(xlog(emit_probs[symbol][state]));
}

