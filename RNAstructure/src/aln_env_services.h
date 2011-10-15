#ifndef _ALN_ENV_SERVICES_
#define _ALN_ENV_SERVICES_

// Alignment envelope services, includes probabilistic and M heuristic and extendible in future.
void get_alignment_envelopes(double threshold_prob, bool** aln_env, int* M);
char** get_aln_env(char* seq1_name, char* seq2_name);
bool check_connection(int cnt1, int cnt2, bool** aln_env, char** corr_aln_env);
bool check_forward_connection(int cnt1, int cnt2, char** pruned_aln_env, char** corr_aln_env);
void prune_aln_env(bool** aln_env);
void copy_aln_env(bool** aln_env);
void calculate_M_from_aln_env(int* M);
void dump_M_alignment_plane(int M);
//////////////////////////////////////////////////////////////////////////////////

int get_n_M_const(int M, int N1, int N2);
bool check_M_range(int M, int i, int k, int N1, int N2);

#endif
