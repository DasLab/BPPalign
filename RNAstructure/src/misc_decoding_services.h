#ifndef _MISC_DECODING_SERVICES_
#define _MISC_DECODING_SERVICES_

// Misc. decoding services.
void decode_posterior();
void decode_ML(short**);
void get_windowed_ins_probs(int win_length);
void get_windowed_ins_probs2(int win_length);
////////////////////////////////////////////////////////////////////////////

void free_ML_arrays(double*** score_array, char*** ML_path_array);

#endif
