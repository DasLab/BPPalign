#include "aln_env_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double get_aln_similarity(char* aln_line1, char* aln_line2)
{
	//printf("%s\n%s\n", aln_line1, aln_line2);

	if(strlen(aln_line1) != strlen(aln_line2))
	{
		printf("alignment lines are not of same length, exiting at %s(%d)\n", __FILE__, __LINE__);
		exit(0);
	}

	int n_match_pos = 0;

	for(int cnt = 0; cnt < strlen(aln_line1); cnt++)
	{
		if(aln_line1[cnt] != '.' && aln_line1[cnt] == aln_line2[cnt])
		{
			n_match_pos++;
		}
	}

	int n_total_pos = 0;

	for(int cnt = 0; cnt < strlen(aln_line1); cnt++)
	{
		if(aln_line1[cnt] == '.' && aln_line2[cnt] == '.')
		{
		}
		else
		{
			n_total_pos++;
		}
	}

	return( (double)n_match_pos / n_total_pos );
}

// Returns fraction of aligned positions by total alignment positions.
double get_aln_percentage(char* aln_line1, char* aln_line2)
{
	//printf("%s\n%s\n", aln_line1, aln_line2);

	if(strlen(aln_line1) != strlen(aln_line2))
	{
		printf("alignment lines are not of same length, exiting at %s(%d)\n", __FILE__, __LINE__);
		exit(0);
	}

	int n_match_pos = 0;

	for(int cnt = 0; cnt < strlen(aln_line1); cnt++)
	{
		if(aln_line1[cnt] != '.' && aln_line2[cnt] != '.')
		{
			n_match_pos++;
		}
	}

	int n_total_pos = 0;

	for(int cnt = 0; cnt < strlen(aln_line1); cnt++)
	{
		if(aln_line1[cnt] == '.' && aln_line2[cnt] == '.')
		{
		}
		else
		{
			n_total_pos++;
		}
	}

	return( (double)n_match_pos / n_total_pos );
}
