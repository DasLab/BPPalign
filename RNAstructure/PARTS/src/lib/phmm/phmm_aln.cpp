#include "phmm_aln.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../structure/structure.h"
#include "../utils/xmath/log/xlog_math.h"
#include "phmm.h"
#include "phmm_array.h"
#include "p_alignment.h"
#include <ctype.h>

bool _DUMP_PHMM_ALN_MESSAGES_ = false;

t_phmm_aln::t_phmm_aln(char* seq1_fp, char* seq2_fp, int _phmm_band_constraint_size)
{
	// Must use a copy constructor here!
	this->phmm = NULL;
	this->phmm_band_constraint_size = _phmm_band_constraint_size;

	this->seq1 = new t_structure(seq1_fp);
	this->seq2 = new t_structure(seq2_fp);
	this->check_set_seqs();
}

t_phmm_aln::t_phmm_aln(t_structure* _seq1, t_structure* _seq2, int _phmm_band_constraint_size)
{
	// Must use a copy constructor here!
	this->phmm = NULL;
	this->phmm_band_constraint_size = _phmm_band_constraint_size;

	this->seq1 = new t_structure(_seq1);
	this->seq2 = new t_structure(_seq2);
	this->check_set_seqs();
}

t_phmm_aln::t_phmm_aln(t_structure* _seq1, t_structure* _seq2)
{
	// Must use a copy constructor here!
	this->phmm = NULL;
	this->phmm_band_constraint_size = DEFAULT_PHMM_BAND_CONSTRAINT;

	this->seq1 = new t_structure(_seq1);
	this->seq2 = new t_structure(_seq2);
	this->check_set_seqs();
}


t_phmm_aln::t_phmm_aln(char* seq1_fp, char* seq2_fp)
{
	// Must use a copy constructor here!
	this->phmm = NULL;
	this->phmm_band_constraint_size = 0x1fffff;

	this->seq1 = new t_structure(seq1_fp);
	this->seq2 = new t_structure(seq2_fp);

	this->seq1 = new t_structure(seq1_fp);
	this->seq2 = new t_structure(seq2_fp);
	this->check_set_seqs();
}

void t_phmm_aln::check_set_seqs()
{
	// Check and verify the sequences, this includes replacing unknown nucleotides by random nucleotides.
	srand((int)time(0));
	
	// Take care of unknown nucleotides by randomly selecting them.
	for(int cnt = 1; cnt <= seq1->numofbases; cnt++)
	{
		if(toupper(seq1->nucs[cnt]) != 'A' &&
		toupper(seq1->nucs[cnt]) != 'C' &&
		toupper(seq1->nucs[cnt]) != 'G' &&
		toupper(seq1->nucs[cnt]) != 'T' &&
		toupper(seq1->nucs[cnt]) != 'U')
		{
			seq1->nucs[cnt] = generate_random_nuc();
		}

		seq1->numseq[cnt] = this->nuc2num(seq1->numseq[cnt]);
	}

	for(int cnt = 1; cnt <= seq2->numofbases; cnt++)
	{
		if(toupper(seq2->nucs[cnt]) != 'A' &&
		toupper(seq2->nucs[cnt]) != 'C' &&
		toupper(seq2->nucs[cnt]) != 'G' &&
		toupper(seq2->nucs[cnt]) != 'T' &&
		toupper(seq2->nucs[cnt]) != 'U')
		{
			seq2->nucs[cnt] = generate_random_nuc();
		}

		seq2->numseq[cnt] = this->nuc2num(seq2->numseq[cnt]);
	}
}

int* t_phmm_aln::get_seq2_aln_const(int* seq1_aln_const)
{
	int* seq2_aln_const = NULL;
	if(seq1_aln_const != NULL)
	{
		seq2_aln_const = (int*)malloc(sizeof(int) * (this->l2() + 2));
		for(int i2 = 0; i2 <= this->l2(); i2++)
		{
			// Init the constraint for this position as 0, then check seq1 alignment constraints.
			seq2_aln_const[i2] = 0;

			for(int i1 = 0; i1 <= this->l1(); i1++)
			{
				if(seq1_aln_const[i1] != 0 &&
					seq1_aln_const[i1] == i2)
				{
					seq2_aln_const[i2] = i1;
				}
			} // i1 loop
		} // i2 loop
	}

	return(seq2_aln_const);
}

void t_phmm_aln::get_aln_permissions(int* seq1_aln_const, 
									int* seq2_aln_const, 
									bool& forbid_STATE_ALN, 
									bool& forbid_STATE_INS1, 
									bool& forbid_STATE_INS2, 
									int i, 
									int k)
{
	if(seq1_aln_const == NULL)
	{
		// No constraints!
		forbid_STATE_ALN = false;
		forbid_STATE_INS1 = false;
		forbid_STATE_INS2 = false;			
	}
	else
	{ 
		// Is there a constrained for i?
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
		else // There is no constraint.
		{
			forbid_STATE_ALN = false;
			forbid_STATE_INS1 = false;
			forbid_STATE_INS2 = false;				
		}
	} // NULL check for alignment constraints.
} // get_aln_permissions function.

// Following function generates random nucleotide in place of an unknown nucleotide.
char t_phmm_aln::generate_random_nuc()
{
	int rand_nuc = rand() % 4;

	switch(rand_nuc)
	{
	case 0:
		return('A');
		break;
	case 1:
		return('C');
		break;
	case 2:
		return('G');
		break;
	case 3:
		return('U');
		break;
	default:
		printf("Invalid random nuc!!!\n");
		exit(0);
		break;
	};
}

t_phmm_aln::~t_phmm_aln()
{
	delete(this->seq1);
	delete(this->seq2);

	//if(this->phmm != NULL)
	//{
	//	delete(this->phmm);
	//}
}

int t_phmm_aln::l1()
{
	return(this->seq1->numofbases);
}

int t_phmm_aln::l2()
{
	return(this->seq2->numofbases);
}

int t_phmm_aln::nuc2num(char nuc)
{
	if(nuc == 'A' || nuc == 'a')
	{
		return(0);
	}
	else if(nuc == 'C' || nuc == 'c')
	{
		return(1);
	}
	else if(nuc == 'G' || nuc == 'g')
	{
		return(2);
	}
	else if(nuc == 'U' || nuc == 'u' || nuc == 'T' || nuc == 't')
	{
		return(3);
	}
	else
	{
		return(4);
	}
}

double t_phmm_aln::get_trans_emit_prob(int prev_state, int current_state, int i, int k)
{
	double current_trans_prob = this->phmm->get_trans_prob(prev_state, current_state);

	int i_sym;
	int k_sym;

	// Fix symbols to gaps in case of of insertions.
	if(current_state == STATE_INS1 || k == 0)
	{
		// Gap is coded into value 4 in the emission table.
		k_sym = 4; 
	}
	else
	{
		k_sym = this->nuc2num(seq2->nucs[k]);
	}

	if(current_state == STATE_INS2 || i == 0)
	{
		// Gap is coded into value 4 in the emission table.
		i_sym = 4;
	}
	else
	{
		i_sym = this->nuc2num(seq1->nucs[i]);
	}

	// Compute the symbol index into emission table using the coded nucleotide values:
	// A->0, C->1, G->2, U->3, T->3, .->4
	// This defines a counting system in base of 5. (25 values.)
	// There are also emission of start and end symbols. These correspond to 25th and 26th indices in the emission probability table.
	int sym_index = i_sym * 5 + k_sym;

	// Check for exceptional cases of start and end symbols.
	// The indices correspond to the start symbol?
	if(i == 0 && k == 0)
	{
		sym_index = 25;
	}

	// The indices correspond to the end symbol?
	if(i == (l1() + 1) && k == (l2() + 1))
	{
		sym_index = 26;
	}

	double current_emission_prob = this->phmm->get_emit_prob(sym_index, current_state);

	return(xlog_mul(current_emission_prob, current_trans_prob));
}
