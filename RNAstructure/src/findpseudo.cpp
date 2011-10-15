

#include "findpseudo.h"


//find pseudoknots, break them, and return the number of pseudoknots and their location
void findpseudo(structure *ct, int structnum, int *npseudo, int *pseudo) {
	register int i,j,k,L;
	int cons,incons;


   *npseudo = 0;
  for (i=1;i<=ct->numofbases;i++) {
   	if (ct->basepr[structnum][i]>i) {//basepair found
		for (j=ct->numofbases;j>ct->basepr[structnum][i];j--) {

			if ((ct->basepr[structnum][j]<ct->basepr[structnum][i])&&(ct->basepr[structnum][j]>i)) {
				//a pseudoknot has been found
				cons=1;
				incons=1;
				for (k=i+1;k<j;k++) {
					if (ct->basepr[structnum][k]>k) {
						
						if(ct->basepr[structnum][k]>ct->basepr[structnum][i]) incons++;
						else cons++;

					

					}


				}
				if (cons>incons) {
					L = ct->basepr[structnum][i];
					for (k=i;k<L;k++) {

						if (ct->basepr[structnum][k]>ct->basepr[structnum][i]) {
							(*npseudo)++;
							if ((*npseudo)<maxpseudo) pseudo[*npseudo]=k;
							(*npseudo)++;
							if ((*npseudo)<maxpseudo) pseudo[*npseudo]=ct->basepr[structnum][k];
							ct->basepr[structnum][ct->basepr[structnum][k]]=0;
							ct->basepr[structnum][k]=0;
				

						}

					}


				}
				else {
					L = ct->basepr[structnum][i];
					for (k=i;k<L;k++) {

						if (ct->basepr[structnum][k]<=ct->basepr[structnum][i]&&ct->basepr[structnum][k]>0) {
							(*npseudo)++;
							if ((*npseudo)<maxpseudo) pseudo[*npseudo]=k;
							(*npseudo)++;
							if ((*npseudo)<maxpseudo) pseudo[*npseudo]=ct->basepr[structnum][k];
							ct->basepr[structnum][ct->basepr[structnum][k]]=0;
							ct->basepr[structnum][k]=0;
				

						}

					}


				}


			}
		}
   }
}

}

//return true if the structure #StructureNumber has a pseudoknot or false otherwise
bool HasPseudo(structure *ct, int StructureNumber) {
	int i,j;

	//use the standard definition of a pseudoknot to find pseudoknots
	for (i=1;i<=ct->numofbases;i++) {
		if (ct->basepr[StructureNumber][i]>i) {
			//found a pair
			for (j=i+1;j<ct->basepr[StructureNumber][i];j++) {
				if (ct->basepr[StructureNumber][j]>0&&ct->basepr[StructureNumber][j]<i) return true;
				else if (ct->basepr[StructureNumber][j]>ct->basepr[StructureNumber][i]) return true;

			}
		}
	}


	return false;//no pseudoknot found

}