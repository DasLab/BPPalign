/*
 * This is the command-line interface to NAPSS, 
 * written by James Hart
 *
 * Modified from the command-line interface to Dynalign, 
 * written by David Mathews; Copyright 2002, 2003, 2004, 2005, 2006
 *
 *----------------------------------------------------------------
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cstring>
#include <cstdlib>

//#include "stdafx.h"
#include "../src/configfile.h"
#include "../src/defines.h"
#include "../src/platform.h"
#include "../src/rna_library.h"
#include "../src/structure.h"
#include "../src/algorithm.h"

#define BETA_1 96
#define BETA_1PM 150
#define BETA_2 1
#define BETA_3 1

// Flags for text output
//#undef VERBOSE_MODE
#define VERBOSE_MODE
#undef VERY_VERBOSE_MODE
//#define VERY_VERBOSE_MODE

// Flags for debug features
#undef DEBUG_MODE
//#define DEBUG_MODE
//#undef ALREADY_USED_CHECK
#define ALREADY_USED_CHECK


using namespace std;

void getdat(char *loop, char *stackf, char *tstackh, char *tstacki,
	char *tloop, char *miscloop, char *danglef, char *int22,
	char *int21,char *coax, char *tstackcoax,
    char *coaxstack, char *tstack, char *tstackm, char *triloop,
    char *int11, char *hexaloop, char *tstacki23, char *tstacki1n,
	char *datapath, bool isRNA);

struct conMatch {
	short* coords;
};

void firstDotMatch(short**,bool*,bool*,short,short,short**,short**,vector<conMatch>*,short*,short*,short**,short**,int,short);
			// Loops through x,y of entire dotplot searching for matches to the first basepair of each constraint

void recursiveMatch(short**,bool*,short,short,short**,short**,vector<conMatch>*,short*,short*,short**,short**,int,short,bool*);
			// Takes over for firstDotMatch to find the matches to the remaining basepairs in a constraint

void pseudodptrim(structure*,int*, int*); // Searches dotplot for complicated pseudoknot folds to exclude

void efn2mod(datatable *data,structure *ct, int structnum, bool simplemb, structure *ctbroken);
			// Pseudoknot-capable energy calculation

int pseudoenergy(datatable *data, structure *ct, int count, bool *nopseudos, structure *ctbroken, int *init);
			// Searches ct for pseudoknots, returns energy for breaking apart the pseudoknot (including Beta penalties)

int ctout2(structure *ct, int cutoff, const char *ctoutfile); 
			// Modified version of ctout that can truncate unstable structures; returns number of output structures

void pairout(structure *ct, int cutoff, const char* pairsoutfile); // Outputs PseudoViewer3 structural data

// The main entry point for NAPSS using a text interface:
int main(int argc, char* argv[]) {
	string inseq;
	string indotplot;
	string inconstraints;
	string outct;
	string outpairs;

	char loop[maxfil],stackf[maxfil],tstackh[maxfil],tstacki[maxfil],
		tloop[maxfil],miscloop[maxfil],danglef[maxfil],int22[maxfil],
		int21[maxfil],coax[maxfil],tstackcoax[maxfil],
		coaxstack[maxfil],tstack[maxfil],tstackm[maxfil],triloop[maxfil],int11[maxfil],hexaloop[maxfil],
		tstacki23[maxfil], tstacki1n[maxfil],datapath[maxfil],*pointer;

	int maxtracebacks,percent,windowsize,cutoff;
	   
	structure ct;
	datatable data;

	int i,j,k,iter,iter2;
	string a;


	// Check if we have a command-line parameter; if not, fall back to
	// interactive mode.
	if( argc == 2 ) {
		ConfigFile config(argv[1]);

		// Check for mandatory config file settings.
		bool valid_config =
			config.contains("inseq") &&
			config.contains("indotplot") &&
			config.contains("inconstraints") &&
			config.contains("outct");
	
		if (valid_config) {
			// Read all mandatory config file settings.
			inseq         = config.getOption<string>("inseq");
			indotplot     = config.getOption<string>("indotplot");
			inconstraints = config.getOption<string>("inconstraints");
			outct         = config.getOption<string>("outct");
		}

		// Read optional config file settings
		if (valid_config) {
			if (config.contains("outpairs")) outpairs = config.getOption<string>("outpairs");
			else outpairs = ""; // Indicates that pairs output file is not desired
			if (config.contains("maxtracebacks")) maxtracebacks = config.getOption<int>("maxtracebacks");
			else maxtracebacks = 100;
			if (config.contains("percent")) percent = config.getOption<int>("percent");
			else percent = 25;
			if (config.contains("windowsize")) windowsize = config.getOption<int>("windowsize");
			else windowsize=0;
			if (config.contains("cutoff")) cutoff = config.getOption<int>("cutoff");
			else cutoff = 0; // Indicates not to truncate output ct file
		}
		
		if (!valid_config) {
			cerr << "ERROR: At least one parameter could not be read from the "
				<< "configuration file. Aborting." << endl;
			cin>>i;
			exit(1);
		}
	} 

	else {
		// Interactive mode:
		cout << "Usage: NAPSS [config file]\nNo config file specified, using interactive mode:\n\n";
		
		cout << "Enter the name of the sequence file: ";
		cin >> inseq;
		
		cout << "Enter the name of the dotplot file: ";
		cin >> indotplot;

		cout << "Enter the name of the NMR constraints file: ";
		cin >> inconstraints;

		cout << "Enter the name of the output ct file: ";
		cin >> outct;

		cout << "Enter file name for the Positions Paired style output\n(enter 0 if not desired): ";
		cin >> outpairs;
		if (outpairs == "0") outpairs = "";

		cout << "Enter the maximum number of structures for each\nNMR constraint match combination: ";
		cin >> maxtracebacks;
		
		cout << "Enter the maximum percent energy difference to consider in the dotplot\n(where entering 25 is 25%): ";
		cin >> percent;

		cout << "Enter the base pair window size <recommend 0>: ";
		cin >> windowsize;

		cout << "Enter the maximum percent energy difference for valid output structures\n(enter 0 to output all structures): ";
		cin >> cutoff;
	}

	// Get the location of the data files
	pointer = getenv("DATAPATH");
	if (pointer!=NULL) {
		strcpy(datapath,pointer);
		strcat(datapath,"/");
	}
	else strcpy(datapath,"");

	// Get the location of the data files
	//	strcpy(datapath,"./data/");

	// Open the thermodynamic data tables
	getdat (loop, stackf, tstackh, tstacki,tloop, miscloop, danglef, int22,
		int21,coax, tstackcoax,coaxstack, tstack, tstackm, triloop,
		int11, hexaloop, tstacki23, tstacki1n, datapath, true);//the true indicates RNA parameters
	if (opendat (loop, stackf, tstackh, tstacki,tloop, miscloop, danglef, int22,
		int21,coax, tstackcoax,coaxstack, tstack, tstackm, triloop,
		int11,hexaloop,tstacki23, tstacki1n, &data)==0) {
			cerr << "A data file was lost";
			cin >> i;
	}

	// Open the sequence
	if(!openseq(&ct,inseq.c_str())){cerr<<"ERROR: Could not open sequence file "<<inseq<<"\n";cin>>i;return -1; }

#if defined(VERY_VERBOSE_MODE)
	// Display basic info
	cout << "\n\n";
	for (i = 1; i <= ct.numofbases; i++) {
		cout << ct.nucs[i];
	}
	cout << "\n";
	for (i = 1; i <= ct.numofbases; i++) {
		cout << ct.numseq[i];
	}
#endif

#if defined(VERBOSE_MODE)
	cout << "\n";
	cout << "Length of sequence: " << ct.numofbases << "\n";
#endif


	// Create a basepair type lookup matrix (AU = 5, GC = 6, GU = 7)
	short bpLookup[5][5] = {{0,0,0,0,0},
						    {0,0,0,0,5},
						    {0,0,0,6,0},
						    {0,0,6,0,7},
						    {0,5,0,7,0}};

	// Open original dotplot file
	ifstream inDotplotFile;
	inDotplotFile.open(indotplot.c_str());
	if (!inDotplotFile){cerr<<"ERROR: Could not open dotplot file "<<indotplot<<"\n";cin>>i;return -1;}

	// Create 2D matrices to store converted dotplot, DeltaG, and vmb/vext values
	short** convertedDotplot;
	short** dgArray;
	short** mbDotplot;
	convertedDotplot = new short* [(ct.numofbases)+1];
	dgArray = new short* [(ct.numofbases)+1];
	mbDotplot = new short* [(ct.numofbases)+1];	
	for(i = 0; i <= ct.numofbases; i++) {
		*(convertedDotplot + i) = new short[(ct.numofbases)+1];
		*(dgArray + i) = new short[(ct.numofbases)+1];
		*(mbDotplot + i) = new short[(ct.numofbases)+1];
	}
	
	// Clear the values stored in the matrices
	for (i = 0; i <= ct.numofbases; i++) {
		for (j = 0; j <= ct.numofbases; j++) {
			convertedDotplot[i][j] = 0;
			dgArray[i][j] = 0;
			mbDotplot[i][j] = 0;
		}
	}
	
	// Convert the dotplot and store it in the 2D matrix
	string s;
	int x, y, dpLineCount = 0, dgMin = 0;
	double z, vmb, vext;
	istringstream instream;
	while (getline(inDotplotFile,s)) {
		instream.clear();
		instream.str(s);
		if (instream >> x >> y >> z >> vmb >> vext) {
			convertedDotplot[x][y] = bpLookup[ct.numseq[x]][ct.numseq[y]];
			dgArray[x][y] = short(z*10);
			mbDotplot[x][y] = min(short(vmb*10),short(vext*10));
			// Copy these values to the other half of the matrices
			convertedDotplot[y][x] = bpLookup[ct.numseq[x]][ct.numseq[y]];
			dgArray[y][x] = short(z*10);
			mbDotplot[y][x] = min(short(vmb*10),short(vext*10));
			if(z*10 < dgMin) dgMin = short(z*10); // Finds and stores minimum free energy value from dotplot
		}
	}
	inDotplotFile.close();

	// Store cutoff DG value
	int dgCutoff = dgMin*(100-percent)/100;

	// Loop back through convertedDotplot and remove those values that have a DG value greater than dgCutoff
	for (i = 0; i <= ct.numofbases; i++) {
		for (j = 0; j <= ct.numofbases; j++) {
			if (dgArray[i][j] > dgCutoff) convertedDotplot[i][j] = 0;
		}
	}

#if defined(VERY_VERBOSE_MODE)
	// Display the converted dotplot
	cout << "\nConverted Dotplot:\n";
	for(i = 0; i <= ct.numofbases; i++) {
		for(j = 0; j <= ct.numofbases; j++) {
			cout << convertedDotplot[i][j];
		}
		cout << "\n";
	}
	cout << "\n";
#endif


	// Load NMR constraints
	short maxConLength = 0, numOfCon = 0, totalConLength = 0;
	
	ifstream inConFile;
	inConFile.open(inconstraints.c_str());
	if (!inConFile) {
		cerr << "Unable to open constraints file";
		cin >> i;
		exit(1);
	}
	
	// Loop through each line of the file to get count of constraints and maximum length of a constraint
	while (getline(inConFile,s)) {
		if(s.length() < 2){cerr<<"Constraints cannot be less than two basepairs in length!\n";cin>>i;return -1;}
		numOfCon++;
		totalConLength+=s.length();
		if(s.length() > maxConLength) maxConLength = s.length();
	}
	inConFile.close();
	inConFile.clear();

	
#if defined(VERBOSE_MODE)
	cout << numOfCon << " NMR constraints\n";
	cout << "Total length of NMR constraints: " << totalConLength << "\n";
	cout << "Length of largest NMR constraint: " << maxConLength << "\n\n";
#endif


	// Allocate 2D array for storing constraints (first column contains length of each constraint)
	// Also allocate two more copies of this array which will be used later for storing temporary
	// dotplot coordinates
	short** conArray;
	short** xCoords;
	short** yCoords;
	conArray = new short* [numOfCon+1];
	xCoords = new short* [numOfCon+1];
	yCoords = new short* [numOfCon+1];
	for(i = 0; i <= numOfCon; i++) {
		*(conArray + i) = new short[(maxConLength+1)];
		*(xCoords + i) = new short[(maxConLength+1)];
		*(yCoords + i) = new short[(maxConLength+1)];
	}

	// Re-open NMR constraints file and fill the 2D arrays (only fill Coords arrays with 0's)
	i = 0;
	inConFile.open(inconstraints.c_str());
	while (getline(inConFile,s)) {
		conArray[i][0] = s.length();
		xCoords[i][0] = 0;
		yCoords[i][0] = 0;
		for (j = 1; j <= s.length(); j++) {
			conArray[i][j] = atoi(s.substr(j-1,1).c_str());
			xCoords[i][j] = 0;
			yCoords[i][j] = 0;
		}
		i++;		
	}
	inConFile.close();
	inConFile.clear();

#if defined(VERY_VERBOSE_MODE)
	// Display the constraints array
	cout << "NMR Constraints:\n";
	for (i = 0; i < numOfCon; i++) {
		for (j = 1; j <= conArray[i][0]; j++) {
			cout << conArray[i][j];
		}
		cout << "\n";
	}
#endif


	// Create 1D bool array for storing the nucleotides that are being considered for a match
	bool* alreadyUsed = new bool[ct.numofbases+1];
	for (i = 0; i <= ct.numofbases; i++) {
		alreadyUsed[i] = false;
	}

	// Create storage container for matches
	vector<conMatch> matchVector;

	// Create 1D bool array for storing information about which constraints are symmetric
	bool* isSymmetric = new bool[numOfCon];
	for (i = 0; i < numOfCon; i++) {
		isSymmetric[i] = false;
	}

	// Test to see which constraints are symmetric
	bool tempEquality;
	for (i = 0; i < numOfCon; i++) {
		j = 0;
		tempEquality = true;
		while (j < conArray[i][0] && tempEquality) {
			if (conArray[i][conArray[i][0]-j] != conArray[i][j+1]) {
				tempEquality = false;
				continue;
			}
			j++;
			if (j == conArray[i][0]) {
				isSymmetric[i] = true;

#if defined(VERY_VERBOSE_MODE)
				cout << "Symmetric constraint: " << i << "\n";
#endif

			}
		}
	}
	short currConNum = 0;
	short currConPos = 1;

	// Ready to begin recursive constraint matching
	firstDotMatch(conArray,alreadyUsed,isSymmetric,totalConLength,numOfCon,xCoords,yCoords,&matchVector,
		&currConNum,&currConPos,convertedDotplot,mbDotplot,dgCutoff,ct.numofbases);

	// When finished with constraint matching
	conMatch tempMatchVector;
	cout << "Constraint matching completed!\nTotal number of matches: " << matchVector.size() << "\n";

#if defined(VERY_VERBOSE_MODE)
	// Display the matches
	for (i = 0; i < matchVector.size(); i++) {
		tempMatchVector = matchVector[i];
		for (j = 0; j < totalConLength; j++) {
			cout << tempMatchVector.coords[j*2] << "," << tempMatchVector.coords[j*2+1] << " ";
		}
		cout << "\n";
	}
#endif


	// Determine if refolding should occur
	if (matchVector.size() == 0) {
		cout << "No possible constraint matches - NAPSS will terminate.\n";
		cin >> i;
		return(0);
	}
	if (matchVector.size() > 100000) {
		cout << "Too many matches, refolding will not be performed.\n";
		cin >> i;
		return(0);
	}
	while (true) {
		cout << "\nContinue with refolding? (y/n) ";
		cin >> a;
		if (a == "n") return (0);
		if (a == "y") break;
	}
	cout << "\n";

	// Initialize ct's that will be used for storing all refolded structures and for calculating pseudoknot energies
	structure ct2;
	structure ctbroken;
	ct2.numofstructures = 0;
	ctbroken.numofstructures = 0;
	openseq(&ct2,inseq.c_str());
	openseq(&ctbroken,inseq.c_str());
	int start = 0;
	int count = 0;
	// Also create a 1D array for temporarily storing possible extensions and a control switch for that loop
	int* helixExtend = new int[ct.numofbases+1];
	bool extensionAdded = true;


	// Create 2D boolean array (triangular matrix with sides of length ct.numofbases+1) to store dotplot data
	ct.allocatetem();

	// Loop over all matches
	for (iter = 0; iter < matchVector.size(); iter++) {

#if defined(VERBOSE_MODE)
		if (iter % 100 == 0 && iter > 0) cout << iter << " refoldings completed so far.\n";
#endif

		// Initialize all positions to false:
		for (i = 0; i < ct.numofbases; i++) {
			for (j = i+1; j <= ct.numofbases; j++) {
				ct.tem[j][i] = false;
			}
		}

#if defined(VERY_VERBOSE_MODE)
		cout << "Non-conflicting match " << iter+1 << ":\n";
#endif

		// Copy original dotplot into tem** - change all non-zero convertedDotplot values to true
		for (i = 0; i < ct.numofbases; i++) {
			for (j = i+1; j <= ct.numofbases; j++) {
				if (convertedDotplot[j][i] != 0) ct.tem[j][i] = true;
			}
		}

		// Create new array to temporarily hold matched basepairs (initialize all to 0)
		int* tempbasepr = new int[ct.numofbases+1];
		for (i=0; i<=ct.numofbases; i++) tempbasepr[i]=0;
		
		// Trim the dotplot according to the base pairs of the current match
		int trimming_i, trimming_j;
		
		// Loop over all base pairs in the match combination
		for (iter2 = 0; iter2 < totalConLength; iter2++) {
			trimming_j = matchVector[iter].coords[(2*iter2)];
			trimming_i = matchVector[iter].coords[(2*iter2)+1];

			// Store coordinates in temporary basepair array
			tempbasepr[trimming_i] = trimming_j;
			tempbasepr[trimming_j] = trimming_i;

			// Change all tem values to false where (i or j) = (trimming_i or trimming_j)
			for (j = 0; j <= trimming_i-1; j++) {
				ct.tem[trimming_i][j] = false;
			}
			for (i = trimming_i + 1; i <= ct.numofbases; i++) {
				ct.tem[i][trimming_i] = false;
			}
			for (j = 0; j <= trimming_j-1; j++) {
				ct.tem[trimming_j][j] = false;
			}
			for (i = trimming_j + 1; i <=ct.numofbases; i++) {
				ct.tem[i][trimming_j] = false;
			}
		}

		// Search for regions of dotplot that can be excluded
		pseudodptrim(&ct, tempbasepr, &count);
		delete[] tempbasepr;

		// Do the structure prediction considering only the allowable basepairs from the trimmed dotplot
		dynamic(&ct,&data,maxtracebacks,99,windowsize);

#if defined(VERY_VERBOSE_MODE)
		cout << "Refolding yields " << ct.numofstructures << " structures.\n\n";
#endif

		// Reinsert the basepairs that match the NMR constraints
		for (i = 1; i <= ct.numofstructures; i++) {
			for (iter2 = 0; iter2 < totalConLength; iter2++) {
				trimming_j = matchVector[iter].coords[(2*iter2)];
				trimming_i = matchVector[iter].coords[(2*iter2)+1];
				ct.basepr[i][trimming_j] = trimming_i;
				ct.basepr[i][trimming_i] = trimming_j;
			}
		}

		// Check for helical extension possibilities
		for (i = 1; i <= ct.numofstructures; i++) {
			extensionAdded = true;
			while (extensionAdded) {
				extensionAdded = false;
				// Initialize the temporary array
				for (iter2 = 0; iter2 <= ct.numofbases; iter2++) {
					helixExtend[iter2] = 0;
				}
				
				iter2 = 1;
				while (iter2 <= ct.numofbases) {
					if (ct.basepr[i][iter2] < iter2) {
						iter2++;
						continue;
					}
					else {
						// Check to see if the nucleotides downstream can pair with each other
						if (iter2 > 1 && ct.basepr[i][iter2] < ct.numofbases &&
							ct.basepr[i][iter2-1] == 0 && ct.basepr[i][ct.basepr[i][iter2]+1] == 0) { 
								// Nucleotides are currently unpaired 
								if (convertedDotplot[iter2-1][ct.basepr[i][iter2]+1] !=0) { 
									// They could form a valid, stable pair
									if (helixExtend[iter2-1] == 0 && helixExtend[ct.basepr[i][iter2]+1] == 0) {
										// There are no other possible extensions to these nucleotides yet - store this one
										helixExtend[iter2-1] = ct.basepr[i][iter2]+1;
										helixExtend[ct.basepr[i][iter2]+1] = iter2-1;
									}
									else {
										// One or both nucleotides already have a possible extension - if both do, 
										// make no changes because one new pair shouldn't break up two other possible pairs
										if (helixExtend[iter2-1] != 0 && helixExtend[ct.basepr[i][iter2]+1] == 0) {
											// iter2 already has a possible helical extension to it - determine if the latest
											// pairing is more favorable in the original dotplot
											if (dgArray[iter2-1][helixExtend[iter2-1]] < 
												dgArray[iter2-1][ct.basepr[i][iter2]+1]) {
													helixExtend[helixExtend[iter2-1]] = 0;
													helixExtend[iter2-1] = ct.basepr[i][iter2]+1;
													helixExtend[ct.basepr[i][iter2]+1] = iter2-1;
											}
										}
										if (helixExtend[ct.basepr[i][iter2]+1] != 0 && helixExtend[iter2-1] == 0) {
											// the nucleotide opposite iter2 already has a possible helical extension to it - 
											// determine if the latest pairing is more favorable in the original dotplot
											if (dgArray[ct.basepr[i][iter2]+1][helixExtend[ct.basepr[i][iter2]+1]] < 
												dgArray[ct.basepr[i][iter2]+1][iter2-1]) {
													helixExtend[helixExtend[ct.basepr[i][iter2]+1]] = 0;
													helixExtend[ct.basepr[i][iter2]+1] = iter2-1;
													helixExtend[iter2-1] = ct.basepr[i][iter2]+1;
											}
										}
									}
								}
							}

							// Check to see if the nucleotides upstream can pair with each other
							if (ct.basepr[i][iter2+1] == 0 && ct.basepr[i][ct.basepr[i][iter2]-1] == 0) { 
								// Nucleotides are currently unpaired 
								if (convertedDotplot[iter2+1][ct.basepr[i][iter2]-1] !=0) { 
									// They could form a valid, stable pair
									if (helixExtend[iter2+1] == 0 && helixExtend[ct.basepr[i][iter2]-1] == 0) {
										// There are no other possible extensions to these nucleotides yet - store this one
										helixExtend[iter2+1] = ct.basepr[i][iter2]-1;
										helixExtend[ct.basepr[i][iter2]-1] = iter2+1;
									}
									else {
										// One or both nucleotides already have a possible extension - if both do, 
										// make no changes because one new pair shouldn't break up two other possible pairs
										if (helixExtend[iter2+1] != 0 && helixExtend[ct.basepr[i][iter2]-1] == 0) {
											// iter2 already has a possible helical extension to it - determine if the latest
											// pairing is more favorable in the original dotplot
											if (dgArray[iter2+1][helixExtend[iter2+1]] < 
												dgArray[iter2+1][ct.basepr[i][iter2]-1]) {
													helixExtend[helixExtend[iter2+1]] = 0;
													helixExtend[iter2+1] = ct.basepr[i][iter2]-1;
													helixExtend[ct.basepr[i][iter2]-1] = iter2+1;
											}
										}
										if (helixExtend[ct.basepr[i][iter2]-1] != 0 && helixExtend[iter2+1] == 0) {
											// the nucleotide opposite iter2 already has a possible helical extension to it - 
											// determine if the latest pairing is more favorable in the original dotplot
											if (dgArray[ct.basepr[i][iter2]-1][helixExtend[ct.basepr[i][iter2]-1]] < 
												dgArray[ct.basepr[i][iter2]-1][iter2+1]) {
													helixExtend[helixExtend[ct.basepr[i][iter2]-1]] = 0;
													helixExtend[ct.basepr[i][iter2]-1] = iter2+1;
													helixExtend[iter2+1] = ct.basepr[i][iter2]-1;
											}
										}
									}
								}
							}
							iter2++;
							continue;
					}
				}

				// Now loop through the temporary array and add any remaining potential pairs to the final structure
				for (iter2 = 1; iter2 <= ct.numofbases; iter2++) {
					if (helixExtend[iter2] != 0) {
						// Adding a pair - flip the switch
						extensionAdded = true;
						ct.basepr[i][iter2] = helixExtend[iter2];
					}
				}
			}
		}





		// Copy current ct to new ct file that is used for appending the results from all matches
		for (i = 1; i <= ct.numofstructures; i++) {
			ct2.checknumberofstructures();
			ct2.numofstructures = ct2.numofstructures+1;
			ctbroken.checknumberofstructures();
			ctbroken.numofstructures = ctbroken.numofstructures+1;			
			for (j = 1; j <= ct.numofbases; j++) {
				ct2.basepr[i+start][j] = ct.basepr[i][j];
			}
			strcpy(ct2.ctlabel[i+start], ct.ctlabel[i]);
		}
		start += ct.numofstructures;
	}

#if defined(VERBOSE_MODE)
	cout << "Refolding yields " << ct2.numofstructures << " total structures.\n\n";
#endif
	
	delete[] helixExtend;

	// Re-calculate free energy with modified efn2
	for (i=1;i<=ct2.numofstructures;i++){
		efn2mod(&data, &ct2, i, false, &ctbroken);

#if defined(VERBOSE_MODE)
		if (i % 100 == 0 && i > 0) cout << i << " free energies calculated so far.\n";
#endif

		// Re-insert pseudoknot stems that were broken for energy calculation
		for (j=1; j<=ct2.numofbases; j++){
			if (ctbroken.basepr[0][j]!=0) ct2.basepr[i][j] = ctbroken.basepr[0][j];
		}
	}

	// Re-sort the structures according to efn2mod-calculated energy
	sortstructures(&ct2);

	// Calculate value of cutoff for output structures
	if (cutoff !=0) cutoff = (100-cutoff)*ct2.energy[1]/100;

	// Re-output ct file
	cout << ctout2(&ct2, cutoff, outct.c_str()) << " structures are within the specified cutoff percentage.\n\n";

//  Optional: output paired-positions text file (for structure viewing in PseudoViewer3)
//	Note: this outputs one concatenated file, individual structures must be manually cut from this file
//  and placed in a new text file for PseudoViewer3 to read it.
	if(outpairs != "") pairout(&ct2, cutoff, outpairs.c_str());

	// Insert steps to clean up allocated memory
	matchVector.clear();
	delete[] alreadyUsed;
	delete[] isSymmetric;
	
	for (i=0;i<=ct.numofbases;i++){
		delete[]convertedDotplot[i];
		delete[]dgArray[i];
		delete[]mbDotplot[i];
	}
	delete[]convertedDotplot;
	delete[]dgArray;
	delete[]mbDotplot;
	
	for (i=0;i<=numOfCon;i++){
		delete[]conArray[i];
		delete[]xCoords[i];
		delete[]yCoords[i];
	}
	delete[]conArray;
	delete[]xCoords;
	delete[]yCoords;

	cout << "NAPSS completed.\n";
	return 0;
}

void getdat(char *loop, char *stackf, char *tstackh, char *tstacki,
	char *tloop, char *miscloop, char *danglef, char *int22,
	char *int21,char *coax, char *tstackcoax,
    char *coaxstack, char *tstack, char *tstackm, char *triloop,
    char *int11, char *hexaloop, char *tstacki23, char *tstacki1n,
	char *datapath, bool isRNA)
{
	if( !isRNA) strcat( datapath,"dna");
	strcpy (loop,datapath);
	strcpy (stackf,datapath);
	strcpy (tstackh,datapath);
	strcpy (tstacki,datapath);
	strcpy (tloop,datapath);
	strcpy (miscloop,datapath);
	strcpy (danglef,datapath);
	strcpy (int22,datapath);
	strcpy (int21,datapath);
	strcpy (triloop,datapath);
	strcpy (coax,datapath);
	strcpy (tstackcoax,datapath);
	strcpy (coaxstack,datapath);
	strcpy (tstack,datapath);
	strcpy (tstackm,datapath);
	strcpy (int11,datapath);
	strcpy (hexaloop,datapath);
	strcpy (tstacki23,datapath);
	strcpy (tstacki1n,datapath);

	strcat (loop,"loop.dat");
	strcat (stackf,"stack.dat");
	strcat (tstackh,"tstackh.dat");
	strcat (tstacki,"tstacki.dat");
	strcat (tloop,"tloop.dat");
	strcat (miscloop,"miscloop.dat");
	strcat (danglef,"dangle.dat");
	strcat (int22,"int22.dat");
	strcat (int21,"int21.dat");
	strcat (triloop,"triloop.dat");
	strcat (coax,"coaxial.dat");
	strcat (tstackcoax,"tstackcoax.dat");
	strcat (coaxstack,"coaxstack.dat");
	strcat (tstack,"tstack.dat");
	strcat (tstackm,"tstackm.dat");
	strcat (int11,"int11.dat");
	strcat (hexaloop,"hexaloop.dat");
	strcat (tstacki23,"tstacki23.dat");
	strcat (tstacki1n,"tstacki1n.dat");
}

void errmsg(int err,int erri) {

if (err==30) {
	cout << "End Reached at traceback #"<<erri<<"\n";
   exit(1);
}
if (err==100) {
	cout << "error # "<<erri;
   exit(1);
}
switch (err) {
	case 1:
   	cout << "Could not allocate enough memory";
      break;
   case 2:
   	cout << "Too many possible base pairs";
      break;
   case 3:
   	cout << "Too many helixes in multibranch loop";
   case 4:
   	cout << "Too many structures in CT file";
   default:
   	cout << "Unknown error";
}
cin >> err;
exit(1);
return;
}


void firstDotMatch(short** conArray, bool* alreadyUsed, bool* isSymmetric, short totalConLength,
				   short numOfCon, short** xCoords, short** yCoords, vector<conMatch>* matchVector,
				   short* pCurrConNum, short* pCurrConPos, short** convertedDotplot, short** mbDotplot,
				   int dgCutoff, short sequenceLength) {

	short i, j;
	short currConNum = *(pCurrConNum);
	short currConPos = *(pCurrConPos);

#if defined(DEBUG_MODE)
	cout << currConNum << " " << currConPos << "\n";
	cin >> i;
#endif

	// Loop through x,y of dotplot to find potential match for first position of the current constraint
	// Skip all values that have already been used by matches for previous constraints
	for (j = 1; j < sequenceLength+1; j++) {
		if (alreadyUsed[j] == true) continue;

		for (i = 1; i < sequenceLength+1; i++) {
			if (alreadyUsed[i] == true) continue;

			// If constraint is symmetric, only consider matches from upper-right half of dotplot
			if (isSymmetric[currConNum]==false || (isSymmetric[currConNum]==true && i<j)) {

				// Check the first value only here - if this value matches,
				// then call the recursive method, passing it the proper parameters
				if (conArray[currConNum][1] == convertedDotplot[i][j]) {
					xCoords[currConNum][1] = i;
					yCoords[currConNum][1] = j;
					alreadyUsed[i] = true;
					alreadyUsed[j] = true;
					*(pCurrConPos) = 2;
					recursiveMatch(conArray, alreadyUsed, totalConLength, numOfCon, xCoords, yCoords, 
						matchVector, pCurrConNum, pCurrConPos, convertedDotplot, mbDotplot, dgCutoff, sequenceLength, 
						isSymmetric);
					currConNum = *(pCurrConNum);
					currConPos = *(pCurrConPos);
				}
			} 
		}
	}
	// Once we've looped through all coordinates for this base pair, step back values for constraint 
	// number and position, flip switches for last base pair, and return
	currConNum--;
	if (currConNum > -1) {
		currConPos = conArray[currConNum][0];
		alreadyUsed[xCoords[currConNum][currConPos]] = false;
		alreadyUsed[yCoords[currConNum][currConPos]] = false;
	}
	*(pCurrConNum) = currConNum;
	*(pCurrConPos) = currConPos;
	return;
}

void recursiveMatch(short** conArray, bool* alreadyUsed, short totalConLength, short numOfCon,
					short** xCoords, short** yCoords, vector<conMatch>* matchVector, 
					short* pCurrConNum, short* pCurrConPos, short** convertedDotplot, short** mbDotplot,
					int dgCutoff, short sequenceLength, bool* isSymmetric) {

	short i,j,k,xIter,yIter,start,stop;
	short currConNum = *(pCurrConNum);
	short currConPos = *(pCurrConPos);
	conMatch* tempMatch;
	string a;
	bool searchX;  // Switch to indicate in which direction we're searching
	bool searchExtCoax; // Switch to indicate whether we're currently searching for a coaxial stack in 
						// the external direction
	

subroutine: // Had to use this so that part of the code could return here from a nested "for" loop

	// This function will loop until we reach the end of the current constraint, 
	// or until no more matches are found using the current nucleotides for the
	// first base pair of this constraint
	while (currConPos > 1) {

#if defined(DEBUG_MODE)
		cout << currConNum << " " << currConPos << " " << xCoords[currConNum][currConPos-1] << 
			" " << yCoords[currConNum][currConPos-1] << " " << xCoords[currConNum][currConPos] <<
			" " << yCoords[currConNum][currConPos] << "\n";
		cin >> i;
#endif

		// Check if we've matched all the base pairs for the current constraint
		if (currConPos == conArray[currConNum][0]+1) {

			// Check if this is the final constraint
			if (currConNum == numOfCon-1) {

				// If so, we've completed a match - compile 1D short array of dotplot 
				// coordinates and append this to matchVector
				tempMatch = new conMatch;
				tempMatch->coords = new short[totalConLength*2];
				i = 0; // current base pair to be written
				j = 0; // takes place of currConNum
				k = 1; // takes place of currConPos
				
				while (i < totalConLength) {
					if (k > conArray[j][0]) {j++; k=1; continue;} // reached end of current constraint
					tempMatch->coords[i*2] = xCoords[j][k];
					tempMatch->coords[i*2+1] = yCoords[j][k];
					i++;
					k++;
				}

/*				cout << "Match found!: ";
				for (i = 0; i < totalConLength*2; i++) {
					cout << tempMatch[i] << " ";
				}
				cout << "\n";
*/


#if defined(ALREADY_USED_CHECK)
				// Verify that the number of "true" alreadyUsed values = twice the number of basepairs in all constraints
				j = 0;
				for (i = 1; i<=sequenceLength; i++){
					if (alreadyUsed[i]) j++;
				}
				if (j!=totalConLength*2){cerr<<"ERROR: Discrepancy in alreadyUsed array.\n";}
#endif


				matchVector->push_back(*tempMatch);
				delete[] tempMatch;

				if (matchVector->size() % 50000 == 0) {

#if defined(VERBOSE_MODE)
					// Display a current count of matches
					cout<<matchVector->size()<<" ";
					for (i=0; i<numOfCon; i++){
						cout <<xCoords[i][1]<<","<<yCoords[i][1]<<" ";
					}
					cout << "\n";
#endif

				}				
				if (matchVector->size() == 100000) {
					while (true) {
						cout << "Large number of matches being found - continue matching? (y/n) ";
						cin >> a;
						if (a == "n") exit(0);
						if (a == "y") break;
					}
				}

				// Now step back two base pairs, clear final position in the Coords arrays,
				// flip the alreadyUsed values of last 2 positions, and continue searching
				currConPos=currConPos-2;
				alreadyUsed[xCoords[currConNum][currConPos+1]] = false;
				alreadyUsed[yCoords[currConNum][currConPos+1]] = false;
				xCoords[currConNum][currConPos+1] = 0;
				yCoords[currConNum][currConPos+1] = 0;
				alreadyUsed[xCoords[currConNum][currConPos]] = false;
				alreadyUsed[yCoords[currConNum][currConPos]] = false;
				continue;
			}

			// If this is not the last constraint, we need to move on and match the first position 
			// of the next constraint
			currConNum++;
			currConPos = 1;
			*(pCurrConNum) = currConNum;
			*(pCurrConPos) = currConPos;
			firstDotMatch(conArray,alreadyUsed,isSymmetric,totalConLength,numOfCon,xCoords,yCoords,matchVector,
				pCurrConNum,pCurrConPos,convertedDotplot,mbDotplot,dgCutoff,sequenceLength);
			// Once this returns, need to continue the search for matches to this constraint
			currConNum = *(pCurrConNum);
			currConPos = *(pCurrConPos);
			continue;
		}
		
		// Check if we're still in-bounds on the next position across and down.
		// If we're not, step back, flip switches, and continue
		if (xCoords[currConNum][currConPos-1]-1 < 1 || yCoords[currConNum][currConPos-1]+1 > sequenceLength) {
			currConPos--;
			alreadyUsed[xCoords[currConNum][currConPos]] = false;
			alreadyUsed[yCoords[currConNum][currConPos]] = false;
			continue;
		}

		// If this is the first or last step, or if the previous basepair is unstable as a closing pair for 
		// external and mb loops, only look for non-bulged/non-stacked match on current step
		if (currConPos==2 || currConPos==(conArray[currConNum][0]) || 
			mbDotplot[xCoords[currConNum][currConPos-1]][yCoords[currConNum][currConPos-1]] >= dgCutoff /*||
			xCoords[currConNum][currConPos-2]-xCoords[currConNum][currConPos-1]>1 ||
			yCoords[currConNum][currConPos-1]-yCoords[currConNum][currConPos-2]>1*/) {

				if (xCoords[currConNum][currConPos] != 0 || yCoords[currConNum][currConPos] !=0) {
					// This would indicate that this basepair has already been checked
					// Reset coords, step back, flip switches, and continue
					xCoords[currConNum][currConPos] = 0;
					yCoords[currConNum][currConPos] = 0;
					currConPos--;
					alreadyUsed[xCoords[currConNum][currConPos]] = false;
					alreadyUsed[yCoords[currConNum][currConPos]] = false;
					continue;
				}

				if (conArray[currConNum][currConPos]==convertedDotplot[xCoords[currConNum][currConPos-1]-1]
				[yCoords[currConNum][currConPos-1]+1] 
				&& alreadyUsed[xCoords[currConNum][currConPos-1]-1] == false
					&& alreadyUsed[yCoords[currConNum][currConPos-1]+1] == false) {

						// Match found at lastX-1,lastY+1
						xCoords[currConNum][currConPos] = xCoords[currConNum][currConPos-1]-1;
						yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1]+1;
						alreadyUsed[xCoords[currConNum][currConPos-1]-1] = true;
						alreadyUsed[yCoords[currConNum][currConPos-1]+1] = true;
						currConPos++;
						continue;
				}	

				// If not a match, step back, flip switches, and continue
				currConPos--;
				alreadyUsed[xCoords[currConNum][currConPos]] = false;
				alreadyUsed[yCoords[currConNum][currConPos]] = false;
				continue;
		}

		// Next step is in-bounds - find out where the search left off
	
		if (xCoords[currConNum][currConPos] == 0 && yCoords[currConNum][currConPos] == 0){
			// This is the first time searching forward from this position
			// Set coords values to -1,-1 to act as placeholder
			xCoords[currConNum][currConPos] = -1;
			yCoords[currConNum][currConPos] = -1;
			searchX = true;
			searchExtCoax = false;
			start = xCoords[currConNum][currConPos-1] - 1;
		}

		else if (yCoords[currConNum][currConPos] == yCoords[currConNum][currConPos-1]+1 &&
			xCoords[currConNum][currConPos] <= xCoords[currConNum][currConPos-1]-1) {
				// We left off searching in the X direction

				if (xCoords[currConNum][currConPos-1] -1 < 1){ 
					// Last found a match at left edge of dot plot - start searching in Y direction
					xCoords[currConNum][currConPos] = -1;
					yCoords[currConNum][currConPos] = -1;
					continue;
				}

				// Otherwise, continue searching in X direction
				searchX = true;
				searchExtCoax = false;
				start = xCoords[currConNum][currConPos] - 1;
		}

		else if (yCoords[currConNum][currConPos] == -1) {
			// Search in X direction was completed

			if (yCoords[currConNum][currConPos-1] + 2 > sequenceLength) {
				// Cannot search in Y direction because last match is near bottom edge of dot plot
				// Time to search for external coaxial stacks
				yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1]+1;
				xCoords[currConNum][currConPos] = sequenceLength+1;
				continue;
			}
			
			// Otherwise, start at +2 because +1 was already covered in searchX
			searchX = false;
			searchExtCoax = false;
			start = yCoords[currConNum][currConPos-1] + 2;
		}

		else if (xCoords[currConNum][currConPos] == xCoords[currConNum][currConPos-1]-1 &&
			yCoords[currConNum][currConPos] >= yCoords[currConNum][currConPos-1]+1) {
				// We left off searching in the Y direction

				if (yCoords[currConNum][currConPos] == sequenceLength){ 
					// Last found a match at bottom edge of dot plot - start searching for external coaxial stacks
					yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1]+1;
					xCoords[currConNum][currConPos] = sequenceLength+1;
					continue;
				}

				// Otherwise, continue searching in Y direction
				searchX = false;
				searchExtCoax = false;
				start = yCoords[currConNum][currConPos] + 1;
		}

		else if (yCoords[currConNum][currConPos] == yCoords[currConNum][currConPos-1]+1 &&
			xCoords[currConNum][currConPos] > yCoords[currConNum][currConPos]) {
				// We're currently searching for a external coaxial stack in the X direction
				searchX = true;
				searchExtCoax = true;
				start = xCoords[currConNum][currConPos] - 1;
		}

		else if (xCoords[currConNum][currConPos] == xCoords[currConNum][currConPos-1]-1 &&
			yCoords[currConNum][currConPos] != -1 && 
			yCoords[currConNum][currConPos] < xCoords[currConNum][currConPos]) {
				// We're currently searching for a external coaxial stack in the Y direction
				searchX = false;
				searchExtCoax = true;
				start = yCoords[currConNum][currConPos] + 1;
		}


		else {
			cerr<<"ERROR: Unclassified looping condition in constraint matching.\n"
				<<xCoords[currConNum][currConPos]<<" "<<yCoords[currConNum][currConPos]<<"\n";			
			cin>>i;break;
		}

		// Enter this loop if we need to search in the X direction
		if (searchX && !searchExtCoax) {

			if (alreadyUsed[yCoords[currConNum][currConPos-1]+1]) {
				// Y coordinate has already been used, skip searching in X direction
				xCoords[currConNum][currConPos] = -1;
				yCoords[currConNum][currConPos] = -1;
				continue;
			}
			
			// If we're in the upper-right half of the dot plot...
			// Search the entire subspace of dotplot where next coordinate is 
			// x = {start,...,lastY + 1}, y = lastY+ 1
			if (xCoords[currConNum][currConPos-1] > yCoords[currConNum][currConPos-1]) {
				stop = (yCoords[currConNum][currConPos-1]);
			}

			// If we're in the lower-left half of the dot plot...
			// Search the entire subspace of dotplot where next coordinate is 
			// x = {start,...,1}, y = lastY + 1
			else stop = 1;

#if defined(DEBUG_MODE)
			cout << start << " " << stop << "\n";
			cin >> i;
#endif
			
			for (xIter=0;(start-xIter) >= stop; xIter++) {
				// Don't include certain values of xCoord
				if (alreadyUsed[start-xIter] ||
					(xCoords[currConNum][currConPos-1]-(start-xIter) > 2 &&
					xCoords[currConNum][currConPos-1]-(start-xIter) < 8)) continue;

#if defined(DEBUG_MODE)
				cout << start-xIter;
#endif

				if (conArray[currConNum][currConPos]==
					convertedDotplot[start-xIter][yCoords[currConNum][currConPos-1]+1] && 
					(start-xIter == xCoords[currConNum][currConPos-1] - 1 ||
					mbDotplot[start-xIter][yCoords[currConNum][currConPos-1]+1] <= dgCutoff)) {
						// Match found: store coords, flip switches, step forward, and continue
						xCoords[currConNum][currConPos] = start-xIter;
						yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1]+1;
						alreadyUsed[start-xIter] = true;
						alreadyUsed[yCoords[currConNum][currConPos-1]+1] = true;
						currConPos++;
						// Can't simply continue the "while" loop because we're in a nested "for" loop
						goto subroutine;
				}
			}

			// If we reach this point, we've iterated through all X values at this Y coordinate
			// Time to continue the search in the Y direction
			xCoords[currConNum][currConPos] = -1;
			yCoords[currConNum][currConPos] = -1;
			continue;
		}

		// Enter this loop if we need to search in the y direction 
		else if (!searchX && !searchExtCoax) {

			if (alreadyUsed[xCoords[currConNum][currConPos-1]-1]) {
				// X coordinate has already been used
				
				// Check if we're in the lower-left half of the dot plot - if so, need to search for external coaxial stack
				if (xCoords[currConNum][currConPos-1] < yCoords[currConNum][currConPos-1]) {
					yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1] + 1;
					xCoords[currConNum][currConPos] = sequenceLength + 1;
					continue;
				}
				
				// We're in the upper-right half of the dot plot - reset the coords, step back, flip switches, and continue
				xCoords[currConNum][currConPos] = 0;
				yCoords[currConNum][currConPos] = 0;
				currConPos--;
				alreadyUsed[xCoords[currConNum][currConPos]] = false;
				alreadyUsed[yCoords[currConNum][currConPos]] = false;
				continue;
			}

			// If we're in the upper-right half of the dot plot...
			// Search the entire subspace of dotplot where next coordinate is 
			// x = lastX - 1, y = {start,...,lastX - 1}
			if (xCoords[currConNum][currConPos-1] > yCoords[currConNum][currConPos-1]) {
				stop = (xCoords[currConNum][currConPos-1]-1);
			}

			// If we're in the lower-left half of the dot plot...
			// Search the entire subspace of dotplot where next coordinate is 
			// x = lastX - 1, y = {start,...,sequenceLength}
			else stop = sequenceLength;

#if defined(DEBUG_MODE)
			cout << start << " " << stop << "\n";
			cin >> i;
#endif

			for (yIter=0;(start+yIter) <= stop; yIter++) {
				// Don't include certain values of yCoord
				if (alreadyUsed[start+yIter] ||
					((start+yIter)-yCoords[currConNum][currConPos-1] > 2 &&
					(start+yIter)-yCoords[currConNum][currConPos-1] < 8)) continue;

#if defined(DEBUG_MODE)
				cout << start+yIter;
#endif

				if (conArray[currConNum][currConPos]==
					convertedDotplot[xCoords[currConNum][currConPos-1]-1][start+yIter] &&
					mbDotplot[xCoords[currConNum][currConPos-1]-1][start+yIter] <= dgCutoff) {

						// Match found: store coords, flip switches, step forward, and continue
						xCoords[currConNum][currConPos] = xCoords[currConNum][currConPos-1]-1;
						yCoords[currConNum][currConPos] = start+yIter;
						alreadyUsed[xCoords[currConNum][currConPos-1]-1] = true;
						alreadyUsed[start+yIter] = true;
						currConPos++;
						// Can't simply continue the "while" loop because we're in a nested "for" loop
						goto subroutine;
				}

			}

			// If we reach this point, we've iterated through all x and y values for this base pair

			// Check if we're in the lower-left half of the dot plot - if so, need to search for external coaxial stack
			if (xCoords[currConNum][currConPos-1] < yCoords[currConNum][currConPos-1]) {
				yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1] + 1;
				xCoords[currConNum][currConPos] = sequenceLength + 1;
				continue;
			}

			// We're in the upper-right half of the dot plot - reset the coords, step back, flip switches, and continue
			xCoords[currConNum][currConPos] = 0;
			yCoords[currConNum][currConPos] = 0;
			currConPos--;
			alreadyUsed[xCoords[currConNum][currConPos]] = false;
			alreadyUsed[yCoords[currConNum][currConPos]] = false;
			continue;
		}

		// Enter this loop if we're searching for an external coaxial stack in the X direction
		else if (searchX && searchExtCoax) {

			if (alreadyUsed[yCoords[currConNum][currConPos-1]+1]) {
				// Y coordinate has already been used, skip searching in X direction
				xCoords[currConNum][currConPos] = xCoords[currConNum][currConPos-1]-1;
				yCoords[currConNum][currConPos] = 0;
				continue;
			}
			
			stop = (yCoords[currConNum][currConPos-1]+1);
			
#if defined(DEBUG_MODE)
			cout << start << " " << stop << "\n";
			cin >> i;
#endif
			
			for (xIter=0;(start-xIter) >= stop; xIter++) {
				// Don't include certain values of xCoord
				if (alreadyUsed[start-xIter]) continue;

#if defined(DEBUG_MODE)
				cout << start-xIter;
#endif

				if (conArray[currConNum][currConPos]==
					convertedDotplot[start-xIter][yCoords[currConNum][currConPos-1]+1] && 
					mbDotplot[start-xIter][yCoords[currConNum][currConPos-1]+1] <= dgCutoff) {
						// Match found: store coords, flip switches, step forward, and continue
						xCoords[currConNum][currConPos] = start-xIter;
						yCoords[currConNum][currConPos] = yCoords[currConNum][currConPos-1]+1;
						alreadyUsed[start-xIter] = true;
						alreadyUsed[yCoords[currConNum][currConPos-1]+1] = true;
						currConPos++;
						// Can't simply continue the "while" loop because we're in a nested "for" loop
						goto subroutine;
				}
			}

			// If we reach this point, we've iterated through all x values at this y coordinate
			// Time to continue the external coaxial stack search in the y direction
			xCoords[currConNum][currConPos] = xCoords[currConNum][currConPos-1]-1;
			yCoords[currConNum][currConPos] = 0;
			continue;
		}

		// Enter this loop if we need to search for external coaxial stacks in the y direction 
		else if (!searchX && searchExtCoax) {

			if (alreadyUsed[xCoords[currConNum][currConPos-1]-1]) {
				// X coordinate has already been used, time to reset, step back, flip, and continue 
				xCoords[currConNum][currConPos] = 0;
				yCoords[currConNum][currConPos] = 0;
				currConPos--;
				alreadyUsed[xCoords[currConNum][currConPos]] = false;
				alreadyUsed[yCoords[currConNum][currConPos]] = false;
				continue;
			}
			
			stop = (xCoords[currConNum][currConPos-1]-1);

#if defined(DEBUG_MODE)
			cout << start << " " << stop << "\n";
			cin >> i;
#endif

			for (yIter=0;(start+yIter) <= stop; yIter++) {
				// Don't include certain values of yCoord
				if (alreadyUsed[start+yIter]) continue;

#if defined(DEBUG_MODE)
				cout << start+yIter;
#endif

				if (conArray[currConNum][currConPos]==
					convertedDotplot[xCoords[currConNum][currConPos-1]-1][start+yIter] &&
					mbDotplot[xCoords[currConNum][currConPos-1]-1][start+yIter] <= dgCutoff) {

						// Match found: store coords, flip switches, step forward, and continue
						xCoords[currConNum][currConPos] = xCoords[currConNum][currConPos-1]-1;
						yCoords[currConNum][currConPos] = start+yIter;
						alreadyUsed[xCoords[currConNum][currConPos-1]-1] = true;
						alreadyUsed[start+yIter] = true;
						currConPos++;
						// Can't simply continue the "while" loop because we're in a nested "for" loop
						goto subroutine;
				}
			}

			// If we reach this point, we've iterated through all x and y values for this base pair
			// Time to reset the coords, step back, flip switches, and continue
			xCoords[currConNum][currConPos] = 0;
			yCoords[currConNum][currConPos] = 0;
			currConPos--;
			alreadyUsed[xCoords[currConNum][currConPos]] = false;
			alreadyUsed[yCoords[currConNum][currConPos]] = false;
			continue;
		}
	}

	// The current position has gone back to the first base pair of a constraint
	// Flip switches and return
	alreadyUsed[xCoords[currConNum][currConPos]] = false;
	alreadyUsed[yCoords[currConNum][currConPos]] = false;
	*(pCurrConNum) = currConNum;
	*(pCurrConPos) = currConPos;
	return;
}

// Search dotplot for complicated pseudoknot folds to exclude
void pseudodptrim(structure *ct, int *tempbp, int * count) {
	int i, iter, iter2, iter3, init=1, a=0, b=0, c=0, d=0, e=0, f=0, j=0;
	bool rightSingle, findCF, leftSingle, findAD;
	while(init <= ct->numofbases) { 
		if(tempbp[init] > init) {
			for(i = 1; i < init; i++) {
				if(tempbp[i] > init && tempbp[i] < tempbp[init]) {
					// Define stems of pseudoknotted region
					j = tempbp[init];
					e = tempbp[i];
					b = tempbp[j];

					rightSingle = false;
					findCF = false;
					f = e+1; // Initial guess
					while( findCF == false && rightSingle == false) {
						if( tempbp[f] != 0 && tempbp[f] < e) {
							c = tempbp[f];
							findCF = true;
						}
						else {
							f++;
							if( f >= j) {
								rightSingle = true;
								f = j;
								c= b;
							}
						}
					}

					leftSingle = false;
					findAD = false;
					a = b - 1; // Initial guess
					while( findAD == false && leftSingle == false) {
						if( tempbp[a] != 0 && tempbp[a] >  c) {
							d = tempbp[a];			
							findAD = true;
						}
						else {
							a--;
							if( a <= i) {
								leftSingle = true;
								a = i;
								d = e;
							}
						}
					}

#if defined(VERY_VERBOSE_MODE)
					cout << "Pseudoknot detected: " << 
						i << "-" << e << "," << a << "-" << d << " ; " <<
						c << "-" << f << "," << b << "-" << j << "\n";
					}
#endif


					// Remove all possible base pairs that would involved pairing between the 2 gap regions
					int count = 0;
					for (iter = a+1; iter < b; iter++) {
						for (iter2 = e+1; iter2 < f; iter2++) {
							if (ct->tem[iter2][iter] == true) count++;
							ct->tem[iter2][iter] = false;
						}
					}

#if defined(VERY_VERBOSE_MODE)
					if (count > 0) cout << "Excluded pairings between gaps: " << count << "\n";
#endif

					// Scan gaps for contained secondary structures - forbid this region from pairing outside the gap
					// a-b gap:
					count = 0;
					iter = a+1;
					while (iter < b ) {
						if (tempbp[iter] > iter && tempbp[iter] > a && tempbp[iter] < b) { 
							// Gap has self-contained secondary structure
							// First condition is to prevent endless looping scenario
							// Forbid all pairings for iter-tempbp[iter] to anything outside the gap
							for (iter2 = iter; iter2 <= tempbp[iter]; iter2++) {
								for (iter3 = 1; iter3 <= a; iter3++) {
									if (ct->tem[iter2][iter3] == true) count++;
									ct->tem[iter2][iter3] = false;
								}
								for (iter3 = b; iter3 <= ct->numofbases; iter3++) {
									if (ct->tem[iter3][iter2] == true) count++;
									ct->tem[iter3][iter2] = false;
								}
							}
							// Continue loop at end of contained secondary structure
							iter = (tempbp[iter] + 1);
							continue;
						}
						iter++;
						continue;
					}
#if defined(VERY_VERBOSE_MODE)
					if (count > 0) cout << "Pairings forbidden from structures in a-b gap: " << count << "\n";
#endif
					
					// c-d gap:
					count = 0;
					iter = c+1;
					while (iter < d ) {
						if (tempbp[iter] > iter && tempbp[iter] > c && tempbp[iter] < d) { 
							// Gap has self-contained secondary structure
							// First condition is to prevent endless looping scenario
							// Forbid all pairings for iter-tempbp[iter] to anything outside the gap
							for (iter2 = iter; iter2 <= tempbp[iter]; iter2++) {
								for (iter3 = 1; iter3 <= c; iter3++) {
									if (ct->tem[iter2][iter3] == true) count++;
									ct->tem[iter2][iter3] = false;
								}
								for (iter3 = d; iter3 <= ct->numofbases; iter3++) {
									if (ct->tem[iter3][iter2] == true) count++;
									ct->tem[iter3][iter2] = false;
								}
							}
							// Continue loop at end of contained secondary structure
							iter = (tempbp[iter] + 1);
							continue;
						}
						iter++;
						continue;
					}
#if defined(VERY_VERBOSE_MODE)
					if (count > 0) cout << "Pairings forbidden from structures in c-d gap: " << count << "\n";
#endif
					
					// e-f gap:
					count = 0;
					iter = e+1;
					while (iter < f ) {
						if (tempbp[iter] > iter && tempbp[iter] > e && tempbp[iter] < f) { 
							// Gap has self-contained secondary structure
							// First condition is to prevent endless looping scenario
							// Forbid all pairings for iter-tempbp[iter] to anything outside the gap
							for (iter2 = iter; iter2 <= tempbp[iter]; iter2++) {
								for (iter3 = 1; iter3 <= e; iter3++) {
									if (ct->tem[iter2][iter3] == true) count++;
									ct->tem[iter2][iter3] = false;
								}
								for (iter3 = f; iter3 <= ct->numofbases; iter3++) {
									if (ct->tem[iter3][iter2] == true) count++;
									ct->tem[iter3][iter2] = false;
								}
							}
							// Continue loop at end of contained secondary structure
							iter = (tempbp[iter] + 1);
							continue;
						}
						iter++;
						continue;
					}
#if defined(VERY_VERBOSE_MODE)
					if (count > 0) cout << "Pairings forbidden from structures in e-f gap: " << count << "\n";
#endif
				
					// Continue loop at 3' end of pseudoknot
					init = j;
					continue;
				}
			}				
		}
		// If no pseudoknot detected, step up one position and recheck
		init++;
		continue;
	}
	return;
}
/*	Function efn2mod - pseudoknot-capable version of efn2

   Calculates the free energy of each structure in a structure called structure.

   Structures can contain certain types of pseudoknots (those allowed in the Dirks and Pierce algorithm)

   Structnum indicates which structure to calculate the free energy
   	the default, 0, indicates "all" structures

	simplemb = true means that the energy function for mb loops is linear,
		so identical to the dynamic programming algorithms.
		The default is false, so logarithmic in the number of unpaired nucleotides.

	This new version of efn2 is now stacking aware -- if ct->stacking==true, it
		will use the stacking information stored in ct, rather than inferring stacking
		by free energy minimization.
*/
void efn2mod(datatable *data,structure *ct, int structnum, bool simplemb, structure *ctbroken) {
	int i,j,k,open,null,stz,count,sum,ip,jp;
	stackstruct stack;
	forceclass fce(ct->numofbases);
	int start, stop;

	ofstream out;

#if defined(debugmode)
	char filename[maxfil],temp[maxfil];
	int tempsum;
#endif
	
	/*	stack = a place to keep track of where efn2 is located in a structure
	inter = indicates whether there is an intermolecular interaction involved in
	a multi-branch loop
	*/

	stack.sp = 0;  //set stack counter
	
	if (ct->intermolecular) {//this indicates an intermolecular folding
		for (i=0;i<3;i++) {
			forceinterefn(ct->inter[i],ct,&fce);
		}
	}
	
	if (structnum!=0) {
		start = structnum;
		stop = structnum;
	}
	else {
		start = 1;
		stop = ct->numofstructures;
	}
	
	for (count=start;count<=stop;count++){//one structure at a time

#if defined(debugmode) //open an output file for debugging
		strcpy(filename,"efn2dump");
		itoa(count,temp,10);
		strcat(filename,temp);
		strcat(filename,".out");
		out.open(filename);
#endif

		// Reset the zero-order ctbroken basepr array for storing pairing information of all broken pseudoknot stems
		for (j = 1; j <=ctbroken->numofbases; j++){
			ctbroken->basepr[0][j] = 0;
		}

		// Scan the structure to determine if pseudoknot(s) are present
		bool nopseudos = false;
		int init = 1;
		ct->energy[count] = 0;
		while (!nopseudos) {
			ct->energy[count] += pseudoenergy(data, ct, count, &nopseudos, ctbroken, &init);
		}

		// Original efn2 takes over from here
		ct->energy[count] += ergexterior(count, ct, data);

#if defined(debugmode)
		gcvt((float (ct->energy[count]))/conversionfactor,6,temp);
		out << "Exterior loop = "<<temp<<"\n";	
#endif
		
		i=1;
		while (i<ct->numofbases) { 	
			if (ct->basepr[count][i]!=0) {
				push(&stack,i,ct->basepr[count][i],1,0);
				i = ct->basepr[count][i];
			}
			i++;
		}

subroutine://loop starts here (I chose the goto statement to make the code more readable)
		pull(&stack,&i,&j,&open,&null,&stz);//take a substructure off the stack

		while (stz!=1){
			while (ct->basepr[count][i]==j) { //are i and j paired?

#if defined(debugmode) 
				tempsum=0;
#endif

				while (ct->basepr[count][i+1]==j-1) {//are i,j and i+1,j-1 stacked?
					ct->energy[count] = ct->energy[count] + erg1(i,j,i+1,j-1,ct,data);

#if defined(debugmode)
					tempsum = tempsum + erg1(i,j,i+1,j-1,ct,data); 
					gcvt((float (erg1(i,j,i+1,j-1,ct,data)))/conversionfactor,6,temp);
					out << "\tStack = "<<temp<<"  for stack of "<<i<<"-"<<j<<"\n";
#endif

					i++;
					j--;
				}

#if defined(debugmode)
				gcvt((float (tempsum))/conversionfactor,6,temp);
				out << "Helix total = "<<temp<<"\n";	
#endif

				sum = 0;
				k = i + 1;

				// now efn2 is past the paired region, so define the intervening non-paired segment
				while (k<j) {
					if (ct->basepr[count][k]>k)	{
						sum++;
						ip = k;
						k = ct->basepr[count][k] + 1;
						jp = k-1;
					}
					else if (ct->basepr[count][k]==0) k++;
				}

				if (sum==0) {//hairpin loop
					ct->energy[count] = ct->energy[count] + erg3(i,j,ct,data,fce.f(i,j-i));

#if defined(debugmode)
					gcvt((float (erg3(i,j,ct,data,fce.f(i,j-i))))/conversionfactor,6,temp);
					out << "Hairpin = "<<temp<<"  for closure of "<<i<<"-"<<j<<"\n";	
#endif
					
					goto subroutine;
				}

				else if (sum==1) {//bulge/internal loop
					ct->energy[count] = ct->energy[count] +	erg2(i,j,ip,jp,ct,data,fce.f(i,ip-i),fce.f(jp,j-jp));

#if defined(debugmode)
					gcvt((float (erg2(i,j,ip,jp,ct,data,fce.f(i,ip-i),fce.f(jp,j-jp))))/conversionfactor,6,temp);
					out << "Internal/bulge = "<<temp<<"  for closure of "<<i<<"-"<<j<<"\n";	
#endif

					i = ip;
					j = jp;
				}

				else {//multi-branch loop
					ct->energy[count] = ct->energy[count]+ ergmulti(count, i, ct, data, simplemb);
					
#if defined(debugmode)
					gcvt((float (ergmulti(count, i, ct, data)))/conversionfactor,6,temp);
					out << "Multiloop = "<<temp<<"  for closure of "<<i<<"-"<<j<<"\n";	
#endif

					//put the exiting helices on the stack:
					sum++;//total helixes = sum + 1
					i++;
					for (k=1;k<sum;k++) {
						while (ct->basepr[count][i]==0) i++;
						push (&stack,i,ct->basepr[count][i],1,0);
						i = ct->basepr[count][i]+1;
					}
					
					goto subroutine;
				}
			}    
		}
		
#if defined(debugmode)
		gcvt((float (ct->energy[count]))/conversionfactor,6,temp);
		out << "\n\nTotal energy = "<<temp<<"\n";
		out.close();
#endif

		return;
}
}

// Scan structure ct[count] for pseudoknots, returning the energy for breaking pseudoknot stem(s) * 10
int pseudoenergy(datatable *data, structure *ct, int count, bool *nopseudo, structure *ctbroken, int *pinit) {
	int i, iter, a=0, b=0, c=0, d=0, e=0, f=0, j=0;
	int stem1bps = 0, stem2bps = 0;
	int init = *(pinit);
	bool rightSingle, findCF, leftSingle, findAD;
	for(init; init <= ct->numofbases; init++) {
		if(ct->basepr[count][init] > init) {
			for(i = 1; i < init; i++) {
				if(ct->basepr[count][i] > init && ct->basepr[count][i] < ct->basepr[count][init]) {
					// Update pointer so that we remember where the pseudoknot search left off
					*(pinit) = init;
					// Define stems of pseudoknotted region
					j = ct->basepr[count][init];
					e = ct->basepr[count][i];
					b = ct->basepr[count][j];

					rightSingle = false;
					findCF = false;
					f = e+1; // Initial guess
					while( findCF == false && rightSingle == false) {
						if( ct->basepr[count][f] != 0 && ct->basepr[count][f] < e) {
							c = ct->basepr[count][f];
							findCF = true;
						}
						else {
							f++;
							if( f >= j) {
								rightSingle = true;
								f = j;
								c= b;
							}
						}
					}

					leftSingle = false;
					findAD = false;
					a = b - 1; // Initial guess
					while( findAD == false && leftSingle == false) {
						if( ct->basepr[count][a] != 0 && ct->basepr[count][a] >  c) {
							d = ct->basepr[count][a];			
							findAD = true;
						}
						else {
							a--;
							if( a <= i) {
								leftSingle = true;
								a = i;
								d = e;
							}
						}
					}
#if defined(VERY_VERBOSE_MODE)
					cout << "Pseudoknot detected in structure " << count << ": " << 
						i << "-" << e << "," << a << "-" << d << " ; " <<
						c << "-" << f << "," << b << "-" << j << "\n";
#endif

					// Reset the current ctbroken basepr array for storing pairing information of the latest 
					// broken pseudoknot stem
					for (j = 1; j <=ctbroken->numofbases; j++){
						ctbroken->basepr[count][j] = 0;
					}

					// Determine which pseudoknot formation penalty to use
					int beta1 = 0;
					bool pseudoOrMulti = false;
					for (iter = 1; iter <= i; iter++) {
						if (ct->basepr[count][iter] > j) pseudoOrMulti = true;
					}
					if (pseudoOrMulti) beta1 = BETA_1PM;
					else beta1 = BETA_1;

					// Calculate RNAstructure hairpin energy of each pseudoknot stem
					ctbroken->basepr[count][a] = d;
					ctbroken->basepr[count][d] = a;
					efn2(data, ctbroken, count, false);
					int hairpin1Energy = ctbroken->energy[count];
					int stackingEnergy1 = data->tstack[ctbroken->numseq[d]][ctbroken->numseq[a]]
					[ctbroken->numseq[d+1]][ctbroken->numseq[a-1]];
					ctbroken->basepr[count][a] = 0;
					ctbroken->basepr[count][d] = 0;

					ctbroken->basepr[count][c] = f;
					ctbroken->basepr[count][f] = c;
					efn2(data, ctbroken, count, false);
					int hairpin2Energy = ctbroken->energy[count];					
					int stackingEnergy2 = data->tstack[ctbroken->numseq[f]][ctbroken->numseq[c]]
					[ctbroken->numseq[f+1]][ctbroken->numseq[c-1]];
					ctbroken->basepr[count][c] = 0;
					ctbroken->basepr[count][f] = 0;

#if defined(VERY_VERBOSE_MODE)
					cout<<"Hairpin Energies: "<<hairpin1Energy<<" "<<hairpin2Energy<<"\n";
#endif

					// Calculate paired and unpaired nucleotide penalties for the interior of the pseudoknot
					int beta2 = 0, beta3 = 0;

					// Add 2 beta2 penalties for a-d and c-f pairs
					beta2 += 2*(BETA_2);

					// Walk through gap region 1
					iter = a+1;
					while (iter < b){
						// Check for closing bp of a contained secondary structure in the gap
						if (ct->basepr[count][iter] > iter && ct->basepr[count][iter] < b){
							beta2 += BETA_2;
							iter = ct->basepr[count][iter] + 1;
						}
						else if (ct->basepr[count][iter] != 0 && 
							(ct->basepr[count][iter] <= a || ct->basepr[count][iter] >= b)) {
							// If nt(iter) pairs outside the gap, add large penalty to the structure
							beta2 += 1000;
							iter++;
						}
						else {
							// Handles cases where nt(iter) is unpaired
							beta3 += BETA_3;
							iter++;
						}
					}

					// Walk through gap region 2
					iter = c+1;
					while (iter < d){
						// Check for closing bp of a contained secondary structure in the gap
						if (ct->basepr[count][iter] > iter && ct->basepr[count][iter] < d){
							beta2 += BETA_2;
							iter = ct->basepr[count][iter] + 1;
						}
						else if (ct->basepr[count][iter] != 0 && 
							(ct->basepr[count][iter] <= c || ct->basepr[count][iter] >= d)) {
							// If nt(iter) pairs outside the gap, add large penalty to the structure
							beta2 += 1000;
							iter++;
						}
						else {
							// Handles cases where nt(iter) is unpaired
							beta3 += BETA_3;
							iter++;
						}
					}

					// Walk through gap region 3
					iter = e+1;
					while (iter < f){
						// Check for closing bp of a contained secondary structure in the gap
						if (ct->basepr[count][iter] > iter && ct->basepr[count][iter] < f){
							beta2 += BETA_2;
							iter = ct->basepr[count][iter] + 1;
						}
						else if (ct->basepr[count][iter] != 0 && 
							(ct->basepr[count][iter] <= e || ct->basepr[count][iter] >= f)) {
							// If nt(iter) pairs outside the gap, add large penalty to the structure
							beta2 += 1000;
							iter++;
						}
						else {
							// Handles cases where nt(iter) is unpaired or paired to something outside the gap
							beta3 += BETA_3;
							iter++;
						}
					}

					// Scans both stems to determine which one is simpler to break
					for (iter = i; iter <= a; iter++) {
						if (ct->basepr[count][iter] > iter) stem1bps++;
					}
					for (iter = d; iter <= e; iter++) {
						if (ct->basepr[count][iter] > iter) stem1bps++;
					}
					for (iter = b; iter <= c; iter++) {
						if (ct->basepr[count][iter] > iter) stem2bps++;
					}
					for (iter = f; iter <= j; iter++) {
						if (ct->basepr[count][iter] > iter) stem2bps++;
					}

					if (stem1bps <= stem2bps) {
#if defined(VERY_VERBOSE_MODE)
						cout << "Breaking stem 1...";
#endif

						for (iter = i; iter <= a; iter++) {
							if (ct->basepr[count][iter] >= d && ct->basepr[count][iter] <= e) {
#if defined(VERY_VERBOSE_MODE)
								cout << iter << "-" << ct->basepr[count][iter] << " ";
#endif

								ctbroken->basepr[count][iter] = ct->basepr[count][iter];
								ctbroken->basepr[count][ct->basepr[count][iter]] = 
									ct->basepr[count][ct->basepr[count][iter]];
								ct->basepr[count][ct->basepr[count][iter]] = 0;
								ct->basepr[count][iter] = 0;
							}
						}
#if defined(VERY_VERBOSE_MODE)
						cout << "\n";
#endif
					}

					else {
#if defined(VERY_VERBOSE_MODE)
						cout << "Breaking stem 2...\n";
#endif

						for (iter = b; iter <= c; iter++) {
							if (ct->basepr[count][iter] >= f && ct->basepr[count][iter] <= j) {
#if defined(VERY_VERBOSE_MODE)
								cout << iter << "-" << ct->basepr[count][iter] << " ";
#endif

								ctbroken->basepr[count][iter] = ct->basepr[count][iter];
								ctbroken->basepr[count][ct->basepr[count][iter]] = ct->basepr[count][ct->basepr[count][iter]];
								ct->basepr[count][ct->basepr[count][iter]] = 0;
								ct->basepr[count][iter] = 0;
							}
						}
#if defined(VERY_VERBOSE_MODE)
						cout << "\n";
#endif
					}
					
					// **Calculate energy of broken stem
					efn2(data, ctbroken, count, false);

#if defined(VERY_VERBOSE_MODE)
					cout << "Energy from broken stem = " << ctbroken->energy[count] << "\n";
#endif

					int uncorrectedEnergy = beta1+beta2+beta3+ctbroken->energy[count];

					// Update zero-order ctbroken basepr array with the latest broken pseudoknot stem
					for (j = 1; j <=ctbroken->numofbases; j++){
						if (ctbroken->basepr[count][j]!=0) ctbroken->basepr[0][j] = ctbroken->basepr[count][j];
					}

					return (uncorrectedEnergy-(hairpin1Energy-stackingEnergy1)-(hairpin2Energy-stackingEnergy2));
				}
			}				
		}
	}
	// no more pseudoknots found
	*(nopseudo) = true;
	return 0;
}

int ctout2 (structure *ct,int cutoff,const char *ctoutfile) {
	int count,i,j=0;
	char line[2*ctheaderlength],number[2*numlen];

	FILE *ctfile;
	ctfile=fopen(ctoutfile,"w");

	for (count=1;count<=(ct->numofstructures);count++) {

		if (cutoff !=0 && ct->energy[count] > cutoff) break;
		j++;

		strcpy(line,"");
		sprintf(line,"%5i",ct->numofbases);

		
		if (ct->energy[count]!=0) {
   			strcat(line,"  ENERGY = ");

			//gcvt((float (ct->energy[count]))/conversionfactor,6,number);
			if (conversionfactor==10)
				sprintf(number,"%.1f",(float (ct->energy[count]))/conversionfactor);
			else if (conversionfactor==100)
				sprintf(number,"%.2f",(float (ct->energy[count]))/conversionfactor);
			else sprintf(number,"%f",(float (ct->energy[count]))/conversionfactor);
	
   			strcat(line,number);
   			strcat(line,"  ");
		}
		else strcat(line,"  ");
		strcat(line,ct->ctlabel[count]);
		fputs (line,ctfile);
		for (i=1;i<ct->numofbases;i++) {
			if (ct->stacking) {
				sprintf(line,"%5i%2c%8i%5i%5i%5i%5i\n",
					i,ct->nucs[i],(i-1),(i+1),ct->basepr[count][i],ct->hnumber[i],ct->basepr[count][i+ct->numofbases]);
			}
			else {
				sprintf(line,"%5i%2c%8i%5i%5i%5i\n",
					i,ct->nucs[i],(i-1),(i+1),ct->basepr[count][i],ct->hnumber[i]);
			}
			fputs(line,ctfile);
		}
		
   
		//last nucleotide not connected--
		i = ct->numofbases;
		if (ct->stacking) {
			sprintf(line,"%5i %1c%8i%5i%5i%5i%5i\n",
				i,ct->nucs[i],(i-1),0,ct->basepr[count][i],ct->hnumber[i],ct->basepr[count][i+ct->numofbases]);
		}
		else {
			sprintf(line,"%5i %1c%8i%5i%5i%5i\n",
				i,ct->nucs[i],(i-1),0,ct->basepr[count][i],ct->hnumber[i]);
		}
		fputs(line,ctfile);

	}

	fclose (ctfile);
	return(j);
}

void pairout(structure *ct, int cutoff, const char* pairsoutfile) {
	int iter, iter2, a, b, c, d, bpCount;
	ofstream outFile(pairsoutfile);
	if (outFile.bad()) {
		cerr << "Error opening file " << pairsoutfile << "\n";
		return; 
	}
	for (iter = 1; iter <= ct->numofstructures; iter++) {
		if (cutoff !=0 && ct->energy[iter] > cutoff) break;
		outFile << "#  Structure #" << iter << "; ENERGY = " << double(ct->energy[iter])/10 << "; " << 
			ct->ctlabel[iter] << endl;
		outFile << "1" << endl;
		for (iter2 = 1; iter2 <= ct->numofbases; iter2++) outFile << ct->nucs[iter2];
		outFile << endl << endl;
		outFile << "Positions paired." << endl;
		iter2 = 1;
		while (iter2 <= ct->numofbases) {
			if (ct->basepr[iter][iter2] < iter2) {
				iter2++;
				continue;
			}
			else {
				a = iter2;
				b = ct->basepr[iter][iter2];
				bpCount = 1;
				while (ct->basepr[iter][iter2+bpCount] == (b-bpCount)) bpCount++;
				c = iter2 + bpCount - 1;
				d = ct->basepr[iter][iter2+(bpCount-1)];
				outFile << a << "-" << c << "; " << d << "-" << b << endl;
				iter2 = c + 1;
				continue;
			}
		}
		outFile << endl << "---------------------------------------" << endl << endl;
	}
	if (outFile.bad()) {
		cerr << "Error writing to file " << pairsoutfile << "\n";
		return; 
	}
	outFile.close();
	return;
}
