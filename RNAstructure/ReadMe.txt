
This is the source code for RNAstructure, which includes tools for
RNA secondary structure prediction and analysis.

It is provided to you under the GNU GPL, which is spelled out
in gpl.txt.

----

requirements.txt spells out the current requirements for compiling
the components.

building.txt provides instruction for compiling and linking the tools
in a *nix environment.

----

A note on "make install"

If "make install" is run, an error may occur in permissions. This is because
"make install" can only be run as root. It copies the contents of the exe
directory to /usr/local/RNAstructure, so "make all" must be run first.

Note that it is not required to run "make install" to use RNAstructure; it is
simply a convenience target.

----

NOTE: Most calculations require that a set of nearest neighbor
folding parameters be read from disk.  To do this, the programs 
read the DATAPATH environment variable, which needs to be set to
the location of the data_tables directory.

For example, in BASH, this is accomplished with:
export DATAPATH=[directory in which RNAstructure resides]/RNAstructure/data_tables/

This should probably be in your login script.

----

Organization of RNAstructure project directories:

Code:
RNA_class	A class library which wraps most back end code for easy use
src		Location of back end code

GUIs:
RNAstructure_java_interface	Cross-platform Java interface
RNAstructure_windows_interface	GUI specifically for Windows

Text Interfaces:
AllSub              Generate a group of suboptimal structures
bifold              Bimolecular folding
bipartition         Bimolecular partition function for predicting pair probabilities
CircleCompare 	    Compare two structures for the same sequence, with the nucleic acid 
                        backbone arranged around a circle to facilitate easy comparisons
ct2dot              Convert ct notation to bracket notation
dot2ct              Convert bracket notation to ct notation
draw                Generate postscript images of secondary structures
dynalign            Predict the lowest free energy structure common to two sequences
DuplexFold          Bimolecular folding without intramoleculr pairs
efn2                Predict folding free energy change of a given structure
EnergyPlot          Genereate Postscript images of energy dot plots
Fold                Single sequence secondary structure prediction by free energy minimization
MaxExpect           Secondary structure prediction by maximizing expected accuracy
NAPSS               Use NMR data to improve the prediction of an RNA secondary structure
oligoscreen         Screen a set of oligonuceotides for hybridization stability and self-structure
OligoWalk           Calculate thermodynamic features of sense-antisense hybridization and predict 
                         free energy changes of oligonucleotides binding to target RNA
PARTS               Predict pair probabilities for folding of two sequences to a common structure
pfunction           Partition function calculation of base pairing probabilities
ProbabilityPlot	    Generate a base pairing probabilities dot plot from a partition function save file, 
                          exported as a Postscript file	
ProbablePair        Generate structures composed of only highly probable pairs as determined by a threshold
refold	     	    Refold a previously folded structure
RemovePseudoknots   Remove pseudoknots from a structure
scorer              Calculate sensitivity and positive predictive value for two structures being compared  
stochastic          Stochastic sampling from the Boltzman ensemble

RNAstructure_windows_text_interfaces  Interfaces specifically for Windows


Resources:
data_tables	Location of thermodynamic parameters
exe		Holds executables for text interfaces, GUIs, and any libraries
tests		Project tests for each text interface

Help:
examples	Files with known output, can be used to illustrate functions
manual		Help and some usage documentation

----
