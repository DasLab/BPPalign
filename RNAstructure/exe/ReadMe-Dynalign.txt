Dynalign Help.

This readme accompanies Dynalign, an algorithm for
simultaneously predicting the lowest free energy RNA
secondary structure common to two sequences and the
alignment of the sequences.  It was described in
detail in: Mathews & Turner, Journal of Molecular
Biology, 317:191-203 (2002) and Mathews,
Bioinformatics, 21:2246-2253 (2005).

This version of Dynalign includes suboptimal
structure and alignment prediction.

This version also uses the accelerations described
in Uzilov, Keegan, & Mathews, BMC Bioinformatics, 
7:173 (2006) and those described in Harmanci, 
Sharma, & Mathews, BMC Bioinformatics, 8:130 (2007).

This version is capable of being run on a
multi-processor (SMP) machine for faster processing,
thanks to the work of Chris Connett, Andrew Yohn,
and Paul Tymann of the Rochester Institute of 
Technology.  A POSIX-compliant threading library 
(e.g. Pthreads) is required by this feature.



Dynalign can be run from the command line using a
configuration file, or in an interactive mode:

dynalign [configfile]

If a single command-line argument is given, it is
read as the configuration filename.  See
ReadMe_configuration.txt for more information on the
configuration file format.

With no command-line arguments, Dynalign will
interactively prompt for the required information.

Here is a summary of required input:
inseq1 and inseq22 are the two input sequence
files. 

outct1 and outct2 are the output files for structures and
are in the connect table format.

alignment is the output file for the alignments.

The following are the optional inputs:
M is the maximum separation parameter.  This is 
largely deprecated.  We advise entering -99 so that
the alignment constraint as described by Harmanci
et al. is used.  
This version of dynalign uses the following scheme 
for M (if M is used):  For 
nucleotide i from sequence 1 to align to nucleotide
k in sequence 2:
| i * N2/N1 - k | <= M
where N1 is the length of sequence 1 and N2 is the
length of sequence 2.  M = 6 works for 5S rRNA and
tRNA.  A larger M may be required for
longer sequences.  

The gap_cost (fgap) is used to discourage the introduction
of gaps in the alignemnt.  0.4 kcal/mol/gap is the
recommended value.

max_percent_diff is the maximum percent difference
in free energy in the suboptimal structures.  20
(for 20%) is a good starting value.

bp_window specifies how different the suboptimal
structures must be from each other.  2 is probably a
good starting point and smaller values will result
in more suboptimal structures.

Similarly, the align_window specifies how different
the suboptimal alignments must be from each other.
1 is a good starting point.

single_bp_inserts specifies whether single base pair
inserts are allowed; 0 = no; 1 = yes.  A savefile
can be specified.  These are required to generate
dot plot information.

Files can be specified that contain structure and/
or alignment constraints.  These are described in
ReadMe_constraint.txt.  0=no constraints; 1=read
constraint file.

Local alignments can be performed by using a 
configuration file and setting local = 1.  Note that
local alignment only applies to the calculation by
Dynalign.  The alignment pre-filter always runs 
in global mode. By default, Dynalign runs in global
mode.

In the SMP version of Dynalign, num_processors
should be set to the total number of processors (or
processor cores on multi-core CPUs) available on the
system.

Dot plots are described in ReadMe_dotplot.txt.

Dynalign is also available as a user-friendly GUI
for Microsoft Windows in RNAstructure.  It can be
downloaded at http://rna.urmc.rochester.edu .
RNAstructure also works using WINE on Linux.

Note that on Macintosh OS X, there is a small
default stack limit.  To run Dynalign, the stack
limit needs to be increased.  On the default shell,
bash, use:
ulimit 4096
or, if you are using tcsh, use:
limit stack 4096

If you have any trouble with Dynalign, please 
email David Mathews:
David_Mathews@urmc.rochester.edu
