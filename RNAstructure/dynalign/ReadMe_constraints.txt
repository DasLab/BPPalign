Folding constraints are loaded from plain text files.

The file format for structural constraints is as 
follows:

DS:
XA
-1
SS:
XB
-1
Mod:
XC
-1
Pairs:
XD1 XD2
-1 -1
FMN:
XE
-1
Forbids:
XF1 XF2
-1 -1
 
XA are nucleotides that will be double-stranded.  
XB are nucleotides that will be single-stranded 
(unpaired).  
XC are nucleotides accessible to chemical
modification.  
XD1 and XD2 are base pairs.  
XE are nucleotides accessible to FMN cleavage 
(U's in GU pairs).  
XF1 and XF2 are base pairs that are prohibited.  
	
For multiple entries of a specific type of 
constraint, entries are listed on a separate line.  
When there is no constraint of a type, there are 
no lines required.  Note that all specifiers, 
followed by -1 or -1 -1 are expected.
Note: for all specifiers that take two arguments, it is assumed that the first
argument is the lower base pair number. 
	
	A sample file with constraints might contain:

DS:
15
25
76
-1
SS:
-1
Mod:
2
15
-1
Pairs:
16 26
-1 -1
FMN:
-1
Forbids:
15 27
-1 -1

This would indicate that nucleotides 15, 25, and 76
must be double stranded.  No nucleotides are forced
single stranded.  Nucleotides 2 and 15 are accessible
to chemical modification.  Nucleotides 16 and 26 are
base paired.  No nucleotides are FMN cleaved.  A 
base pair between nucleotides 15 and 27 is prohibited.
Note: for all specifiers that take two bases, it is 
assumed that the first base listed is the one with the 
lower number.

Alignment constraints are also loaded from plain text
files.  The format is:

S1 S2
-1 -1

where S1 is a nucleotide from sequence one that is aligned
to S2 from nucleotide 2.  For example:

10 12
11 13
-1 -1

would align nucleotides 10 and 12 and 11 and 13, from
sequences 1 and 2, respectively.
The -1 -1 is a required specifier.


