########################
    BPPalign  v0.5
########################
   (c) Rhiju Das, Stanford University 2011
   rhiju@stanford.edu

Used in publication: Kladwang, W, Chou, F.C., Das, R., "Automated RNA structure prediction uncovers a missing link in double glycine riboswitches", under review. (2011).

This tool provides basic scripts to run the partition executable RNAstructure on hundreds of homolog sequences from RFAM, and to then visualize the results by averaging across these matrices. Examples from purine-binding and double glycine-binding riboswitches are included.

Requirements:  
1. Python 2.x or later.
2. MATLAB
3. This SVN archive includes a version of the RNAstructure executable 'partition' that you need to compile for UNIX machines. This executable has been modiifed to output the base pair probability matrix into a text file bpp.txt.
4. Internet connection (for Entrez queries to get nucleotide sequences)

BPPalign has been tested on Mac and some Linux variants.

#################################
EXAMPLE 1   [purine riboswitch]
#################################

Go into the directory add/.

An example starting alignment, taken from the RFAM database, is there:

  add_rfam.txt 

First, we want to get sequence files that correspond these sequences. One of the main uses of bpp_align is to look for structure in flanking sequences that may have been missed in the RFAM alignment. So we retrieve the original sequences along with 5' and 3' sequences, through queries to Entrez Pubmed:

../scripts/get_original_sequences.py  add_rfam.txt add_seqs 10 70

The numbers 10 and 70 mean look 10 nucleotides upstream (5') of the RFAM alignment boundary and 70 nucleotides downstream (3') of the boundary. This command should retrieve the sequences from the Entrez nucleotide database into the folder add_seqs/  as .fasta files and also as .seq files appropriate for input to the RNAstructure executable for predicting secondary structure.  There is an additional file 'new_align.txt' that includes alignment information for these sequences.

Then, let's run RNAstructure.

../scripts/run_partition.py  add_seqs/*seq 

This requires having a working executable in RNAstructure/exe/partition -- please check directions for compiling RNAstructure. Note that the version included with BPPalign has a modification to automatically output the base pairing probability matrix bpp to a text file.

 To visualize the data go into the scripts directory, and start MATLAB. Run the commands in 

   plot_bpp_ADD_script.m 

Check out that script for comments on what input is needed to get the plot. The plot should look like BPPalign_add_example.eps


##########################################
EXAMPLE 2   [double glycine riboswitch]
##########################################
Go into directory gly/. Run 

  get_double_sequences.py

This looks at the alignment of single-riboswitch domains in gly_rfam.txt and finds pairs of single aptamers that can go together into double riboswitches. Some nucleotides in the linker are tagged 'X'. The output is in gly2_rfam.txt. 

Then, again. retrieve sequences with downstream and upstream flanking regions -- also fillin 'X' nucleotides in the linker:

../scripts/get_original_sequences.py  gly2_rfam.txt gly_seqs 100 100

Here, 100 100 means get 100 nucleotides both upstream and downstream of the double-riboswitch domain. Then run RNAstructure on all these sequences:

../scripts/get_original_sequences.py  gly2_rfam.txt gly_seqs 100 100

Go back up and into the scripts/ directory. In MATLAB, run the commands in 

   plot_bpp_GLY_script.m 

The default commands in this script select out sequences 1-360 (the rest are marine metagenome sequences and appear pretty redundant), and align them to the F. nucleatum sequence.

