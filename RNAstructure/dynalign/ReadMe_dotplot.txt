
To generate data for dot plots, first
run Dynalign and generate a save file.

Then run the dotplot utility, 
dynalign_dotplot.  The command line
parameters are:

dynalign_dotplot input.dsv max_%_difference output.out

where input.dsv is a dynalign savefile,
max_%_difference is the maximum percent 
difference in free energy from the
lowest free energy structure for which
dots will be recorded (for example, use 
20 for 20%).  output.out is the name 
of the file that will be written with 
the dot plot data.

The output format is:
i j k l energy*20, with tab delimiters.
i and j are the nucleotides paired
in sequence 1, k and l are the
nucleotides paired in sequence 2,
and energy*10 is the lowest free
energy structure that contains these
pairs, where the free energy is in
kcal/mol and the energy is multiplied
by ten.


