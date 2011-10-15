% sequence of 'reference' sequence -- here the adenine riboswitch add from Vibrio vulnificus
% This is actually not necessary below, just put here as a check.
%sequence  = 'ACAGGUGCUAGACUUUCGGCGAUCAACGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUCUGUCGCUUUAUCCGAAAUUUUAUAAAGAGAAGACUCAUGAAUUACUUUGACCUGCCGAAGAUCGAUUUGCACUGCCACCUAGAUGGUAGCGUACGCCCGCAGACCAUUAUUGACCUUGCUGACGAACAGAACCUCACUUUGCCAUCACGUGACAUUAACGUGA'

% If we input structures in dot-bracket notation we can get them plotted.
% 'native structure' -- includes a putative translation repressor stem
structure1= '......................................((((((.........))))))........((((((.......)))))).((((((.((((.((.((((((..........)))))).))))))))))))...........................................................................................................................'; 
% additional reference structure -- derived from crystallographic (ligand-bound) model 1y26
structure2= '..........................(((((((((...((((((.........))))))........((((((.......))))))..)))))))))...................................................................................................................................................................';
structures = { structure1, structure2 };

% read in the sequences. They will be aligned based on the new_align.txt file which should be in outpath. 
outpath = '../add/add_seqs/';  ref = 'AE016796.1'; % Genbank ID for the reference sequence.
[nres, all_bpp ] = read_bpps( outpath );

% how many residues upstream and downstream compared to 'conventional' boundary -- i.e., input parameters to get_original_sequences.py
box_bounds = [10 70 ]; 
bpp_mean = make_bpp_plot( nres, all_bpp, outpath, ref, structures, box_bounds );
