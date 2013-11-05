% sequence of 'reference' sequence -- here the adenine riboswitch add from Vibrio vulnificus
% This is actually not necessary below, just put here as a check.
sequence   = 'GGUUUGGAUAGUUGCUAUGUCAUUUUUAAUGGGGAAUGGGCGUUGAGACGUUUUUGAAGUUAGGUUUCGAUAAUUAGCUAUCCUGUUAGCUAUACUUAAGUGUGUUACACAUUGAUUAGGUGCAGCAUAGAACUGCUGCAUGGGAAUCUGGUGAAAAUCCAGAGCUGACGCGCAGCGGUGAAGGUGCAAGUGAGUGCUUCAAUGUAGCCACUGAGAGUAUAAAAACUCUUGGGAAGGUGAGGCAAUUACUCUCGCGUAGAGCACCCCAGUCCGAAGACCGGCCUAAUCAGAAACAUGUGCUGCUGUCCUCGCGUUCUAGGUAAGCAGACCGAUAUGUAGUUGAAGUUGGU';
structure1 = '...................................................................................................................(((((((((............))).(((...(((((.......)))))..(((...))).(((...((((.........((((((......(((.((((((..........))))))...)))))))))..................))))....)))....)))))))))................................................................';

structures = { structure1 };

% read in the sequences. They will be aligned based on the new_align.txt file which should be in outpath. 
outpath = '../b12/b12_seqs/';  ref = 'BX248357.1/28472-28872';
[nres, all_bpp ] = read_bpps( outpath );

% how many residues upstream and downstream compared to 'conventional' boundary -- i.e., input parameters to get_original_sequences.py
box_bounds = [100 100 ]; 
bpp_mean = make_bpp_plot( nres, all_bpp, outpath, ref, structures, box_bounds );





      
%:::::::::::::::[[[[[[<<<____________>>>,(((,,,<<<<<_______>>>>>,,<<<___>>>,<<<---<<<<...------<<<<<<-----.<<<-<<<<<<_____...._>>>>>>--->>>>>>>>>----------........>>>>---->>>,,,,)))]]]]]]:::::::::::::::           
%uaaaauuauaccgguuauGGUcCcccaaaaagaaaggGgxAGGGAAucCGGUGaAAauCCGgaGCuGCCCCGCaACUGUAAgCGg...gaagcagcccccAaaau.gCCACUGgcccguaa....gggcCGGGAAGGcgggggcaagcgaugAc........cCGcgAGcCAGGAGACCUGCCauaaaaaauaaaauaacc 190       
%UGUGUUACACAUUGAUUAGGUGCAGCAUAGAACUGCUGCAUGGGAAUCUGGUGAAAAUCCAGAGCUGACGCGCAGCGGUGAAGGUgcaAGUGAGUGCUUCAAUGUaGCCACUGAGAGUAUAaaaaCUCUUGGGAAGGUGAGGCAAUUACUCUCGcguagagcACCCCAGUCCGAAGACCGGCCUAAUCAGAAACAUGUGCU 301       

%UGUGUUACACAUUGAUUAGGUGCAGCAUAGAACUGCUGCAUGGGAAUCUGGUGAAAAUCCAGAGCUGACGCGCAGCGGUGAAGGUGCAAGUGAGUGCUUCAAUGUAGCCACUGAGAGUAUAAAAACUCUUGGGAAGGUGAGGCAAUUACUCUCGCGUAGAGCACCCCAGUCCGAAGACCGGCCUAAUCAGAAACAUGUGCUGCUGUCCUCGCGUUCUAGGUAAGCAGACCGAUAUGUAGUUGAAGUUGGUCAUGUGAUCAGACGACUUACCACGCCCGUACUAUUUAUUGCACUUCUUAUU




 :::::::::::::::[[[[[[.<<<____________>>>,,,,,(((,,,<.<<<<_______>>>>.>,,<<<____>>>,<<<---<<<<......------<<<<<<-----<<<-<<<<<<______>>>>>>--->>>>>>>>>----------.........>>>>---->>>,,,,)))]]]]]]:::::::::::::::           
#CM            1 uaaaauuauaccgguuauGGU.cCcccaaaaagaaaggGguuAAaAGGGAAu.cCGGUGaAAauCCGg.aGCuGCCCCCGCaACUGUAAgCGg......gaagcagcccccAaaaugCCACUGgcccguaagggcCGGGAAGGcgggggcaagcgaugAc.........cCGcgAGcCAGGAGACCUGCCauaaaaaauaaaauaacc 190       
#MATCH           +AA+   A  + GGU    GU ::   + +        :: U A :GGGAAU CC:GUG+AA+UC:GG AGC:G CCCC C:ACUGUAA:CG:       +  CA:::C::      CCACUG:::    + :::CGGGAAGG ::G:::+A   A   +         :CG: AGCCAGGAGACC:GC    A+A++U+A+A+ ++            
#SEQ           7 AAAUCGCACCAAGGUGUUAGUgUUGGGUGUGUU----GAAGUGAUCGGGAAUgCCAGUGCAAUUCUGGcAGCGGACCCCGCCACUGUAACCGCaaauuaUGCCCAUUUCUA-----GCCACUGCAG----UCCUGCGGGAAGGGUAGAAAGACACACUUUuuacugggcGCGGUAGCCAGGAGACCGGCUUCAACAUUUUAUACGUAG 201       
#PP              


CUGUACAAAUCGCACCAAGGUGUUAGUGUUGGGUGUGUUGAAGUGAUCGGGAAUGCCAGUGCAAUUCUGGCAGCGGACCCCGCCACUGUAACCGCAAAUUAUGCCCAUUUCUAGCCACUGCAGUCCUGCGGGAAGGGUAGAAAGACACACUUUUUACUGGGCGCGGUAGCCAGGAGACCGGCUUCAACAUUUUAUACGUAG