outpath = '../gly/EXAMPLE_OUTPUT/gly_seqs/';
[nres, all_bpp ] = read_bpps( outpath );

% sequence of 'reference' sequence -- here the adenine riboswitch add from Vibrio vulnificus
% This is actually not necessary below, just put here as a check.
ref = 'AE009951.2'; % F. nucleatum sequence.
sequence =  'AAUGUCAAUAAAUAAAAUUUAUGUUAUCUAAAUUUUUCUUAUUGACAAAAAUAUAAAAAAGUGAUAAUAUUUCCAUCAUAAAAACUGAAUAAAUAAUCGGAUGAAGAUAUGAGGAGAGAUUUCAUUUUAAUGAAACACCGAAGAAGUAAAUCUUUCAGGUAAAAAGGACUCAUAUUGGACGAACCUCUGGAGAGCUUAUCUAAGAGAUAACACCGAAGGAGCAAAGCUAAUUUUAGCCUAAACUCUCAGGUAAAAGGACGGAGUAAUUGUGCAAUUUAUAUUUUUAAUAUAUUUUAAAAUAUAUUUUUUAUUUGUAUAAUUCUUUUUUAAUUUCAGAGGUCUUAAUAAAGACCUUUUUUUAUUCAGUAAAA';
structure = '.........................................................................................................((((((((......((((((....)))))).(((...((((.....))))..)))........))))))))........(((((......(((((.....))))).(((...((((....((((....)))).....))))..))).......)))))............................................................................................................';
%structure = '..(((((((((..((((((((.......))))))))..))))))))).............(((((.........)))))....((((((((((...(((......((((((((......((((((....)))))).((((((((......)))))..)))........))))))))...)))..(((....))).(((((.....))))).(((...((((....((((....)))).....))))..))).((((((......((((((((((...(((.....((((((((....))))))))...))))))))))))).))))))......(((((((((....))))))))).))))))))))....';


%ref = 'AE003852.1' % V. cholera
%sequence =  'GCGCAACAUUUCAUCAGUUAAUUUCUUUUGUAGGUGACAUCACAUUUUCUUUUCGCUAGUAUCCGCCUUGCAAAUCGUUUUAUUCAAGACGAUUGUUCCGUUGAAGACUGCAGGAGAGUGGUUGUUAACCAGAUUUUAACAUCUGAGCCAAAUAACCCGCCGAAGAAGUAAAUCUUUCAGGUGCAUUAUUCUUAGCCAUAUAUUGGCAACGAAUAAGCGAGGACUGUAGUUGGAGGAACCUCUGGAGAGAACCGUUUAAUCGGUCGCCGAAGGAGCAAGCUCUGCGCAUAUGCAGAGUGAAACUCUCAGGCAAAAGGACAGAGGAGUGAAAGGCCAAUCUUUUAGUGAGCUCGCUAGAGCUCUGCUGUGCAUUUUUCGCACCCUUUCUCUCUUCUCCUUAUGUUUAUUUUUUAAGGGGAAAUCAUGAAUAA'
%structure = '.........................................................................................................((((((((......(((((((...(((((......))))).....))))))).(((...((((.....))))..)))((.((((((...((((.....))))...)))))))).....)))))))).......((((((......((((......)))).(((...((((...((((((((....))))))))....))))..))).......))))))...........................................................................................................';


%ref = 'CP000431.1';
%sequence = 'UUAACACGCCGGAACGGCGAAGUCAUGCCUUCGUAAUACGUCUGCGGAACAGUACGUCCCGAGUCAGUUCGAUUCGGAACCAUCACGGUGGUUUCUCCCUAUGAACCUGGCUGGAGAGUUCCGGGACCGAGUCCCGGACGCCGAAGGAGCAAGUACCCACAAACUCUCAGGCACCAGGACAGUCAGGGUCGGGAACUCUGGAGAGCAGCGUGUCGGCUGUAACGGCACGUUCACCGACGGGGCACGAUGACCACUCCCACGAGUGUGUCUCCAAACUCUCAGGUUACGGAACAGAGCGGGUGUAGCCGAGGCGAUCCACCGGAUCACCCUCACAUCGCUCGAGCGAUCGGUUCGCCGGUUACUCACGCUCCCUGGGAGAUCCAUGAAUGUGACUGCUAGUGG'
%structure = '';

%whichseq = [];
%box_bounds = [+100 100];
%bpp_mean = make_bpp_plot( nres, all_bpp, outpath, ref, structure, box_bounds, [], 0.85 );

% no metagenome
whichseq = [1:360];
box_bounds = [+100 100];
bpp_mean = make_bpp_plot( nres, all_bpp, outpath, ref, structure, box_bounds, whichseq );

%whichseq = [35 361:1453]; % just metagenome
%box_bounds = [+100 100];
%bpp_mean = make_bpp_plot( nres, all_bpp, outpath, ref, structure, box_bounds, whichseq );


%show_AUG_and_polyU_script( outpath, ref )







