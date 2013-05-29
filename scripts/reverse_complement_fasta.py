#!/usr/bin/python

from sys import argv
from os import system, chdir
from os.path import exists
import string

input_fasta_file = argv[1]
lines = map( lambda x:x[:-1], open( input_fasta_file ).readlines() )

# parse into a bunch of sequences -- should make
# use of a universal FASTA reader, but anyway...
tags = []
sequences = []
sequence = ''
for line in lines:
    if len( line ) > 1 and line[0] == '>':
        tags.append( line[1:] )
        if len( tags ) > 1: sequences.append( sequence )
        sequence = ''
        continue
    sequence += line
sequences.append( sequence )

outfile = input_fasta_file.replace( '.fasta','_REVERSE_COMPLEMENT.fasta' )
print outfile

RC_map = { 'A':'U', 'U':'A', 'G':'C', 'C':'G', '.':'.'};
def RC( seq ):
    seq_out = ''
    seq_rev = seq[::-1]
    for m in seq_rev:
        seq_out += RC_map[ m ]
    return seq_out

print 'Creating outfile: ', outfile
fid_fasta =  open( outfile, 'w' )
for i in range( len(sequences) ):

    fid_fasta.write( '>%s\n' % tags[i] )
    seq = sequences[i].replace('.','').replace( 'T','U' )
    seq = RC( seq )
    fid_fasta.write( '%s\n\n' % seq )

fid_fasta.close()

