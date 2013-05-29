#!/usr/bin/python

from sys import argv
from os import system, chdir
from os.path import exists
import string

input_fasta_file = argv[1]
outpath = argv[2]

lines = map( lambda x:x[:-1], open( input_fasta_file ).readlines() )

if not exists( outpath ): system( 'mkdir -p '+outpath )
print 'MOVING TO ', outpath
chdir( outpath )

align_file = 'new_align.txt'
fid_align = open( align_file, 'w' )

# parse into a bunch of sequences
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


for i in range( len(sequences) ):

    # go ahead and create fasta files, to match original
    # bpp_align workflow which pulled sequences from pubmed.
    fasta_file = 'seq%d.fasta' % ( i )
    print fasta_file
    fid_fasta = open( fasta_file, 'w' )
    fid_fasta.write( '>%s\n' % tags[i] )
    seq_align = string.upper(sequences[i]).replace( 'T','U' )
    seq = seq_align.replace('.','').replace('-','')
    fid_fasta.write( '%s\n' % seq )
    fid_fasta.close()

    # convert to seq.
    seq_file = fasta_file.replace('.fasta','.seq')
    fasta_lines = open( fasta_file ).readlines()
    fid = open( seq_file, 'w' )
    fid.write(';\n')
    fid.write( fasta_lines[0][1:] )

    complete_sequence = ''
    for line in fasta_lines[1:]:
        rna_line = line.replace( 'T', 'U' )
        fid.write( rna_line )
        complete_sequence += rna_line[:-1]
    fid.write( '1\n' )
    fid.close()

    header = ''
    tail = ''
    fid_align.write( '%s %s%s%s\n' % (tags[i], header, seq_align, tail) )

fid_align.close()

print 'Outputted full align file to : %s/%s' % ( outpath, align_file)
