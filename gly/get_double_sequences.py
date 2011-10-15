#!/usr/bin/python

from sys import argv
from os import system, chdir
from os.path import exists
import string

rfam_file = argv[1]

lines = map( lambda x:x[:-1], open( rfam_file ).readlines() )

ids = []
seq_starts = []
seq_stops = []
seqs = []
for i in range( len(lines) ):
    line = lines[i]

    cols = string.split( line )
    tag = cols[0]
    seq = cols[1]
    seqs.append(  seq )

    cols = tag.split( '/' )
    id = cols[0]
    cols = cols[1].split('-')
    seq_start = int( cols[0] )
    seq_stop = int( cols[1] )

    ids.append( id )
    seq_starts.append( seq_start )
    seq_stops.append( seq_stop )


outfile = 'gly2_rfam.txt'
fid = open( outfile, 'w' )

match_ids = []
for i in range( len( ids ) ):
    id1 = ids[ i ]
    for j in range( i+1, len( ids ) ):
        id2 = ids[ j ]
        if ( id1 == id2 ):

            pair = [i,j]
            if seq_starts[ i ] < seq_stops[ i ]: # 5' to 3'
                if seq_stops[i] > seq_starts[j]:
                    pair = [j,i]
            else: # 3' to 5'
                if seq_stops[i] < seq_starts[j]:
                    pair = [j,i]

            seq_sep = abs( seq_stops[pair[0]] - seq_starts[pair[1]] )

            print 'MATCH! ', id1, i, j, seq_sep

            MAX_SEQ_SEP = 9
            if seq_sep > MAX_SEQ_SEP: continue

            seq_sep_string = ''
            for m in range( seq_sep-1 ): seq_sep_string += 'X' # unknown nucleotide.
            for m in range( seq_sep, MAX_SEQ_SEP ): seq_sep_string += '.' # gap

            print seq_sep, seq_sep_string

            id_tag = '%s/%d-%d' % (id1,seq_starts[ pair[0] ],  seq_stops[ pair[1] ])

            for i in range( len(id_tag), 40 ): id_tag += ' '

            fid.write( '%s %s%s%s\n' % ( id_tag,\
                                               seqs[ pair[0] ], seq_sep_string, seqs[ pair[1] ] ) )


fid.close()
print len( match_ids )
print 'Made: ', outfile

