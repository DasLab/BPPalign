#!/usr/bin/python

from sys import argv
from os import system, chdir
from os.path import exists
import string

rfam_file = argv[1]
outpath = argv[2]
LOOK_BACK    = int( argv[3] )
LOOK_FORWARD = int( argv[4] )

lines = map( lambda x:x[:-1], open( rfam_file ).readlines() )

if not exists( outpath ): system( 'mkdir -p '+outpath )
print 'MOVING TO ', outpath
chdir( outpath )

align_file = 'new_align.txt'
fid_align = open( align_file, 'w' )

#for i in range( 200 ):
for i in range( len(lines) ):
    line = lines[i]

    cols = string.split( line )
    if len( cols ) == 2:
        tag = cols[0]
        seq = cols[1][:-1]
    else:
        assert( len( cols ) == 3 )
        tag = cols[1]
        seq = cols[2][:-1]

    cols = tag.split( '/' )
    id = cols[0]
    cols = cols[1].split('-')
    seq_start = int( cols[0] )
    seq_stop = int( cols[1] )

    strand_tag = '&strand=1'
    seq_start_query = seq_start - LOOK_BACK
    seq_stop_query  = seq_stop + LOOK_FORWARD

    if seq_start > seq_stop:
        strand_tag = '&strand=2'
        seq_start_query = seq_start + LOOK_BACK
        seq_stop_query  = seq_stop  - LOOK_FORWARD

    if ( seq_start_query < 1 ): seq_start_query = 1
    if ( seq_stop_query  < 1 ): seq_stop_query = 1

    #retrieve from pubmed
    fasta_file = 'seq%d.fasta' % i

    get_fasta_file = 0
    if exists( fasta_file ):
        if len( open( fasta_file ).readlines()[1] ) < 2: get_fasta_file = 1
    else:
        get_fasta_file = 1

    if get_fasta_file:
        command = 'wget "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=''%s''&rettype=fasta&seq_start=%d&seq_stop=%d%s" -O %s' % ( id, min(seq_start_query, seq_stop_query),  max(seq_start_query, seq_stop_query),  strand_tag,fasta_file )

        print command
        system( command )

    # convert to seq.
    seq_file = fasta_file.replace('.fasta','.seq')
    fasta_lines = open( fasta_file ).readlines()
    fid = open( seq_file, 'w' )
    fid.write(';\n')
    fid.write( fasta_lines[0][1:] )

    complete_sequence = ''
    for line in fasta_lines[1:-1]:
        rna_line = line.replace( 'T', 'U' )
        fid.write( rna_line )
        complete_sequence += rna_line[:-1]
    fid.write( '1\n' )
    fid.close()


    #Need to also create sequence_alignment.
    print i, abs(seq_start-seq_stop)+1, abs(seq_start_query-seq_stop_query)+1, len( complete_sequence)

    # weird special case?

    full_tag = '%s/%d-%d' % ( id, seq_start_query, seq_stop_query )
    for m in range( len(full_tag), 40 ): full_tag+=' '

    header = ''
    tail = ''
    extra_nts_at_head = seq_start - seq_start_query
    extra_nts_at_tail = seq_stop_query - seq_stop

    if ( seq_start > seq_stop ):
        extra_nts_at_head *= -1
        extra_nts_at_tail *= -1

        #print 'EXTRA: ', extra_nts_at_head, extra_nts_at_tail

    for i in range( LOOK_BACK - extra_nts_at_head ): header += '.'
    header += complete_sequence[0: extra_nts_at_head ]

    tail += complete_sequence[ -extra_nts_at_tail : -1]
    for i in range( LOOK_FORWARD - extra_nts_at_tail ): tail += '.'


    # fix X.
    seq_from_align = ''
    count = 0
    for m in range( len( seq ) ):
        if seq[m] == 'X':
            seq = seq[0:m] + complete_sequence[ extra_nts_at_head+count ] +seq[m+1:]
        if seq[m] != '.':
            count+= 1
            seq_from_align += seq[m]
    # consistency check
    print seq_from_align
    print complete_sequence[ extra_nts_at_head: -extra_nts_at_tail-1]

    fid_align.write( '%s %s%s%s\n' % (full_tag, header, seq.replace('T','U'), tail) )

fid_align.close()

print 'Outputted full align file to : %s/%s' % ( outpath, align_file)
