#!/usr/bin/python

from sys import argv
from os import system, chdir
from os.path import exists,basename
import string

seq_files = argv[1:]
EXE = './run_partition.py'
bsub_file = 'bsubMINI'

fid = open( bsub_file,'w')

tot_jobs = 0

seq_files_to_do = []

for i in range( len( seq_files ) ):
    seq_file = seq_files[ i ]
    bpp_file = seq_file.replace( '.seq','.bpp' )

    if not exists( bpp_file ):

        seq_files_to_do.append( seq_file )


outfile = '/dev/null'
errfile = '/dev/null'

NCHUNK = 10

for i in range( int( len( seq_files_to_do )/NCHUNK )+1 ):
    command =  'bsub -W 16:0 -o %s -e %s %s ' % (outfile, errfile, EXE )

    found_job = 0
    for q in range( NCHUNK*i, NCHUNK*(i+1) ):
        if q < len( seq_files_to_do ):
            command += ' ' + seq_files_to_do[q]
            found_job = 1

    if found_job: fid.write( command + '\n' )

    tot_jobs += 1

print 'Created bsub submission file ',bsub_file,' with ',tot_jobs, ' jobs queued. To run, type: '
print '>source',bsub_file
print
