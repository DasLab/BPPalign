#!/usr/bin/python

from sys import argv
from os import system, chdir,getcwd
from os.path import exists,basename, abspath
import string

seq_files = argv[1:]

for seq_file in seq_files:
    EXE = '~/projects/rdat/external/RNAstructure/exe/partition'

    bpp_file = abspath( seq_file ).replace( '.seq','.bpp' )

    ORIGINALDIR = getcwd()

    workdir =  seq_file.replace( '.seq', 'TEMP' )
    system( 'mkdir -p '+workdir )
    system( 'cp '+seq_file+' '+workdir )
    chdir( workdir )

    command = '%s %s blah.ct' % (EXE, basename(seq_file) )
    system( command )

    command = 'mv bpp.txt '+bpp_file
    system( command )

    chdir( ORIGINALDIR )
    system( 'rm -rf '+workdir )
