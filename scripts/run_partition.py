#!/usr/bin/python

from sys import argv
from os import system, chdir
from os.path import exists,basename,dirname
import string

seq_files = argv[1:]
EXE = dirname( argv[0] ) + '/../RNAstructure/exe/partition'

for i in range( len( seq_files ) ):
    seq_file = seq_files[ i ]
    bpp_file = seq_file.replace( '.seq','.bpp' )

    if not exists( bpp_file ):
        command = '%s %s blah.ct' % (EXE, seq_file )
        print command
        system( command )

        command = 'mv bpp.txt '+bpp_file
        system( command )
