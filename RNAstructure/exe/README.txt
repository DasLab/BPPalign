This is a set of text interfaces for RNAstructure.

The shared libraries (.so or .dynlib) are required 
for function.  So are the thermodynamic parameters
in the data_tables/ directory.

For proper function, set an environment variable
called DATAPATH to indicate the location of 
data_tables.  For example, in BASH:
export DATAPATH=/home/dhm/RNAstructure/data_tables.

Note the format of input sequences to 
RNAstructure.  There are examples in the
examples/ directory.  Also, FASTA formatted 
sequences are additionally allowed for the text
interfaces.  A sample FASTA format sequence is also
in the examples/ directory.

Most programs provide brief usage information
if invoked with no input and extensive
help if invoked with a -h parameter.

Please report any bug to David Mathews:
David_Mathews@urmc.rochester.edu.


