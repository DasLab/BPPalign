/*
 * (c) 2010  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.io.BufferedWriter;
import java.io.FileWriter;

/**
 * A class that writes helix structure files.
 * @author Jessica Reuter
 */
public class HelixStructureFileWriter extends StructureFileWriter {
    /**
     * Constructor
     */
    private HelixStructureFileWriter() {
	type = "Helix Files,txt";
    }

    /**
     * Prepare to export the structure to a helix file.
     * @return  a new HelixStructureFileWriter
     */
    public static HelixStructureFileWriter exportToHelixFile() {
        return new HelixStructureFileWriter();
    }

    /**
     * Determine a pair for a particular nucleotide.
     * @param counter  the nucleotide to determine pair for
     * @param structure  the structure this nucleotide is found in
     * @param strand  the strand this structure is found in
     * @return  the nucleotide the given index is paired to
     * @throws Exception  if pair determination runs into an error
     */
    private int findPair( int counter, int structure, RNA strand )
                          throws Exception {
        int pair = strand.GetPair( counter, structure );
        int error = strand.GetErrorCode();

        if( error != 0 ) {
            String message = strand.GetErrorMessage( error );
            RNAstructureInfoDialog.error( message );
            throw new Exception();
        }

        return pair;
    }

    /**
     * Write the structure to a helix file.
     */
    protected void writeFile() {
        try {
            FileSelector selector = new FileSelector( false, "save", type );
            String file = selector.showSelector();

            if( !file.equals( "" ) ) {
                int struct =
                    ((DrawStructure)((Imager)RNAstructure.getCurrentWindow())
                    .getSketcher()).getCurrentStructure();

                BufferedWriter writer =
                    new BufferedWriter( new FileWriter( file ) );
                RNA strand = RNAstructure.getCurrentStrand();

                int length = strand.GetSequenceLength();
                int counter = 1;

                while( counter <= length ) {
                    int pair = findPair( counter, struct, strand );

                    // Check for helices
                    if( pair > counter ) {
                        int count = 1;

                        while( findPair( counter + 1, struct, strand ) ==
                               findPair( counter, struct, strand ) - 1 ) {
                            counter++;
                            count++;
                        }

                        writer.write(
                            (counter - count + 1) + " " + pair + " " + count );
                        writer.newLine();
                        writer.flush();
                    }

                    counter++;
                }

                writer.close();
                RNAstructureInfoDialog.message( "Helix File Written." );
            }
        } catch( Exception e ) {
            RNAstructureInfoDialog.error( "Error Writing Helix File." );
        }
    }
}
