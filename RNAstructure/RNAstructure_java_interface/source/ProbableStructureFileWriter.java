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
 * A class that writes probable pair structure files.
 * @author Jessica Reuter
 */
public class ProbableStructureFileWriter extends StructureFileWriter {
    /**
     * Constructor
     */
    private ProbableStructureFileWriter() {
        type = "CT Files,ct";
    }

    /**
     * Prepare to export the structure to a probable structures
     * file.
     * @return  a new ProbableStructureFileWriter
     */
    public static ProbableStructureFileWriter exportToProbableFile() {
        return new ProbableStructureFileWriter();
    }

    /**
     * Write structures to a probable pairs ct file.
     */
    protected void writeFile() {
	String errorString = "Error Writing Probable Structures File.";

        try {
            FileSelector selector = new FileSelector( false, "save", type );
            String file = selector.showSelector();

            if( !file.equals( "" ) ) {
                RNA strand = RNAstructure.getCurrentStrand();

                // Calculate structures
                int error = strand.PredictProbablePairs();
		if( error != 0 ) {
		    errorString = strand.GetErrorMessage( error );
		    throw new Exception();
		}

		// Write CT file
		error = strand.WriteCt( file );
		if( error != 0 ) {
		    errorString = strand.GetErrorMessage( error );
                    throw new Exception();
                }

                RNAstructureInfoDialog.message(
                    "Probable Structures File Written." );
            }
        } catch( Exception e ) { RNAstructureInfoDialog.error( errorString ); }
    }
}
