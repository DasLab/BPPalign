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
 * A class that writes dot bracket structure files.
 * @author Jessica Reuter
 */
public class BracketStructureFileWriter extends StructureFileWriter {
    /**
     * Constructor
     */
    private BracketStructureFileWriter() {
	type = "Dot Bracket Files,bracket";
    }

    /**
     * Prepare to export the structure to a dot bracket file.
     * @return  a new BracketStructureFileWriter
     */
    public static BracketStructureFileWriter exportToDotBracketFile() {
        return new BracketStructureFileWriter();
    }

    /**
     * Write the structure to a dot bracket file.
     */
    protected void writeFile() {
        FileSelector selector = new FileSelector( false, "save", type );
        String file = selector.showSelector();

        if( !file.equals( "" ) ) {
            try {
                int code =
                    RNAstructure.getCurrentStrand().WriteDotBracket( file );
                boolean ok = !((Imager)RNAstructure.getCurrentWindow())
		    .getSketcher().getDataHolder().isErrorStatus( code );

                if( ok ) {
                    String message = "Dot Bracket File Written.";
                    RNAstructureInfoDialog.message( message );
                }
            } catch( Exception e ) {
                String error = "Error Writing Dot Bracket File.";
                RNAstructureInfoDialog.error( error );
            }
        }
    }
}
