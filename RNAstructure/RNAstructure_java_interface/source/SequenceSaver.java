/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.io.BufferedWriter;
import java.io.FileWriter;

/**
 * A class that handles saving of a nucleotide sequence into a sequence file.
 * Regardless of the original format of the sequence, the output of this class
 * is always a sequence file.
 * @author Jessica Reuter
 */
public class SequenceSaver extends ActionBase {
    /**
     * Whether an error occurred during sequence saving
     */
    private boolean error;

    /**
     * The name of the file being saved
     */
    private String saveFile;

    /**
     * Whether the sequence should be saved with its existing name or a new one
     */
    private boolean saveWithNew;

    /**
     * Constructor
     * @param newName  true if the sequence should be saved with a new name,
     *                 false if not
     */
    public SequenceSaver( boolean newName ) { saveWithNew = newName; }

    /**
     * Save a sequence.
     */
    public void act() {
        try {
            RNAstructure.setLabel( "For help, press F1." );

            error = false;
            saveFile = getCurrentFileName();

            // If the file is to be saved under a new name, select a new name.
            // Otherwise, get the current file name to save under.
            if( saveWithNew || saveFile.equals( "New Sequence" ) ) {
                FileSelector selector =
                    new FileSelector( true, "Save", "Sequence Files,seq" );
                saveFile = selector.showSelector();
            }

            // Write the sequence file
            if( !saveFile.equals( "" ) ) {
                writeSequenceFile( saveFile );
		NewSequenceWindow window =
		    (NewSequenceWindow)RNAstructure.getCurrentWindow();

                window.setTitles( saveFile );
                window.setEdited( false );
            } else { error = true; }
        } catch( Exception ex ) {
            RNAstructureInfoDialog.error( "Error saving sequence." );
            error = true;
        }
    }

    /**
     * Get the current file name that must be saved. Since this action is only
     * used when a NewSequenceWindow is the most recent window, we can simply
     * get the title off that window.
     * @return  the current file name
     */
    private String getCurrentFileName() {
        return ((NewSequenceWindow)RNAstructure.getCurrentWindow()).getTitle();
    }

    /**
     * Get the file that this saver is saving.
     * @return  the file name
     */
    public String getSaveFile() { return saveFile; }

    /**
     * Check if an error occurred during sequence saving.
     * @return  true if an error occurred, false if not
     */
    public boolean isError() { return error; }

    /**
     * Write a sequence file.
     * @param file  the file name to write to
     * @throws Exception  all exceptions are handled by actionPerformed
     */
    private void writeSequenceFile( String file ) throws Exception {
        // Get variables
        BufferedWriter writer = new BufferedWriter( new FileWriter( file ) );
        NewSequenceWindow window =
	    (NewSequenceWindow)RNAstructure.getCurrentWindow();

        // Write comments
        String[] comments = window.getArea( 2 ).getText().split( "\n" );
        for( String element: comments ) {
            writer.write( "; " + element );
            writer.newLine();
            writer.flush();
        }

        // If there are no comments, write an empty comment to make sure
        // the file format is still maintained.
        if( comments.length == 0 ) {
            writer.write( "; " );
            writer.newLine();
            writer.flush();
        }

        // Write identifier
        String identifier = window.getArea( 1 ).getText();
        if( identifier.length() == 0 ) { identifier = " "; }
        writer.write( identifier );
        writer.newLine();
        writer.flush();

        // Write sequence
        String[] sequence = window.getArea( 3 ).getText().split( "\n" );
        for( String element: sequence ) {
            writer.write( element );
            writer.newLine();
            writer.flush();
        }

        // Write termination mark and close the writer
        writer.write( "1" );
        writer.close();
    }
}
