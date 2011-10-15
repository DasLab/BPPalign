/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.io.BufferedReader;
import java.io.FileReader;

import javax.swing.text.StyledDocument;

/**
 * A class that handles opening of a sequence in a NewSequenceWindow.
 * @author Jessica Reuter
 */
public class SequenceOpener extends ActionBase {
    /**
     * The NewSequenceWindow a file is being opened into
     */
    private NewSequenceWindow window;

    /**
     * Open a sequence for analysis.
     * @param e  the action event
     */
    public void act() {
        RNAstructure.setLabel( "For help, press F1." );

        // Open a file chooser to select the file the user wants
        FileSelector selector = new FileSelector( false, "Open",
            "Genbank Files,gen", "Plain Text Files,txt",
            "Sequence Files,seq" );
        String file = selector.showSelector();

        // If the user selects a file (rather than cancels), open it
        if( !file.equals( "" ) ) {
            try {
                // Initialize the window and the file reader
                window = new NewSequenceWindow( file );
                BufferedReader reader =
                    new BufferedReader( new FileReader( file ) );

                // Read the file into the window
                if( file.endsWith( "seq" ) ) { readSequenceFile( reader ); }
                else if( file.endsWith( "gen" ) ) {
                    readGenbankFile( reader );
                } else { readTextFile( reader ); }
                reader.close();

                // Show the window
                window.pack();
                window.setLocationRelativeTo( RNAstructure.getFrame() );
                window.setVisible( true );
            } catch( Exception ex ) {
                RNAstructureInfoDialog.error( "Error opening sequence." );
            }
        }
    }

    /**
     * Read a Genbank file into a NewSequenceWindow
     * @param reader  the BufferedReader that reads this file
     * @throws Exception  all exceptions are handled by actionPerformed
     */
    private void readGenbankFile( BufferedReader reader ) throws Exception {
        String line = null;

        // Skip the beginning of the file, then read definition as identifier
        while( ( !( line = reader.readLine() ).startsWith( "DEFINITION" ) ) ) {
            // As long as definition line hasn't been reached, discard input
        }

        String name = line.substring( line.indexOf( 'N' ) + 1 ).trim();
        while( !( line = reader.readLine() ).startsWith( "ACCESSION" ) ) {
            name = name.concat( line.trim() );
        }
        window.getArea( 1 ).getStyledDocument().insertString( 0, name, null );

        // Skip through the file until version reached, read version as comment
        while( ( !( line = reader.readLine() ).startsWith( "VERSION" ) ) ) {}
        String detail = line.substring( line.indexOf( ' ' ) ).trim();
        window.getArea( 2 ).getStyledDocument()
            .insertString( 0, detail, null );

        // Skip over rest of the file until origin reached, then read sequence
        while( !reader.readLine().trim().equals( "ORIGIN" ) ) {}

        StyledDocument sequenceDoc = window.getArea( 3 ).getStyledDocument();
        while( ( line = reader.readLine().trim() ) != "//" ) {
            int length = sequenceDoc.getLength();
            String nextLine = line.substring( 1 ) + "\n";
            sequenceDoc.insertString( length, nextLine, null );
        }
    }

    /**
     * Read a sequence file into a NewSequenceWindow.
     * @param reader  the BufferedReader that reads this file
     * @throws Exception  all exceptions are handled by actionPerformed
     */
    private void readSequenceFile( BufferedReader reader ) throws Exception {
        String line = null;

        // Read comments
        StyledDocument commentDoc = window.getArea( 2 ).getStyledDocument();
        while( ( line = reader.readLine().trim() ).startsWith( ";" ) ) {
            int length = commentDoc.getLength();
            String nextLine = line.substring( 1 ).trim() + "\n";
            commentDoc.insertString( length, nextLine, null );
        }

        // When broken out of the comment loop, we now have the name
        window.getArea( 1 ).getStyledDocument().insertString( 0, line, null );

        // Read the actual sequence
        StyledDocument sequenceDoc = window.getArea( 3 ).getStyledDocument();
        while( ( line = reader.readLine() ) != null ) {
            int length = sequenceDoc.getLength();
            String nextLine = line.replace( '1', ' ' ).trim() + "\n";
            sequenceDoc.insertString( length, nextLine, null );
        }
    }

    /**
     * Read a plain text file into a NewSequenceWindow.
     * @param reader  the BufferedReader that reads this file
     * @throws Exception  all exceptions are handled by actionPerformed
     */
    private void readTextFile( BufferedReader reader ) throws Exception {
        String line = null;

        StyledDocument sequenceDoc = window.getArea( 3 ).getStyledDocument();
        while( ( line = reader.readLine() ) != null ) {
            int length = sequenceDoc.getLength();
            String nextLine = line.replace( '1', ' ' ).trim() + "\n";
            sequenceDoc.insertString( length, nextLine, null );
        }
    }
}