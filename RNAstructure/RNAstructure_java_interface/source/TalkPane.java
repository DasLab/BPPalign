/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.Color;
import java.util.LinkedList;
import java.util.List;

import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.Clip;
import javax.sound.sampled.DataLine;
import javax.swing.JTextPane;
import javax.swing.Timer;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.text.Document;
import javax.swing.text.StyledDocument;

/**
 * A text area with the ability to read its contents out loud. It can also
 * format its contents and control its acceptable input.
 * @author Jessica Reuter
 */
public class TalkPane extends JTextPane implements DocumentListener {
    /**
     * Whether the pane should speak when a character is entered into it
     */
    private boolean insert;

    /**
     * The position of the cursor in this pane.
     */
    private int position;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = 1156635818054833980L;

    /**
     * A timer that handles progressive reading of characters in the text pane
     */
    private Timer timer;

    /**
     * Constructor
     */
    public TalkPane() {
        // Set pane attributes
        setSelectionColor( Color.BLACK );
        setSelectedTextColor( Color.WHITE );
        position = 0;
        insert = false;

        // Create the timer and its action
        ActionBase timerAction = new ActionBase() {
            public void act() {
                try {
                    // If the caret position is still within the length of the
                    // sequence, highlight and read the next nucleotide
                    if( position != getText().length() ) {
                        getCaret().setSelectionVisible( false );
                        select( position, position + 1 );

                        Character base =
                            getSelectedText().toUpperCase().charAt( 0 );
                        if( base != ' ' ) {
                            getCaret().setSelectionVisible( true );
                            playClip( base );
                        }
                        position++;
                    }

                    // Otherwise if end has been reached, terminate the read
                    else {
                        Thread.sleep( 600 );
                        stopRead();
                        ((NewSequenceWindow)getTopLevelAncestor())
                            .switchTalkButton();
                    }
                } catch( Exception ex ) {
                    RNAstructureInfoDialog.error(
                        "Error during sequence playback." );
                    ((Timer)event.getSource()).stop();
                }
            }
        };
        timer = new Timer( 750, timerAction );

        // Add a key listener so only certain input types are accepted
        addKeyListener( new NucleotideSelector() );

        // Add a document listener for speaking out loud
        getDocument().addDocumentListener( this );
    }

    /**
     * Not used in a TalkPane
     */
    public void changedUpdate( DocumentEvent e ) {}

    /**
     * Format the sequence to make it easier to read.
     */
    public void format() {
        try {
            // Get unformatted text from this pane
            StringBuilder text = new StringBuilder(
                getText().replaceAll( "[\\s]+", "" ) );
            int length = text.length();

            // Create a list of lines 50 nucleotides long
            List<String> lines = new LinkedList<String>();
            for( int i = 0; i < length; i += 50 ) {
                boolean fullLine = length - i >= 50;
                int end = ( fullLine ) ? i + 50 : i;
                String line = text.substring( i, end );
                if( fullLine ) { lines.add( line ); }
                else { lines.add( text.substring( i ) ); }
            }

            // For each line, split into blocks of 5 nucleotides, and set each
            // in the text pane document
            setText( "" );

            int size = lines.size();
            for( int i = 0; i < size; i++ ) {
                StringBuilder builder = new StringBuilder( lines.get( i ) );
                for( int j = 5; j < builder.length(); j+=6 ) {
                    builder.insert( j, ' ' );
                }

                StyledDocument doc = getStyledDocument();
                String builderLine = builder.toString() + "\n";
                doc.insertString( doc.getLength(), builderLine, null );
            }

            // Reset the pane text as the formatted text, and set the caret at
            // the end of the text
            String fullText = getText().trim();
            setText( fullText );
            setCaretPosition( 0 );
        } catch( Exception e ) {
            RNAstructureInfoDialog.error( "Error formatting sequence." );
        }
    }

    /**
     * Insert a nucleotide into the documents, and say what it is. This only
     * happens if speaking while inserting is currently enabled.
     * @param e  the document event
     */
    public void insertUpdate( DocumentEvent e ) {
        if( insert ) {
            try {
                Document doc = e.getDocument();
                Character base = doc.getText( doc.getLength() - 1, 1 )
                    .toUpperCase().charAt( 0 ); 

                if( base != ' ' ) { playClip( base ); }
            } catch( Exception ex ) {
                RNAstructureInfoDialog.error(
                    "Error during sequence playback." );
            }
        }
    }

    /**
     * Get whether this pane is talking while typing.
     * @return  true if talking while typing, false if not
     */
    public boolean isTalkingWhileTyping() { return insert; }

    /**
     * Play a clip, depending on the nucleotide.
     * @param base  the nucleotide to play
     */
    public void playClip( char base ) {
        try {
            String file = Character.toString( base ).concat( ".wav" );

            AudioInputStream input =
                AudioSystem.getAudioInputStream(
                    TalkPane.class.getResourceAsStream( "/sounds/" + file ) );
            DataLine.Info info =
                new DataLine.Info( Clip.class, input.getFormat() );

            Clip clip = (Clip)AudioSystem.getLine( info );
            clip.open( input );
            clip.loop( 0 );
            Thread.sleep( 500 );
            clip.close();
        } catch( Exception e ) {
            RNAstructureInfoDialog.error( "Error during sequence playback." );
        }
    }

    /**
     * Not used in a TalkPane
     */
    public void removeUpdate( DocumentEvent e ) {}

    /**
     * Start reading the sequence programmatically
     */
    public void startRead() { timer.start(); }

    /**
     * Stop reading the sequence programmatically
     */
    public void stopRead() {
        position = 0;
        getCaret().setSelectionVisible( false );
        timer.stop();
    }

    /**
     * Set whether the pane speaks when a character is inserted. A call to this
     * method just switches "insert" from its current state to the opposite.
     */
    public void switchInsert() { insert = ( insert ) ? false : true; }
}