/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that formats a nucleotide sequence in a TalkPane.
 * @author Jessica Reuter
 */
public class Formatter extends ActionBase {
    /**
     * Format a sequence
     */
    protected void act() {
        NewSequenceWindow window =
	    (NewSequenceWindow)RNAstructure.getCurrentWindow();
        TalkPane pane = (TalkPane)window.getArea( 3 );
        boolean speaking = pane.isTalkingWhileTyping();
        if( speaking ) { pane.switchInsert(); }
        pane.format();
        if( speaking ) { pane.switchInsert(); }
    }
}
