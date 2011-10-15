/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that reads a nucleic acid sequence out loud.
 * @author Jessica Reuter
 */
public class ReadAloud extends ActionBase {
    /**
     * A boolean, true if reading is automatic (programmatic) or manual (read
     * while typing)
     */
    private boolean automatic;

    /**
     * Constructor
     * @param type  the type of reading being done
     *              true = automatic
     *              false = manual
     */
    private ReadAloud( boolean automatic ) { this.automatic = automatic; }

    /**
     * Activate or deactivate reading of a sequence out loud.
     */
    protected void act() {
        NewSequenceWindow window =
	    (NewSequenceWindow)RNAstructure.getCurrentWindow();

        if( automatic ) {
            window.switchTalkButton();
            if( event.getActionCommand().equals( "Read Sequence" ) ) {
                ((TalkPane)window.getArea( 3 )).startRead();
            } else { ((TalkPane)window.getArea( 3 )).stopRead(); }
        } else { ((TalkPane)window.getArea( 3 )).switchInsert(); }
    }

    /**
     * Read aloud programmatically (automatically)
     * @return  the automatic ReadAloud
     */
    public static ReadAloud readAutomatic() { return new ReadAloud( true ); }

    /**
     * Read aloud manually (while typing in data)
     * @return  the manual ReadAloud
     */
    public static ReadAloud readWhileTyping() {
        return new ReadAloud( false );
    }
}
