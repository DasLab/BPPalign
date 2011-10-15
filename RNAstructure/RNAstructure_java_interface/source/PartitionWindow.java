/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JButton;

/**
 * A class that creates an input window for executing the partition function on
 * one or multiple strands of nucleic acids.
 * @author Jessica Reuter
 */
public abstract class PartitionWindow extends InputWindow {
    /**
     * The input panel holding file names
     */
    protected TextInputPanel inputPanel;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -4293715023645273479L;

    /**
     * Constructor
     * @param title  the title of the window
     */
    protected PartitionWindow( String title ) { setTitles( title ); }

    /**
     * Create the main input panel.
     * @param amount  the number of buttons in the panel
     * @param names  the list of button labels
     */
    protected void createFilePanel( int amount, String... names ) {
        inputPanel = new TextInputPanel();
        setBorder( inputPanel, 1, null );
        inputPanel.create( amount, 1, 30, names );
    }

    /**
     * Check the input panel to see if all required files are present. This
     * method is a wrapper to call the TextInputPanel method checkValues; see
     * TextInputPanel.java for more information.
     * @return  true if all file names are defined, false if not
     */
    protected boolean isReady() { return inputPanel.checkValues(); }

    /**
     * Set components into their proper places
     * @param action  the action to place on the button
     */
    protected void setComponents( StrandStarter action ) {
        setGrid( 2, 1 );
        placeComponent( 0, 0, inputPanel );

        JButton startButton = createStartButton( action );

        setGrid( 1, 1 );
        setPad( 25, 25 );
        setInsets( 0, 0, 10, 0 );
        setAnchor( 2 );
        placeComponent( 1, 1, startButton );
    }

    /**
     * Get all data required for partitioning
     * @throws Exception  if setting of data runs into a problem
     */
    protected abstract void setData() throws Exception;
}
