/*
 * (c) 2010  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JButton;
import javax.swing.JTextField;

/**
 * A class that opens a window allowing for identification of pseudoknots in
 * a nucleic acid sequence.
 * @author Jessica Reuter
 */
public class ProbKnotWindow extends InputWindow {
    /**
     * The file name input panel
     */
    private TextInputPanel panel;

    /**
     * The option values input panel
     */
    private TextInputPanel panel2;

    /**
     * A serialized ID (required when extending GUI classes)
     */
    private static final long serialVersionUID = 7481463058653521351L;

    /**
     * Constructor
     */
    public ProbKnotWindow() {
        // Set titles
        setTitles( "Identify Pseudoknots" );

        // Set the menu bar
        RNAstructure.getBar().refresh();

        // Create the main panel and its actions
        panel = new TextInputPanel();
        panel.create( 2, 1, 30, "Partition Function Save File", "CT File" );
        setBorder( panel, 1, null );

        // Create and set actions for buttons
        // First: Get the file to identify pseudoknots in and set its name
        // Second: Create a strand from the pseudoknotted file
	// Third: Fill the option fields with default values
        // Fourth: Fill in the default name for the output file
        // Fifth: Select an output file name other than the default and set it
        OpenSaveSetter first =
	    OpenSaveSetter.open( 0, "Partition Function Save Files", "pfs" );
        StrandCreator second = StrandCreator.createPFSStrand( 0, true );
	FoldingParamFiller third = FoldingParamFiller.setConstant( "1", "3" );
        OutputFiller fourth = OutputFiller.fillSingle( 1, "ct" );
        OpenSaveSetter fifth = OpenSaveSetter.save( 1, "CT Files", "ct" );

        ChainedAction pseudoSequence =
            new ChainedAction( first, second, third, fourth );

        panel.setAction( 0, pseudoSequence );
        panel.setAction( 1, fifth );

        // Create the options input panel.
        panel2 = new TextInputPanel();
        panel2.create( 2, 2, 15, "Iterations", "Minimum Helix Length" );
	panel2.getFields()[0].setText( "1" );
	panel2.getFields()[1].setText( "3" );
        setBorder( panel2, 1, null );

        // Add components in their proper places
        setGrid( 2, 1 );
        placeComponent( 0, 0, panel );

        setGrid( 1, 1 );
        placeComponent( 0, 1, panel2 );

        setPad( 25, 25 );
        setInsets( 0, 0, 10, 0 );
        setAnchor( 2 );
        placeComponent( 1, 1, createStartButton( new FindPseudoknots() ) );
    }

    /**
     * Get the options panel
     * @return  the options panel
     */
    public TextInputPanel getProbKnotOptionsPanel() { return panel2; }

    /**
     * Check to make sure all required data has been entered into the window.
     * @return  true if all data is there, false if not
     */
    public boolean isReady() {
        return panel.checkValues() && panel2.checkValues();
    }

    /**
     * Get all data required for breaking pseudoknots
     * @throws Exception  if setting of data runs into a problem
     */
    public void setData() throws Exception {
        JTextField[] fields = panel.getFields();
	JTextField[] fields2 = panel2.getFields();

        removeData( 0 );
        addData( fields[1].getText() );
	addData( Integer.parseInt( fields2[0].getText() ) );
	addData( Integer.parseInt( fields2[1].getText() ) );
    }

    /**
     * An inner class which allows a ProbKnotWindow to be constructed from the
     *  more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * Create a ProbKnotWindow.
         * @return  a new ProbKnotWindow
         */
        public ProbKnotWindow createWindow() {
            return new ProbKnotWindow();
        }
    }
}
