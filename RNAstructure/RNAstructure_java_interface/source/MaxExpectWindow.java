/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.BorderFactory;
import javax.swing.JTextField;

/**
 * A class that extends the basic FoldWindow, to generate structures based on a
 * maximum expected accuracy
 * @author Jessica Reuter
 */
public class MaxExpectWindow extends FoldWindow {
    /**
     * The TextInputPanel that holds the gamma value
     */
    private TextInputPanel gammaPanel;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -8160909075883113480L;

    /**
     * Constructor
     * @param acid  the type of nucleic acid this window is folding
     */
    public MaxExpectWindow( String acid ) {
        // Set window defaults
        super( acid + " Maximium Expected Accuracy" );
        boolean isRNA = acid.equals( "RNA" );

        // Create the main panel
        createMainPanel( 2, "Partition Function Save File", "CT File" );
        inputPanel.disableField( 1 );

	// Create the gamma text box
	gammaPanel = new TextInputPanel();
        gammaPanel.create( 1, 2, 10, "Gamma" );
        gammaPanel.setBorder( BorderFactory.createEmptyBorder( 0, 10, 0, 0 ) );
	gammaPanel.getFields()[0].setText( "1" );

        // Create the parameter panel
        createFoldingPanel();
	inputPanel2.getFields()[0].setText( "50" );
	inputPanel2.getFields()[1].setText( "1000" );
	inputPanel2.getFields()[2].setText( "5" );
	inputPanel2.getLabels()[0].setText( "Max % Score Difference" );

        // Create and set actions for buttons
        // First: Select the partition function save file and set its name.
        // Second: Create a strand from the partition function file.
        // Third: Set the name of the default output file.
        // Fourth: Select an output file other than default and set its name.
        // "First", "Second", "Third" execute in sequence on button 0.
        // "Fourth" executes alone on button 1.
        OpenSaveSetter first =
            OpenSaveSetter.open( 0, "Partition Function Save Files", "pfs" );
        StrandCreator second = StrandCreator.createPFSStrand( 0, isRNA );
        OutputFiller third = OutputFiller.fillSingle( 1, "ct" );
        OpenSaveSetter fourth = OpenSaveSetter.save( 1, "CT Files", "ct" );
        ChainedAction foldSequence =
            new ChainedAction( first, second, third );

        inputPanel.setAction( 0, foldSequence );
        inputPanel.setAction( 1, fourth );

        // Set components in place and designate the start button action
        setComponents( new PredictMaxAccuracy() );

	setFill( 1 );
	placeComponent( 0, 1, gammaPanel );
    }

    /**
     * Check if all required data has been placed in the window.
     * @return  true if all data is present, false if not
     */
    public boolean isReady() {
	return super.isReady() && gammaPanel.checkValues();
    }

    /**
     * Set the data for this window.
     * @throws Exception  if an error occurs in setting data
     */
    public void setData() throws Exception {
	removeData( 0 );
	addData( inputPanel.getFields()[1].getText() );

	JTextField[] values = inputPanel2.getFields();
	addData( Double.parseDouble( values[0].getText() ) );
	addData( Integer.parseInt( values[1].getText() ) );
	addData( Integer.parseInt( values[2].getText() ) );

	addData( Float.parseFloat( gammaPanel.getFields()[0].getText() ) );
    }

    /**
     * An inner class which allows a MaxExpectWindow to be constructed from the
     * more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * The type of nucleic acid to fold
         */
        private String acid;

        /**
         * Constructor
         * @param type  the type of nucleic acid to fold
         */
        public Factory( String type ) { acid = type; }

        /**
         * Create a MaxExpectWindow.
         * @return  a new MaxExpectWindow
         */
        public MaxExpectWindow createWindow() {
            return new MaxExpectWindow( acid );
        }
    }
}