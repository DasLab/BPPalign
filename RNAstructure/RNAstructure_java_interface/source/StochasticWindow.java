/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JButton;
import javax.swing.JTextField;

/**
 * A class that opens a window allowing for stochastic probability analysis.
 * @author Jessica Reuter
 */
public class StochasticWindow extends InputWindow {
    /**
     * The main input panel
     */
    private TextInputPanel panel;

    /**
     * The secondary input panel for parameters
     */
    private TextInputPanel panel2;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = 478604612860973084L;

    /**
     * Constructor
     * @param acid  the type of nucleic acid being analyzed
     */
    public StochasticWindow( String acid ) {
        // Set titles and defaults
        setTitles( acid + " Stochastic Sampling" );

        // Set the menu bar
        RNAstructure.getBar().refresh();

        // Create the main input panel
        panel = new TextInputPanel();
        panel.create( 2, 1, 30, "Partition Function Save File", "CT File" );
        setBorder( panel, 1, null );
        panel.disableField( 1 );

        // Create and set actions for buttons
        // First: Select a partition function save file and set its name.
        // Second: Create a strand from the partition function save file. RNA
        //         can be hardcoded here because the save file handles the
        //         final nucleic acid type.
        // Third: Set the name of the default output file.
        // Fourth: Select an output file other than default, and set its name.
        // "First", "Second", and "Third" execute in sequence on button 0.
        // "Fourth" executes alone on button 1.
        OpenSaveSetter first =
            OpenSaveSetter.open( 0, "Partition Function Save Files", "pfs" );
        StrandCreator second = StrandCreator.createPFSStrand( 0, true );
        OutputFiller third = OutputFiller.fillSingle( 1, "ct" );
        OpenSaveSetter fourth = OpenSaveSetter.save( 1, "CT Files", "ct" );
        ChainedAction stochasticSequence =
            new ChainedAction( first, second, third );

        panel.setAction( 0, stochasticSequence );
        panel.setAction( 1, fourth );

        // Create the secondary input panel for parameters
        panel2 = new TextInputPanel();
        panel2.create( 2, 2, 15, "Ensemble Size", "Random Seed" );
        panel2.getFields()[0].setText( "1000" );
        panel2.getFields()[1].setText( "1234" );
        setBorder( panel2, 1, null );

        // Create the "START" button
        JButton startButton = createStartButton( new SampleStochastic() );

        // Add components in their proper places
        setGrid( 2, 1 );
        setFill( 1 );
        placeComponent( 0, 0, panel );

        setGrid( 1, 1 );
        setFill( 2 );
        setWeights( 1, 0 );
        placeComponent( 0, 1, panel2 );

        setPad( 25, 25 );
        setWeights( 0, 0 );
        setInsets( 0, 0, 0, 50 );
        placeComponent( 1, 1, startButton );
    }

    /**
     * Check the input panels to see if all required files and parameter values
     * are present. This method is a wrapper to call the TextInputPanel method
     * checkValues on both input panels; see TextInputPanel.java for more
     * information.
     * @return  true if all values are defined, false if not
     */
    public boolean isReady() {
        return ( panel.checkValues() && panel2.checkValues() );
    }

    /**
     * Set the data from this stochastic sampling input window
     * @return  an array of data
     */
    public void setData() {
        removeData( 0 );

        JTextField[] fields1 = panel.getFields();
        addData( fields1[1].getText() );

        JTextField[] fields2 = panel2.getFields();
        addData( Integer.parseInt( fields2[0].getText() ) );
        addData( Integer.parseInt( fields2[1].getText() ) );
    }

    /**
     * An inner class which allows a StochasticWindow to be constructed from
     * the more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * The nucleic acid the created window handles
         */
        private String acid;

        /**
         * Constructor
         * @param acid  the nucleic acid the window handles
         */
        public Factory( String acid ) { this.acid = acid; }

        /**
         * Create a StochasticWindow.
         * @return  a new StochasticWindow
         */
        public StochasticWindow createWindow() {
            return new StochasticWindow( acid );
        }
    }
}
