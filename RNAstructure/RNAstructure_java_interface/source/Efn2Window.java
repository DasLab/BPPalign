/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;

/**
 * A class respnsible for creating and displaying a window that enables a user
 * to input data for free energy calculations
 * @author Jessica Reuter
 */
public class Efn2Window extends InputWindow {
    /**
     * The check box which holds the option of writing a thermodynamic details
     * file if the user requests it
     */
    private JCheckBox box;

    /**
     * The main input panel
     */
    protected TextInputPanel inputPanel;

    /**
     * A serialized ID (required when extending GUI classes)
     */
    private static final long serialVersionUID = 8739104736254987302L;

    /**
     * Constructor
     * @param acid  the type of nucleic acid this window gives input for
     */
    public Efn2Window( String acid ) {
        // Set titles and create the action for the "START" button
        progress = false;
	setTitles( acid + " Efn2" );
        boolean isRNA = acid.equals( "RNA" );

        // Set the menu bar
        setSections( 1, 3, 4, 5, 7 );
        SmartMenuBar bar = RNAstructure.getBar();
        setMenus( bar.getMenuMaker().createTemperatureMenu() );
        bar.refresh( sections, menus );

        // Create the main input panel
        inputPanel = new TextInputPanel();
        inputPanel.create( 2, 1, 20, "CT File", "Output File" );
        setBorder( inputPanel, 1, null );
        inputPanel.disableField( 1 );

        // Create and set actions for buttons
        // First: Select a ct file for a strand, and set its name.
        // Second: Create a strand from the selected ct file
        // Third: Set the name of the strand's default output file.
        // Fourth: Select an output file other than default and set its name.
        // "First", "Second", and "Third" execute in sequence on button 0.
        // "Fourth" executes alone on button 1.
        OpenSaveSetter first = OpenSaveSetter.open( 0, "CT Files", "ct" );
        StrandCreator second = StrandCreator.createCTStrand( 0, isRNA );
        OutputFiller third = OutputFiller.fillSingle( 1, "out" );
        OpenSaveSetter fourth = OpenSaveSetter.save( 1, "OUT Files", "out" );
        ChainedAction efn2Sequence = new ChainedAction( first, second, third );

        inputPanel.setAction( 0, efn2Sequence );
        inputPanel.setAction( 1, fourth );

        // Create a check box for thermodynamic details
        box = new JCheckBox( "Write Thermodynamic Details File" );
        box.setBorder( BorderFactory.createEmptyBorder( 0, 0, 10, 0 ) );

        // Create the "START" button
        JButton startButton = createStartButton( new CalculateFreeEnergy() );

        // Add the components in their proper places
        setFill( 1 ); 
        setGrid( 2, 1 );
        placeComponent( 0, 0, inputPanel );

        setFill( 2 );
        placeComponent( 0, 1, box );

        setAnchor( 1 );
        setGrid( 1, 1 );
        setPad( 25, 25 );
        setInsets( 10, 0, 0, 10 );
        placeComponent( 2, 0, startButton );
    }

    /**
     * Check the input panel to see if all required files are present. This
     * method is a wrapper to call the TextInputPanel method checkValues; see
     * TextInputPanel.java for more information.
     * @return  true if all files are defined, false if not
     */
    public boolean isReady() { return inputPanel.checkValues(); }

    /**
     * Add the data into the data list for the strand this window is creating.
     */
    public void setData() {
        removeData( 0 );
        addData( inputPanel.getFields()[0].getText() );
        addData( inputPanel.getFields()[1].getText() );
        addData( box.isSelected() );
    }

    /**
     * An inner class which allows an Efn2Window to be constructed from the
     * more generic context of the RNAWindowFactory.
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
         * Create a Efn2Window.
         * @return  a new Efn2Window
         */
        public Efn2Window createWindow() { return new Efn2Window( acid ); }
    }
}
