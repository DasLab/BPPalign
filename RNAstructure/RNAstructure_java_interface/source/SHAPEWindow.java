/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that opens one of two SHAPE constraints windows, depending on the
 * type of constraint specified: hard constraints or pseudo-energy constraints
 * @author Jessica Reuter
 */
public class SHAPEWindow extends DialogWindow {
    /**
     * The button panel holding "OK" and "Cancel" buttons
     */
    protected ButtonPanel buttonPanel;

    /**
     * The main input panel
     */
    protected TextInputPanel panel;

    /**
     * The constraints input panel
     */
    protected TextInputPanel panel2;

    /**
     * A serialized ID (required when extending GUI classes)
     */
    private static final long serialVersionUID = 4920712953409128348L;

    /**
     * Constructor
     */
    protected SHAPEWindow() {;
        setTitle( "Read SHAPE Reactivity Data" );

        // Create the data file input panel
        String[] label = { "SHAPE Data File" };
        panel = new TextInputPanel();
        panel.create( 1, 1, 20, label );
        panel.setAction(
            0, OpenSaveSetter.open( 0, "SHAPE Data Files", "shape" ) );

        // Create the constraints input panel
        // The specific contents of this panel are handled in the subclasses.
        panel2 = new TextInputPanel();

        // Add the button panel
        buttonPanel = new ButtonPanel( 2 );
        buttonPanel.setAction( null, 0, "OK" );
        buttonPanel.setAction( new CloseAction(), 1, "Cancel" );
    }

    /**
     * Set the components of this SHAPEWindow in their places.
     */
    protected void setComponents() {
        buttonPanel.addButtons();

        // Add all components
	setFill( 2 );
        setPad( 0, 10 );
        placeComponent( 0, 0, panel );
        placeComponent( 0, 1, panel2 );
        placeComponent( 0, 2, buttonPanel );
    }
}