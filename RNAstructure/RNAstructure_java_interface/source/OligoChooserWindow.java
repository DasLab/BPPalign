/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.Dimension;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;

/**
 * A class that holds a dialog used to choose a particular oligo from the
 * results window after Oligowalk finishes.
 * @author Jessica Reuter
 */
public class OligoChooserWindow extends DialogWindow {
    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -8379039950797960535L;
    
    /**
     * Constructor
     * @param oligos  the number of oligos to choose from
     */
    public OligoChooserWindow( int oligos ) {
        // Create the spinner label and spinner
        JLabel label = new JLabel( "Base Number" );
        
	OligoWalkResultsWindow parent = (OligoWalkResultsWindow)getOwner();
        int index = parent.getVisiblePanel().getCurrentlySelectedIndex();
        SpinnerModel model = new SpinnerNumberModel( index, 1, oligos, 1 );
        JSpinner spinner = new JSpinner( model );

        // Create the top button
        JButton closeButton = new JButton( "Most Stable" );
        closeButton.addActionListener(
            GraphListener.createMaxOligoSpinLabeler( spinner ) );
        closeButton.setPreferredSize( new Dimension( 125, 30 ) );
        
        // Create the lower button panel
        ButtonPanel buttonPanel = new ButtonPanel( 2 );
        buttonPanel.setAction(
            new ChainedAction(
                GraphListener.createSelectedOligoChanger( spinner ),
                new CloseAction() ),
            0, "OK" );
        buttonPanel.setAction( new CloseAction(), 1, "Cancel" );
        buttonPanel.addButtons();

        JPanel closePanel = new JPanel();
        closePanel.add( closeButton );
        closePanel.setBorder( BorderFactory.createEmptyBorder( 5, 5, 5, 5 ) );
        
        // Place components
        placeComponent( 0, 0, label );
        placeComponent( 0, 1, spinner );
        placeComponent( 0, 2, buttonPanel );
        placeComponent( 0, 3, closePanel );
    }
}
