/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.Dimension;

import javax.swing.ButtonGroup;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

/**
 * A class that handles a dialog which allows the user to specify maximum
 * pairing distances between nucleotides
 */
public class MaxPairWindow extends DialogWindow {
    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -9141513326123205143L;

    /**
     * Constructor
     */
    public MaxPairWindow() {
        setTitle( "Maximum Distance For Pairs" );

        // Create radio button panel
        final ButtonGroup group = new ButtonGroup();
        JPanel panel1 = new JPanel();

        JRadioButton[] buttons = new JRadioButton[2];
        for( int i = 0; i < 2; i++ ) {
            JRadioButton button =
                new JRadioButton( ( i == 0 ) ? "Yes" : "No" );
            buttons[i] = button;
            group.add( button );
            panel1.add( button );
        }

        // Select the appropriate button in the panel
        int distance = RNAstructure.getCurrentStrand()
            .GetMaximumPairingDistance();
        int index = ( distance != -1 ) ? 0 : 1;
        group.setSelected( buttons[index].getModel(), true );

        panel1.setPreferredSize( new Dimension( 200, 75 ) );
        setBorder( panel1, 2, "Limit Distance Between Paired Bases:" );

        // Create input field
        final TextInputPanel panel2 = new TextInputPanel();
        String[] label = { "Maximum Distance" };
        panel2.create( 1, 2, 10, label );
        panel2.getFields()[0].setText( Integer.toString( distance ) );

        // Create pairing distance action
        ActionBase pairAction = new ChainedAction(
            new ForceDistance( panel2 ), new CloseAction() );

        // Create lower button panel
        ButtonPanel buttonPanel = new ButtonPanel( 2 );
        buttonPanel.setAction( pairAction, 0, "OK" );
        buttonPanel.setAction( new CloseAction(), 1, "Cancel" );
        buttonPanel.addButtons();

        // Add all componenents
        setFill( 1 );
        placeComponent( 0, 0, panel1 );
        placeComponent( 0, 1, panel2 );
        placeComponent( 0, 2, buttonPanel );
    }

    /**
     * An inner class which allows a MaxPairWindow to be constructed from the
     * more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * Create a MaxPairWindow.
         * @return  a new MaxPairWindow
         */
        public MaxPairWindow createWindow() {
            return new MaxPairWindow();
        }
    }
}