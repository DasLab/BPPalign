/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

/**
 * A class that handles panels with text input fields on them, preceded either
 * by buttons or labels which describe the contents of the text field
 * @author Jessica Reuter
 */
public class TextInputPanel extends JPanel {
    /**
     * An array of buttons that this text panel contains
     * (if applicable)
     */
    private JButton[] buttons;

    /**
     * An array of text fields that this text panel contains
     */
    private JTextField[] fields;

    /**
     * An array of labels that this text panel contains
     * (if applicable)
     */
    private JLabel[] labels;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -7816818805397111948L;

    /**
     * Check to make sure all required values are defined. All values are
     * considered defined if all fields in this panel have a string inside them
     * whose length is greater than 0.
     * Fields which contain only whitespace are considered to have length 0.
     * @return  true if all field values are defined, false if not
     */
    public boolean checkValues() {
        int size = fields.length;
        for( int i = 0; i < size; i++ ) {
            if( fields[i].getText().trim().length() == 0 ) {
                RNAstructureInfoDialog.warning( "<html>One or more input " + 
                    "values are missing.<br/>Please specify all required " +
                    "values to continue." );
                return false;
            }
        }

        return true;
    }

    /**
     * Create the text input panel.
     * This method handles both major aspects of panel creation -- component
     * layout/look and feel and the addition of actions to buttons, where
     * necessary for a particular panel.
     * @param amount  the number of text fields in the panel
     * @param type  type of objects next to text fields (1=buttons, 2=labels)
     * @param size  the size (width) of the text fields
     * @param list  the data array which holds object names
     */
    public void create( int amount, int type, int size, String... list ) {
        // Set data arrays
        JComponent[] components = new JComponent[amount];
        fields = new JTextField[amount];
        if( type == 1 ) { buttons = new JButton[amount]; }
	else { labels = new JLabel[amount]; }

        // Create objects
        for( int i = 0; i < amount; i++ ) {
            JComponent component;
            if( type == 1 ) {
                component = new JButton( list[i] );
                buttons[i] = (JButton)component;
            } else {
		component = new JLabel( list[i] );
		labels[i] = (JLabel)component;
	    }
            components[i] = component;

            JTextField field = new JTextField();
            field.setColumns( size );
            fields[i] = field;
        }

        // Set layout and layout manager
        setLayout( new GridBagLayout() );
        GridBagConstraints constraints = new GridBagConstraints();
        constraints.fill = GridBagConstraints.HORIZONTAL;
        constraints.ipadx = 5;

        // Add components to the text panel
        for( int i = 0; i < amount; i++ ) {
            constraints.gridy = i;

            constraints.gridx = 0;
            add( components[i], constraints );

            constraints.gridx = 1;
            add( fields[i], constraints );
        }
    }

    /**
     * Disable a particular text field in the panel.
     * @param index  the index of the text field (one-indexed)
     */
    public void disableField( int index ) {
        fields[index-1].setEditable( false );
        fields[index-1].setBackground( Color.WHITE );
    }

    /**
     * Get the array of buttons for this panel, if applicable.
     * @return the button array
     */
    public JButton[] getButtons() { return buttons; }

    /**
     * Get the array of text fields for this panel
     * @return the text field array
     */
    public JTextField[] getFields() { return fields; }

    /**
     * Get the array of labels for this panel, if applicable.
     * @return the text label array
     */
    public JLabel[] getLabels() { return labels; }

    /**
     * Set the action on a particular button of this panel.
     * @param index  the index of the button in the panel
     * @param listener  the action listener for this button
     */
    public void setAction( int index, ActionBase listener ) {
        buttons[index].addActionListener( listener );
    }
}