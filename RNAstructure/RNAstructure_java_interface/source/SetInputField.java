/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JTextField;

/**
 * A class that resets a text field to its empty state.
 * @author Jessica Reuter
 */
public class SetInputField extends ActionBase {
    /**
     * The text field to reset
     */
    private JTextField field;

    /**
     * The value to set in the field
     */
    private String text;

    /**
     * Reset constructor
     * @param field  the field to set
     */
    public SetInputField( JTextField field ) {
        this.field = field;
        text = "";
    }

    /**
     * Text Constructor
     * @param field  the field to set
     * @param text  the value to set
     */
    public SetInputField( JTextField field, String text ) {
        this.field = field;
        this.text = text;
    }

    /**
     * Set the field
     */
    protected void act() { field.setText( text ); }
}