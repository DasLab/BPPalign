/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.Dimension;
import java.awt.GridLayout;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JPanel;

/**
 * A class holding a gridded panel which holds one or more choice buttons
 * @author Jessica Reuter
 */
public class ButtonPanel extends JPanel {
    /**
     * The buttons in the panel
     */
    private JButton[] buttons;

    /**
     * A serialized ID (required when extending GUI classes)
     */
    private static final long serialVersionUID = 4920712953409128348L;

    /**
     * Constructor
     * @param number  the amount of buttons on this panel
     */
    public ButtonPanel( int number ) {
        setLayout( new GridLayout( 1, 0 ) );
        buttons = new JButton[number];
    }

    /**
     * Add the buttons to the panel with their proper sizing and gaps
     */
    public void addButtons() {
        for( int i = 0; i < buttons.length; i++ ) {
            JButton button = buttons[i];
            button.setPreferredSize( new Dimension( 125, 30 ) );

            JPanel buttonPanel2 = new JPanel();
            buttonPanel2.add( button );
            buttonPanel2.setBorder( 
                BorderFactory.createEmptyBorder( 10, 20, 10, 20 ) );
            add( buttonPanel2 );
        }
    }

    /**
     * Set an action on a button
     * @param action  the proper action
     * @param index  the index of the button whose action to set
     * @param text  the text of the button to set
     */
    public void setAction( ActionBase action, int index, String text ) {
        JButton button = new JButton( text );
        button.addActionListener( action );
        buttons[index] = button;
    }

    /**
     * Set actions on multiple buttons at once. Note that the i-th button in
     * the panel is assigned the i-th action and i-th text label.
     * @param acts  the array of action listeners
     * @param texts  the array of button labels
     */
    public void setActions( ActionBase[] acts, String[] texts ) {
        int number = acts.length;
        for( int i = 0; i < number; i++ ) {
            setAction( acts[i], i, texts[i] );
        }
    }

    /**
     * Set a particular button in the button array as a unit.
     * @param i  the index of the button to set
     * @param button  the button to set
     */
    public void setButton( int i, JButton button ) { buttons[i] = button; }
}