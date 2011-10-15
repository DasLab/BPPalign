/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.Color;
import java.awt.Container;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Window;

import javax.swing.BorderFactory;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.border.Border;

/**
 * A base dialog class which holds all common components of the RNAstructure
 * GUI dialogs/windows (except RNAstructureInfoDialog and PopupDialog)
 * @author Jessica Reuter
 */
public class BaseWindow extends JDialog {
    /**
     * The GridBagConstraints for this dialog
     */
    private GridBagConstraints constraints;

    /**
     * The pane of the dialog, set to GridBagLayout
     */
    private Container pane;

    /**
     * A serialized ID (required when extending GUI classes)
     */
    private static final long serialVersionUID = 8739104295948287302L;

    /**
     * Constructor
     * @param parent  the parent of this BaseWindow
     */
    protected BaseWindow( Window parent ) {
        super( parent );
        
        pane = getContentPane();
        pane.setLayout( new GridBagLayout() );
        constraints = new GridBagConstraints();

        setDefaultCloseOperation( JDialog.DISPOSE_ON_CLOSE );
        setResizable( false );
    }

    /**
     * Set a component in its place
     * @param x  the x coordinate of the grid
     * @param y  the y coordinate of the grid
     * @param component  the component to set in place
     */
    public void placeComponent( int x, int y, JComponent component ) {
        constraints.gridx = x;
        constraints.gridy = y;
        pane.add( component, constraints );
    }

    /**
     * Set the anchoring area for the GridBagConstraints
     * @param a  the anchor area
     */
    public void setAnchor( int a ) {
        if( a == 1 ) { constraints.anchor = GridBagConstraints.NORTH; }
        else { constraints.anchor = GridBagConstraints.CENTER; }
    }

    /**
     * Create a border for a component in the window
     * @param component  the component whose border to set
     * @param type  the type of border to set
     *              1. Complex TextInputPanel border (surrounds text fields)
     *              2. Simple Titled Border
     * @param title  the title of the border, if applicable
     */
    public void setBorder( JComponent component, int type, String title ) {
        Border border = ( type == 1 ) ?
            BorderFactory.createLineBorder( Color.LIGHT_GRAY, 1 ) :
            BorderFactory.createTitledBorder( title );
        int size = ( type == 1 ) ? 5 : 2;

        component.setBorder(
            BorderFactory.createCompoundBorder(
                BorderFactory.createEmptyBorder( 10, 10, 10, 10 ),
                BorderFactory.createCompoundBorder(
                    border, BorderFactory.createEmptyBorder(
                        size, size, size, size ) ) ) );
    }

    /**
     * Set the fill for the GridBagConstraints.
     * @param fill  the fill type
     */
    public void setFill( int fill ) {
        if( fill == 1 ) { constraints.fill = GridBagConstraints.HORIZONTAL; }
        else { constraints.fill = GridBagConstraints.CENTER; }
    }

    /**
     * Set the amount of space a component takes up in the grid.
     * @param width  the width of the component
     * @param height  the height of the component
     */
    public void setGrid( int width, int height ) {
        constraints.gridwidth = width;
        constraints.gridheight = height;
    }

    /**
     * Set the external padding of the constraints.
     * @param top  the top padding
     * @param left  the left padding
     * @param bottom  the bottom padding
     * @param right  the right padding
     */
    public void setInsets( int top, int left, int bottom, int right ) {
        constraints.insets = new Insets( top, left, bottom, right );
    }

    /**
     * Set padding of components.
     * @param xPad  the padding in the x direction
     * @param yPad  the padding in the y direction
     */
    public void setPad( int xPad, int yPad ) {
        constraints.ipadx = xPad;
        constraints.ipady = yPad;
    }

    /**
     * Set the title of both this window and the main window
     * @param title  the new title
     */
    public void setTitles( String title ) {
        setTitle( title );
        RNAstructure.getFrame().setTitle( "RNAstructure -- " + title );
    }

    /**
     * Set the relative weight of components (how much space they should fill
     * relative to their minimum size)
     * @param weightx  the weight in the x direction
     * @param weighty  the weight in the y direction
     */
    public void setWeights( double weightx, double weighty ) {
        constraints.weightx = weightx;
        constraints.weighty = weighty;
    }
}