/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.Color;
import java.awt.Dimension;

import javax.swing.BorderFactory;
import javax.swing.border.Border;
import javax.swing.border.LineBorder;
import javax.swing.JButton;

/**
 * A class that holds a bar button, used in Oligowalk bar graphs.
 * @author Jessica Reuter
 */
class GraphBar extends JButton {
    /**
     * The color of this button
     */
    private final Color color;

    /**
     * The index of this button in relation to the total number of oligos
     * (one-indexed)
     */
    private final int index;

    /**
     * The length of the bar in kcal/mol
     */
    private final int length;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -1180987435097504823L;
    
    /**
     * Constructor
     * @param index  the index of this bar
     * @param color  the color of this bar
     * @param length  the length of this bar
     */
    public GraphBar( int index, Color color, int length ) {
        this.color = color;
        this.index = index;
        this.length = length;

        setBar( color );
        addActionListener( GraphListener.createBarHighlighter() );
    }
    
    /**
     * Deselect this bar.
     */
    public void deselect() { setBar( color ); }

    /**
     * Get the index of this bar.
     * @return  the index of this bar
     */
    public int getIndex() { return index; }

    /**
     * Select this bar.
     */
    public void select() { setBar( Color.red ); }
    
    /**
     * Set the bar's size, color, and alignment.
     * @param color  the color of the bar.
     */
    public void setBar( Color color ) {
        setPreferredSize( new Dimension( 10, Math.abs( length ) ) );

	Border border = BorderFactory.createCompoundBorder(
                            new LineBorder( Color.black, 1 ),
			    new LineBorder( color, 4 ) );
	setBorder( border );

        setMaximumSize( getPreferredSize() );
        setMinimumSize( getPreferredSize() );
        
        if( length <= 0 ) { setAlignmentY( JButton.TOP_ALIGNMENT ); }
        else { setAlignmentY( JButton.BOTTOM_ALIGNMENT ); }
    }
}
