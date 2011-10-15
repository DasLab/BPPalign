/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.Color;
import java.awt.Container;

import javax.swing.Box;
import javax.swing.BoxLayout;

/**
 * A class which holds a stacked OligoGraphPanel for simultaneous viewing of
 * overall and duplex values.
 * @author Jessica Reuter
 */
public class StackedOligoGraphPanel extends OligoGraphPanel {
    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -8470026029704075420L;

    /**
     * Constructor
     * @param strand  the Oligowalk_object the panel is displaying
     * @param length  the oligo length
     */
    public StackedOligoGraphPanel( Oligowalk_object strand, int length ) {
	super( strand, length );
	setLayout( new BoxLayout( this, BoxLayout.Y_AXIS ) );
    }
    
    /**
     * Create the stacked graph.
     */
    public void createStackedGraph() {
	removeAll();
        int span = 500;
        
        int bars = getNumberBars();
        graph = new Box[bars][];
	for( int i = 0; i < bars; i++ ) {
	    Box[] column =
		{ Box.createVerticalBox(), Box.createVerticalBox() };

	    graph[i] = column;
	}
        
        Color blue = Color.blue;
        Color green = new Color( 17, 167, 77 );
        
        Box topBox = Box.createHorizontalBox();
	Box bottomBox = Box.createHorizontalBox();
	
	topBox.setAlignmentX( Box.LEFT_ALIGNMENT );
	bottomBox.setAlignmentX( Box.LEFT_ALIGNMENT );

	for( int i = 0; i < graph.length; i++ ) {
	    (((Box[])graph[i])[0]).setAlignmentY( Box.BOTTOM_ALIGNMENT );
	    topBox.add( ((Box[])graph[i])[0] );

	    (((Box[])graph[i])[1]).setAlignmentY( Box.TOP_ALIGNMENT );
	    bottomBox.add( ((Box[])graph[i])[1] );
	}
        
        double minValue = Double.MAX_VALUE, maxValue = Double.MAX_VALUE * -1;
	for( int i = 1; i < graph.length; i++ ) {
            double next =
                OligoWalkResultsWindow.getOligoValue( strand, 1, i );
	    double next2 =
		OligoWalkResultsWindow.getOligoValue( strand, 3, i );

	    double min = Math.min( next, next2 );
	    double max = Math.max( next, next2 );

            if( min < minValue ) { minValue = min; }
	    else if( max > maxValue ) { maxValue = max; }
        }

	minAndMax = new double[2];
	if( maxValue < 0 ) { maxValue = 0; }
	minAndMax[0] = minValue;
	minAndMax[1] = maxValue;
        
        for( int i = 0; i < graph.length; i++ ) {
	    int current = i + 1;
	    double value =
		OligoWalkResultsWindow.getOligoValue( strand, 1, current );
	    double value2 =
		OligoWalkResultsWindow.getOligoValue( strand, 3, current );
	    
	    double totalSpan = maxValue - minValue;
	    int minSpan = (int)(span * ( Math.abs( minValue ) / totalSpan ));
	    int maxSpan = (int)(span * ( maxValue / totalSpan ));
	    
	    int length = ( value < 0 ) ?
		-(int)(( value / minValue ) * minSpan) :
		(int)(( value / maxValue ) * maxSpan);
	    
	    int length2 = ( value2 < 0 ) ?
		-(int)(( value2 / minValue ) * minSpan) :
		(int)(( value2 / maxValue ) * maxSpan);

	    if( (length > 0 && length2 > 0) || (length < 0 && length2 < 0) ) {
		boolean bigLen = ( value > 0 ) ?
		    ( length > length2 ) : ( length < length2 );
		
		int base = ( bigLen ) ? length2 : length;
		int clip = ( bigLen ) ? length - length2 : length2 - length;
		
		if( value > 0 ) {
		    ((Box[][])graph)[i][0]
			.add( new GraphBar( i + 1, green, clip ) );
		    ((Box[][])graph)[i][0]
			.add( new GraphBar( i + 1, blue, base ) );
		    ((Box[][])graph)[i][1]
			.add( new GraphBar( i + 1, green, 0 ) );
		} else {
		    ((Box[][])graph)[i][0]
			.add( new GraphBar( i + 1, green, 0 ) );
		    ((Box[][])graph)[i][1]
			.add( new GraphBar( i + 1, blue, base ) );
		    ((Box[][])graph)[i][1]
			.add( new GraphBar( i + 1, green, clip ) );
		}
	    } else {
		((Box[][])graph)[i][0]
		    .add( new GraphBar( i + 1, blue, length ) );
		((Box[][])graph)[i][1]
		    .add( new GraphBar( i + 1, green, length2 ) );
	    }
        }
        
        add( topBox );
	add( bottomBox );
        repaint();
        
        Container top = ((Box[])graph[0])[0].getTopLevelAncestor();
        getBar( getCurrentlySelectedIndex() ).select();
        if( top != null ) {
	    OligoWalkResultsWindow window = (OligoWalkResultsWindow)top;
	    window.setLabels( current );
        }
    }
    
    /**
     * Get a bar from the graph.
     * @param index  the index of the bar.
     */
    public GraphBar getBar( int index ) {
	Box[] column = ((Box[][])graph)[index-1];
	
	if( column[0].getComponentCount() > 1 ) {
	    return (GraphBar)column[0].getComponent( 1 );
	} else { return (GraphBar)column[1].getComponent( 0 ); }
    }
}
