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
import javax.swing.JPanel;

/**
 * A class that holds a bar graph of oligo values for the results window which
 * opens after Oligowalk completes.
 * @author Jessica Reuter
 */
class OligoGraphPanel extends JPanel {
    /**
     * The current index selected on the graph.
     */
    protected int current;

    /**
     * The graph on the panel.
     */
    protected Object[] graph;

    /**
     * Array of minimum and maximum bar lengths on this graph
     */
    protected double[] minAndMax;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -8143419599993674743L;

    /**
     * The Oligowalk_object that holds this graph's data
     */
    protected Oligowalk_object strand;

    /**
     * Constructor
     * @param strand  the Oligowalk_object whose information is shown
     * @param oligoLength   the oligo length
     */
    public OligoGraphPanel( Oligowalk_object strand, int oligoLength ) {
        this.strand = strand;

        setBackground( Color.WHITE );
        current = 1;

        graph = new GraphBar[strand.GetSequenceLength() - oligoLength + 1];
    }

    /**
     * Create a particular bar graph.
     * @param type  the type of bar graph being made
     * @param color  the color of the bars in the graph
     */
    public void createGraph( int type, Color color ) {
        removeAll();
        int span = 500;
        int graphSize = getNumberBars();

        graph = new GraphBar[graphSize];
        setLayout( new BoxLayout( this, BoxLayout.X_AXIS ) );

        minAndMax = findMinAndMax( type );
        double minValue = minAndMax[0];
        double maxValue = minAndMax[1];
        
        for( int i = 1; i <= graphSize; i++ ) {
             double value =
		 OligoWalkResultsWindow.getOligoValue( strand, type, i );
	    
	     double totalSpan = maxValue - minValue;
	     int minSpan = (int)(span * ( Math.abs( minValue ) / totalSpan ));
	     int maxSpan = (int)(span * ( maxValue / totalSpan ));
	     
	     int length = 0;
	     if( value < 0 ) {
		 length = -(int)(( value / minValue ) * minSpan);
	     } else { length = (int)(( value / maxValue ) * maxSpan); }
	    
	     GraphBar bar = new GraphBar( i, color, length );
	     graph[i-1] = bar;
	     add( bar );
        }
        repaint();
        
        Container top = ((GraphBar)graph[0]).getTopLevelAncestor();
        if( top != null ) {
            OligoWalkResultsWindow window = (OligoWalkResultsWindow)top;
            window.setLabels( current );
        }
    }

    /**
     * Find the min and max values for a particular graph.
     * @param type  the type of graph being made
     */
    public double[] findMinAndMax( int type ) {
        double min = Double.MAX_VALUE;
	double max = Double.MAX_VALUE * -1;

	int size = getNumberBars();
	for( int i = 1; i <= size; i++ ) {
	    double value =
		OligoWalkResultsWindow.getOligoValue( strand, type, i );

	    if( value > max ) { max = value; }
	    else if( value < min ) { min = value; }
	}
	if( max < 0 ) { max = 0; }

	double[] bounds = { min, max };
	return bounds;
    }

    /**
     * Get a bar from the graph.
     * @param index  the index of the bar.
     * @return  the bar
     */
    public GraphBar getBar( int index ) {
	return (GraphBar)graph[index-1];
    }
    
    /**
     * Get the index that's currently selected.
     * @return  the selected index
     */
    public int getCurrentlySelectedIndex() { return current; }

    /**
     * Get the bounds of the panel (min and max()
     * @return  an array of bounds
     */
    public double[] getGraphBounds() { return minAndMax; }
    
    /**
     * Get the number of bars in the graph.
     * @return  the number of bars
     */
    public int getNumberBars() { return graph.length; }
    
    /**
     * Set the current selected value.
     * @param  newValue  the new current value
     */
    public void setCurrentValue( int newValue ) { current = newValue; }
}
