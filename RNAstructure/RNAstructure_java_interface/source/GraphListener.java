/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
*/

package source;

import java.awt.Component;
import java.awt.Container;
import java.awt.Rectangle;

import javax.swing.JComponent;
import javax.swing.JSpinner;

/**
 * A class which handles a variety of listeners for the results window shown
 * after Oligowalk completes.
 * @author Jessica Reuter
 */
public class GraphListener extends ActionBase {
    /**
     * The type of listener
     */
    private int type;

    /**
     * A secondary object affected by the listener or gives information to the
     * listener (may or may not be used)
     */
    private Object object;
    
    /**
     * Interact with the graph.
     */
    public void act() {
	// Obtain common objects
        Object source = event.getSource();
        Container root = ((JComponent)source).getTopLevelAncestor();

        OligoWalkResultsWindow window =
	    ( root instanceof OligoWalkResultsWindow ) ?
	        (OligoWalkResultsWindow)root :
	    ( root instanceof OligoChooserWindow ) ?
	        (OligoWalkResultsWindow)((OligoChooserWindow)root)
	        .getParent() :
	    null;

        // If this action is being used to interpret/change the graph directly
        // from the results window
        if( root instanceof OligoWalkResultsWindow ) {
	    OligoGraphPanel graph = window.getVisiblePanel();

            if( type == 0 ) {
		int selected = (Integer)object;
                int value = graph.getCurrentlySelectedIndex() + selected;
                
                if( selected != 0 ) { setAllIndices( window, value ); }
                else {
                    OligoChooserWindow dialog = new OligoChooserWindow(
			graph.getNumberBars() );
                    dialog.pack();
                    dialog.setVisible( true );
                }
            } else if( type == 1 ) {
		setAllIndices( window, ((GraphBar)source).getIndex() );
            }
        }

        // If this action is being used to interpret/change the graph from the
        // oligo chooser dialog
        else if( root instanceof OligoChooserWindow ) {
            JSpinner spin = (JSpinner)object;
            
            if( type == 2 ) { spin.setValue( window.getMostStableOligo() ); }
            else if( type == 3 ) {
		String value = spin.getValue().toString();
		setAllIndices( window, Integer.parseInt( value ) );
            }
        }
    }

    /**
     * Create a listener that highlights a particular bar.
     * @return  the bar highlighting listener
     */
    public static GraphListener createBarHighlighter() {
        GraphListener listener = new GraphListener();
        listener.setType( 1 );
        return listener;
    }

    /**
     * Create a listener that moves a particular increment (number of bars) on
     * the graph.
     * @param increment  the number of bars to move
     * @return  the increment movement listener
     */
    public static GraphListener createIncrementMover( int increment ) {
        GraphListener listener = new GraphListener();
        listener.setObject( increment );
        listener.setType( 0 );
        return listener;
    }

    /**
     * Create a listener that sets the current oligo to the most stable one.
     * @param spin  the spinner that selects an oligo
     * @return  the max oligo setting listener
     */
    public static GraphListener createMaxOligoSpinLabeler( JSpinner spin ) {
        GraphListener listener = new GraphListener();
        listener.setType( 2 );
        listener.setObject( spin );
        return listener;
    }

    /**
     * Create a listener that selects a particular oligo.
     * @param spin  the spinner that selects an oligo
     * @return  the oligo setting listener
     */
    public static GraphListener createSelectedOligoChanger( JSpinner spin ) {
        GraphListener listener = new GraphListener();
        listener.setType( 3 );
        listener.setObject( spin );
        return listener;
    }
    
    /**
     * Set a selected oligo in all relevant places, then scroll the results
     * window to that location
     * @param window  the results window
     * @param index  the index of the oligo
     */
    private void setAllIndices( OligoWalkResultsWindow window, int index ) {
	Component[] components = window.getGraphPanel().getComponents();
	OligoGraphPanel graph = window.getVisiblePanel();
	int length = graph.getNumberBars();
	
	if( index < 1 ) { index = 1; }
        else if( index > length ) { index = length; }
	
	for( Component element: components ) {
	    OligoGraphPanel panel = (OligoGraphPanel)element;
	    panel.getBar( panel.getCurrentlySelectedIndex() ).deselect();
            panel.setCurrentValue( index );
            panel.getBar( index ).select();
	}
	
	window.setOligo( index );
	window.setLabels( index );

        Rectangle visibleRect = window.getVisiblePanel().getVisibleRect();
        visibleRect.x = ( ( index * 10 ) - 5 ) - visibleRect.width / 2;
        window.getVisiblePanel().scrollRectToVisible( visibleRect );
    }
    
    /**
     * Set an object in the listener.
     * @param object  the object to set
     */
    private void setObject( Object object ) { this.object = object; }
    
    /**
     * Set the type of listener.
     * @param type  the type of listener.
     */
    private void setType( int type ) { this.type = type; }
}
