/* 
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.CardLayout;

import javax.swing.JPanel;

/**
 * A class that switches between graph views in an Oligowalk results window.
 * @author Jessica Reuter
 */
public class OligoGraphSwitcher extends ActionBase {
    /**
     * The type of switcher
     */
    private int type;

    /**
     * Constructor
     * @param type  the type of switcher
     */
    private OligoGraphSwitcher( int type ) { this.type = type; }

    /**
     * Switch between graph views.
     */
    public void act() {
	OligoWalkResultsWindow window =
	    (OligoWalkResultsWindow)RNAstructure.getCurrentWindow();

	JPanel graphPanel = window.getGraphPanel();
	CardLayout layout = (CardLayout)graphPanel.getLayout();

	String layer =
	    ( type == 1 ) ? "6" :
	    ( type == 2 ) ? "2" :
	    ( type == 3 ) ? "3" :
	    ( type == 4 ) ? "1" :
	    ( type == 5 ) ? "and" :
	    ( type == 6 ) ? "4" :
	    null;
	
	layout.show( graphPanel, layer );

	OligoGraphPanel visible = window.getVisiblePanel();
	int current = visible.getCurrentlySelectedIndex();

	window.setLabels( current );
	visible.getBar( current ).select();
    }

    /**
     * Create an action that shows the bimolecular oligo graph.
     * @return  the bimolecular oligo action
     */
    public static OligoGraphSwitcher bimolecular() {
	return new OligoGraphSwitcher( 1 );
    }

    /**
     * Create an action that shows the broken target structures oligo graph.
     * @return  the broken target structures oligo action
     */
    public static OligoGraphSwitcher breakTargetStructure() {
	return new OligoGraphSwitcher( 2 );
    }

    /**
     * Create an action that shows the duplex oligo graph.
     * @return  the duplex oligo action
     */
    public static OligoGraphSwitcher duplex() {
	return new OligoGraphSwitcher( 3 );
    }

    /**
     * Create an action that shows the overall oligo graph.
     * @return  the overall oligo action
     */
    public static OligoGraphSwitcher overall() {
	return new OligoGraphSwitcher( 4 );
    }

    /**
     * Create an action that shows the overall/duplex oligo graph.
     * @return  the overall/duplex oligo action
     */
    public static OligoGraphSwitcher overallAndDuplex() {
	return new OligoGraphSwitcher( 5 );
    }

    /**
     * Create an action that shows the unimolecular oligo graph.
     * @return  the unimolecular oligo action
     */
    public static OligoGraphSwitcher unimolecular() {
	return new OligoGraphSwitcher( 6 );
    }
}
