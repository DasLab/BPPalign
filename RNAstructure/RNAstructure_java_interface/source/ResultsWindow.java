/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.BorderLayout;
import java.awt.Color;

/**
 * A class that holds a window which displays results of a calculation. These
 * results normally have menus associated with them that allow for exploration
 * of them.
 * @author Jessica Reuter
 */
public class ResultsWindow extends MenuChangingWindow {
    /**
     * The panel which holds the results
     */
    protected PrintablePanel resultsPanel;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = 3148741838989255765L;

    /**
     * Constructor
     */
    protected ResultsWindow() {
	getContentPane().setBackground( Color.white );
        setLayout( new BorderLayout() );
        setResizable( true );

        resultsPanel = new PrintablePanel();
        resultsPanel.setLayout( new BorderLayout() );
	resultsPanel.setBackground( Color.white );
        add( resultsPanel );
    }

    /**
     * Get the results panel from this window.
     * @return  the results panel
     */
    public PrintablePanel getResultsPanel() { return resultsPanel; }
}