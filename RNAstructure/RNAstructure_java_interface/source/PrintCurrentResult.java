/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that handles the act of printing a result from a ResultsWindow.
 * @author Jessica Reuter
 */
public class PrintCurrentResult extends ActionBase {
    /**
     * Print a result
     */
    protected void act() {
        ResultsWindow result = (ResultsWindow)RNAstructure.getCurrentWindow();
        result.getResultsPanel().printPanel();
    }
}
