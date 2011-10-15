/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that holds all temporary dialogs which are not of the "option pane"
 * variety -- that is, all dialogs that are not covered by the class
 * "PopupDialog" and its wrapper "RNAstructureInfoDialog". See these two former
 * classes for details.
 * @author Jessica Reuter
 */
public class DialogWindow extends BaseWindow {
    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = 2576671861636179694L;

    /**
     * Constructor
     */
    protected DialogWindow() { super( RNAstructure.getCurrentWindow() ); }
}
