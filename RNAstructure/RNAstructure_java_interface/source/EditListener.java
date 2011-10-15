/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that monitors a text area on a NewSequenceWindow for changes.
 * @author Jessica Reuter
 */
public class EditListener extends KeyActionBase {
    /**
     * At a change in the text area, set the NewSequenceWindow's "edited"
     * property to true.
     */
    public void typeKey() {
        ((NewSequenceWindow)RNAstructure.getCurrentWindow()).setEdited( true );
    }
}
