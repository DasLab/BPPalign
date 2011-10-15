/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that opens an input window whose main source of nucleic acid data is
 * an Oligowalk_object.
 * @author Jessica Reuter
 */
public abstract class OligoWindow extends InputWindow {
    /**
     * The Oligowalk_object
     */
    protected Oligowalk_object oligoObject;

    /**
     * The file input panel
     */
    protected TextInputPanel panel;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -4459437214775988738L;

    /**
     * Get the file input panel.
     * @return  the file input panel
     */
    protected TextInputPanel getFileInputPanel() { return panel; }

    /**
     * Get the OligoWalk object.
     * @return  the OligoWalk object
     */
    protected Oligowalk_object getOligoObject() { return oligoObject; }

    /**
     * Set the OligoWalk object.
     * @param obj  the new OligoWalk object
     */
    protected void setOligoObject( Oligowalk_object obj ) {
        oligoObject = obj;
    }
}
