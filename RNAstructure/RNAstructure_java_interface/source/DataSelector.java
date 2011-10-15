/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that handles actions arising from buttons on a TextInputPanel
 * @author Jessica Reuter
 */
public abstract class DataSelector extends ActionBase {
    /**
     * Do an action (wrapper for specific version in this class and subclasses)
     */
    protected void act() { selectInput(); }

    /**
     * Select a type of input.
     */
    protected abstract void selectInput();
}