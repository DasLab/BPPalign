/* 
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that modifies an existing strand, but doesn't perform any complex
 * calculations or visualization.
 * @author Jessica Reuter
 */
public abstract class Modifier extends ActionBase {
    /**
     * Wrapper for modify method
     */
    public void act() { modify(); }

    /**
     * Modify the current strand.
     */
    protected abstract void modify();
}