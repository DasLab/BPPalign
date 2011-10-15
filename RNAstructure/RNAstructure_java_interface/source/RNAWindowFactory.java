/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * An interface holding a single factory method to show dialogs and assorted
 * windows in the RNAstructure GUI. Since some of these windows are implemented
 * as wrappers around very simple dialogs (because some preparation is needed),
 * and the method is used to call constructors, the method returns an object,
 * not a dialog. The name "RNAWindowFactory" refers to the intended use of the
 * objects returned, not the raw type of the objects themselves.
 * @author Jessica Reuter
 */
public interface RNAWindowFactory {
    /**
     * Create the window of choice.
     * @return  the window (as an object)
     */
    public Object createWindow();
}