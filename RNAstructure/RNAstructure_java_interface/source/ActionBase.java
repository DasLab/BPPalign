/* 
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * A base class for all actions in the RNAstructure GUI.
 * @author Jessica Reuter
 */
public abstract class ActionBase implements ActionListener {
    /**
     * The action event spawned for this action.
     */
    protected static ActionEvent event;

    /**
     * Do an action (in ActionBase subclasses)
     */
    protected abstract void act();

    /**
     * Do an action.
     * @param e  the action event
     */
    public void actionPerformed( ActionEvent e ) {
        event = e;
        act();
    }
}