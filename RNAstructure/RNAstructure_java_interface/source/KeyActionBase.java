/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;

/**
 * A base class for all RNAstructure GUI actions that originate from pressing
 * keys on the keyboard
 * @author Jessica Reuter
 */
public abstract class KeyActionBase extends KeyAdapter {
    /**
     * The key event this listener spawns
     */
    protected static KeyEvent event;

    /**
     * Do an action when a key is pressed.
     */
    public void keyPressed( KeyEvent e ) {
        event = e;
        pressKey();
    }

    /**
     * Do an action when a key is released
     */
    public void keyReleased( KeyEvent e ) {
        event = e;
        releaseKey();
    }

    @Override
    public void keyTyped( KeyEvent e ) {
        event = e;
        typeKey();
    }

    /**
     * Press a key and do something.
     * The base class version of this method doesn't do anything.
     */
    protected void pressKey() {}

    /**
     * Release a key and do something.
     * The base class version of this method doesn't do anything.
     */
    protected void releaseKey() {}

    /**
     * Type a key and do something.
     * The base class version of this method doesn't do anything.
     */
    protected void typeKey() {}
}