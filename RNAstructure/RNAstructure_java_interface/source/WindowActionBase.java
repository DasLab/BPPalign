/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;

/**
 * A base class for all RNAstructure GUI actions that originate from windows.
 * That is, those actions originating from windows themselves, not components
 * (ie. buttons) on them
 * @author Jessica Reuter
 */
public class WindowActionBase implements WindowListener {
    /**
     * The window event this listener spawns
     */
    protected static WindowEvent event;

    /**
     * Activate a window and do something.
     * The base class version of this method doesn't do anything.
     */
    public void activateWindow() {}

    /**
     * Close a window and do something afterward
     * The base class version of this method doesn't do anything.
     */
    public void closeWindowDoAfter() {}

    /**
     * Close a window and do something while the window is closing
     * The base class version of this method doesn't do anything.
     */
    public void closeWindowDoDuring() {}

    /**
     * Deactivate a window and do something.
     * The base class version of this method doesn't do anything.
     */
    public void deactivateWindow() {}

    /**
     * Deiconify a window and do something.
     * The base class version of this method doesn't do anything.
     */
    public void deiconifyWindow() {}

    /**
     * Iconify a window and do something.
     * The base class version of this method doesn't do anything.
     */
    public void iconifyWindow() {}

    /**
     * Open a window and do something.
     * The base class version of this method doesn't do anything.
     */
    public void openWindow() {}

    /**
     * Do an action when a window is activated
     * @param e  the window event
     */
    public void windowActivated( WindowEvent e ) {
        event = e;
        activateWindow();
    }

    /**
     * Do an action just after a window is closed
     * @param e  the window event
     */
    public void windowClosed( WindowEvent e ) {
        event = e;
        closeWindowDoAfter();
    }

    /**
     * Do an action while a window is closing.
     * @param e  the window event
     */
    public void windowClosing( WindowEvent e ) {
        event = e;
        closeWindowDoDuring();
    }

    /**
     * Do an action when a window is deactivated
     * @param e  the window event
     */
    public void windowDeactivated( WindowEvent e ) {
        event = e;
        deactivateWindow();
    }

    /**
     * Do an action when a window is deiconified.
     * @param e  the window event
     */
    public void windowDeiconified( WindowEvent e ) {
        event = e;
        deiconifyWindow();
    }

    /**
     * Do an action when a window is iconified.
     * @param e  the window event
     */
    public void windowIconified( WindowEvent e ) {
        event = e;
        iconifyWindow();
    }

    /**
     * Do an action when a window is opened.
     * @param e  the window event
     */
    public void windowOpened( WindowEvent e ) {
        event = e;
        openWindow();
    }
}