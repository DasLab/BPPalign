/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

/**
 * A base class for all RNAstructure GUI actions that originate from the direct
 * action of a mouse.
 * @author Jessica Reuter
 */
public class MouseActionBase implements MouseListener {
    /**
     * The mouse event this listener spawns
     */
    protected static MouseEvent event;

    /**
     * Click the mouse and do something.
     * The base class version of this method doesn't do anything.
     */
    public void clickMouse() {}

    /**
     * Move the mouse into a region or component and do something.
     * The base class version of this method doesn't do anything.
     */
    public void enterMouse() {}

    /**
     * Move the mouse out of a region or component and do something.
     * The base class version of this method doesn't do anything.
     */
    public void exitMouse() {}

    /**
     * Do an action when the mouse is clicked.
     * @param e  the mouse event
     */
    public void mouseClicked( MouseEvent e ) {
        event = e;
        clickMouse();
    }

    /**
     * Do an action when the mouse enters a region or component.
     * @param e  the mouse event
     */
    public void mouseEntered( MouseEvent e ) {
        event = e;
        enterMouse();
    }

    /**
     * Do an action when the mouse exits a region or component.
     * @param e  the mouse event
     */
    public void mouseExited( MouseEvent e ) {
        event = e;
        exitMouse();
    }

    /**
     * Do an action when the mouse is pressed.
     * @param e  the mouse event
     */
    public void mousePressed( MouseEvent e ) {
        event = e;
        pressMouse();
    }

    /**
     * Do an action when the mouse is released.
     * @param e  the mouse event
     */
    public void mouseReleased( MouseEvent e ) {
        event = e;
        releaseMouse();
    }

    /**
     * Press the mouse and do something.
     * The base class version of this method doesn't do anything.
     */
    public void pressMouse() {}

    /**
     * Release the mouse and do something.
     * The base class version of this method doesn't do anything.
     */
    public void releaseMouse() {}
}