/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.event.KeyEvent;

/**
 * A class that zooms by wrapping a ZoomHandler and linking it to a key event.
 * @author Jessica Reuter
 */
public class KeyZoomHandler extends KeyActionBase {
    /**
     * True if the control key is held down, false if not
     */
    private boolean isControl = false;

    /**
     * Activate the control key based on a key press.
     */
    public void pressKey() {
        if( event.isControlDown() ) { isControl = true; }
    }

    /**
     * Zoom based on a key release.
     */
    public void releaseKey() {
        int keyCode = event.getKeyCode();
        boolean enter =
            ( keyCode == KeyEvent.VK_RIGHT ) ||
            ( keyCode == KeyEvent.VK_LEFT );

        if( enter && isControl ) {
            SingleSketcher sketch =
                ((Imager)RNAstructure.getCurrentWindow()).getSketcher();
            double scale = sketch.getScale();

            if( keyCode == KeyEvent.VK_RIGHT ) { scale += 0.05; }
            else { scale -= 0.05; }

            if( scale < 0.05 ) { scale = 0.05; }
            else if( scale > 5 ) { scale = 5; }

            sketch.setScale( scale );
            new ZoomHandler( sketch ).act();
        } else {
            isControl = false;
            event.consume();
        }
    }
}
