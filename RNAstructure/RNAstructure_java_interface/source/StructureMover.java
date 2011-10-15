/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.event.KeyEvent;

/**
 * A class that sets a particular structure and zooms to it.
 * @author Jessica Reuter
 */
public class StructureMover extends KeyActionBase {
    /**
     * True if the control key is pressed down, false otherwise
     */
    private boolean isControl = false;

    /**
     * Set the control key with a key press.
     */
    public void pressKey() {
        if( event.isControlDown() ) { isControl = true; }
    }

    /**
     * Set a structure and move to it on a key release.
     */
    public void releaseKey() {
        int keyCode = event.getKeyCode();
        boolean enter =
            ( keyCode == KeyEvent.VK_UP ) ||
            ( keyCode == KeyEvent.VK_DOWN );

        if( enter && isControl ) {
            DrawStructure sketch =
                (DrawStructure)(((Imager)RNAstructure.getCurrentWindow())
                    .getSketcher());
            int current = sketch.getCurrentStructure();
            int structures =
                RNAstructure.getCurrentStrand().GetStructureNumber();

            if( keyCode == KeyEvent.VK_DOWN ) { current--; }
            else { current++; }

            if( current < 1 ) { current = 1; }
            else if( current > structures ) { current = structures; }

            sketch.setStructureNumber( current );
            new ZoomHandler( sketch ).act();
        } else {
            isControl = false;
            event.consume();
        }
    }
}
