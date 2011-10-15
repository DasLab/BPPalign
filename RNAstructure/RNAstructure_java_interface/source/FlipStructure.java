/* 
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.util.LinkedList;

/**
 * A class that flips a structure so it's drawn clockwise or counterclockwise.
 * @author Jessica Reuter
 */
public class FlipStructure extends ActionBase {
    /**
     * Constructor
     */
    private FlipStructure() {}

    /**
     * Flip a structure.
     */
    public void act() {
	RNAstructure.setLabel( "For help, press F1." );

	LinkedList<MenuChangingWindow> list = RNAstructure.getWindowList();
        for( MenuChangingWindow window: list ) {
            boolean target =
		window instanceof Imager &&
                ((Imager)window).getSketcher() instanceof DrawStructure;

            if( target ) {
		DrawStructure sketch =
		    (DrawStructure)((Imager)window).getSketcher();
                sketch.renderClockwise( !sketch.isClockwise() );
                break;
            }
        }
    }

    /**
     * Create an action that flips a structure.
     * @return  the flip action
     */
    public static FlipStructure flip() { return new FlipStructure(); }
}
