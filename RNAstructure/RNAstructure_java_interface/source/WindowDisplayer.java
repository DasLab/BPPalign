/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JDialog;

/**
 * A class which handles viewing of dialogs directly from a menu.
 * @author Jessica Reuter
 */
public class WindowDisplayer extends ActionBase {
    /**
     * The RNAWindowFactory this listener uses to create a window on the fly
     */
    private RNAWindowFactory factory;

    /**
     * Constructor
     * Allows for windows to be created directly from a menu
     * @param factory  the RNAWindowFactory this listener handles
     */
    public WindowDisplayer( RNAWindowFactory factory ) {
        this.factory = factory;
    }

    /**
     * View an object from a menu.
     */
    public void act() {
        RNAstructure.setLabel( "For help, press F1." );

        final Object window = factory.createWindow();
        if( window instanceof JDialog ) {
            JDialog dialog = (JDialog)window;
            dialog.pack();
            dialog.setLocationRelativeTo( RNAstructure.getFrame() );
            dialog.setVisible( true );
        }

	if( window instanceof Imager ) {
	    if( ((Imager)window).getSketcher() instanceof DrawStructure ) {
		DrawStructure sketcher =
		    (DrawStructure)(((Imager)window).getSketcher());
		double scale = sketcher.getScale();

		ZoomHandler zoomer = new ZoomHandler( sketcher );
		sketcher.setScale( scale - 1 );
		zoomer.act();
		sketcher.setScale( scale );
		zoomer.act();
	    }
	}
    }
}