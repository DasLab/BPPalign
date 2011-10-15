/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Desktop;
import java.awt.Dimension;
import java.awt.image.BufferedImage;

import java.net.URI;

import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;

/**
 * A class that opens information windows for the RNAstructure interface.
 * @author Jessica Reuter
 */
public class HelpWindow extends DialogWindow {
    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -7257135450358148960L;

    /**
     * Constructor
     * @param help  true if opening a help search dialog,
     *              false if opening a window about the application
     */
    public HelpWindow( boolean help ) {
        setLayout( new BorderLayout() );

        if( help ) { help(); }
        else { about(); }

        pack();
        setDefaultCloseOperation( JDialog.DISPOSE_ON_CLOSE );
        setLocationRelativeTo( RNAstructure.getFrame() );
        setVisible( true );
    }

    /**
     * Open a window describing the application.
     */
    private void about() {
        try {
            setTitle( "About RNAstructure" );

            // Add the about image, which is the same as the splash screen
            BufferedImage image = ImageIO.read(
                HelpWindow.class.getResourceAsStream( "/images/Splash.gif" ) );
            add( new JLabel( new ImageIcon( image ) ) );

            // Create the bottom button panel
            ButtonPanel panel = new ButtonPanel( 1 );
            panel.setAction( new CloseAction(), 0, "OK" );
            panel.addButtons();
            add( panel, BorderLayout.SOUTH );
        } catch( Exception e ) {
            RNAstructureInfoDialog.error( "Error creating about screen." );
        }
    }

    /**
     * Open a help search dialog.
     */
    private void help() {
	try {
	    String page =
		"http://rna.urmc.rochester.edu/GUI/html/Contents.html";
	    Desktop.getDesktop().browse( new URI( page ) );
	} catch( Exception e ) {
	    RNAstructureInfoDialog.error( "Error initializing help." );
	}
    }

    /**
     * An inner class which allows a HelpWindow to be constructed from the more
     * generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * Boolean telling whether to open a help window or not
         */
        private boolean help;

        /**
         * Constructor
         * @param help  true if opening a help window, false if not
         */
        public Factory( boolean help ) { this.help = help; }

        /**
         * Create a HelpWindow.
         * @return  a new HelpWindow
         */
        public HelpWindow createWindow() { return new HelpWindow( help ); }
    }
}