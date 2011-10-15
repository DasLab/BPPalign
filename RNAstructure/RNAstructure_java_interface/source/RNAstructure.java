/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Font;
import java.awt.SplashScreen;
import java.awt.image.BufferedImage;
import java.awt.Window;

import java.util.LinkedList;

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.UIDefaults;
import javax.swing.UIManager;
import javax.swing.plaf.metal.MetalLookAndFeel;

/**
 * A class that starts the RNAstructure interface, creates the initial frame
 * @author Jessica Reuter
 */
public class RNAstructure {
    /**
     * The menu bar at the top of the main frame
     */
    private static SmartMenuBar bar;

    /**
     * The main application frame
     */
    private static JFrame frame;

    /**
     * The window that is currently focused.
     * This window may or may not be the most recent one.
     */
    private static Window focusedWindow;

    /**
     * The information label at the bottom of the frame
     */
    private static JLabel label;

    /**
     * The toolbar at the top of the main frame
     */
    private static QuickAccessBar toolBar;

    /**
     * The list of visible windows
     */
    private static LinkedList<MenuChangingWindow> windowList;

    /**
     * Constructor
     * Creates the first frame
     */
    public RNAstructure() {
        try {
            // Create frame
            frame = new JFrame( "RNAstructure" );
            frame.setLayout( new BorderLayout() );

            // Create window list
            windowList = new LinkedList<MenuChangingWindow>();

            // Set icon image
            BufferedImage image = ImageIO.read(
                RNAstructure.class
		.getResourceAsStream( "/images/Icon.gif" ) );
            frame.setIconImage( image );

            // Create back panel
            JPanel panel = new JPanel();
            panel.setBackground( Color.GRAY );
            frame.add( panel, BorderLayout.CENTER );

            // Add label that shows infotips
            label = new JLabel( " For help, press F1." );
            label.setFont( new Font( label.getFont().getFontName(), 0, 16 ) );
            frame.add( label, BorderLayout.SOUTH );

            // Set the menu bar and the toolbar
            bar = new SmartMenuBar();
            bar.refresh();
            frame.setJMenuBar( bar );

            toolBar = new QuickAccessBar();
            toolBar.setButtonsEnabled( false, false );
            frame.add( toolBar, BorderLayout.PAGE_START );

            // Set general frame properties and visiblity
            frame.setSize( 1024, 768 );
            frame.setDefaultCloseOperation( JFrame.EXIT_ON_CLOSE );
            frame.setLocationRelativeTo( null );
            frame.setVisible( true );
        } catch( Exception e ) {
            System.err.println( "Error creating RNAstructure GUI." );
            e.printStackTrace();
            System.exit( -1 );
        }
    }

    /**
     * Get the main RNAstructure main menu bar.
     * @return  the main menu bar
     */
    public static SmartMenuBar getBar() { return bar; }

    /**
     * Get the strand of nucleic acids the top input window is working with.
     * In this case, the top window is the one that's currently focused, not
     * the most recent one.
     * @return  the strand of nucleic acids linked to the top window.
     */
    public static RNA getCurrentStrand() {
        if( focusedWindow instanceof DynalignWindow ) {
            return ((DynalignWindow)focusedWindow).getWorkingStrand();
        } else if( focusedWindow instanceof OligoWindow ) {
            return ((OligoWindow)focusedWindow).getOligoObject();
        } else if( focusedWindow instanceof InputWindow ) {
	    return ((InputWindow)focusedWindow).getDataHolder().getStrand();
        } else if( focusedWindow instanceof Imager ) {
            return ((Imager)focusedWindow).getSketcher().getDataHolder()
		.getStrand();
        } else return null;
    }

    /**
     * Get the currently focused window
     * @return  the currently focused window
     */
    public static Window getCurrentWindow() { return focusedWindow; }

    /**
     * Get the main RNAstructure frame.
     * @return  the frame
     */
    public static JFrame getFrame() { return frame; }

    /**
     * Get the main RNAstructure toolbar.
     * @return  the main toolbar
     */
    public static QuickAccessBar getToolBar() { return toolBar; }

    /**
     * Get the window list.
     * @return  the window list
     */
    public static LinkedList<MenuChangingWindow> getWindowList() {
        return windowList;
    }

    /**
     * Set the currently focused window.
     * @param window  the new currently focused window
     */
    public static void setCurrentWindow( Window window ) {
	focusedWindow = window;
    }

    /**
     * Set the text of the main label to a particular value.
     * @param text  the new text
     */
    public static void setLabel( String text ) {
	label.setText( " " + text );
    }

    /**
     * The main method.
     * Set the look and feel according to the current system's native look and
     * feel, and read necessary data for GUI back end connection and display.
     * @param args  the command line arguments (ignored)
     */
    public static void main( String[] args ) {
        try {
            // Set the system look and feel to be as close to the current
            // computer as possible.

            // Since all icons aren't defined on all systems, loop through the
            // required icons for this program, and if they aren't defined set
            // them to the Java default.
            UIManager.setLookAndFeel( 
                UIManager.getSystemLookAndFeelClassName() );

            UIDefaults def = new MetalLookAndFeel().getDefaults();
            String[] iconList = {
              "OptionPane.errorIcon",
              "OptionPane.informationIcon",
              "OptionPane.warningIcon",
              "MenuItem.checkIcon"
            };

            for( String next: iconList ) {
                if( UIManager.getIcon( next ) == null ) {
                    UIManager.put( next, def.getIcon( next ) );
                }
            }

            // Call up the splash screen while the shared library loads
            SplashScreen splash = SplashScreen.getSplashScreen();
            System.loadLibrary( "RNAstructure_GUI" );
            Thread.sleep( 3000 );
            splash.close();
        } catch( Exception e ) {
            System.err.println( "Error loading RNAstructure." );
            e.printStackTrace();
            System.exit( -1 );
        }

        // Start the GUI
        new RNAstructure();
    }
}
