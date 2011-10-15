/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.Vector;

import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JToolBar;

/**
 * A class which holds the toolbar unique to the RNAstructure main window.
 * @author Jessica Reuter
 */
public class QuickAccessBar extends JToolBar {
    /**
     * A vector of buttons which holds the contents of this bar.
     */
    private Vector<JButton> buttons;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -8961355029870939220L;

    /**
     * Constructor
     * @throws IOException  if resource acquisition for the toolbar failed
     */
    public QuickAccessBar() throws IOException {
        // Initialize button vector
        buttons = new Vector<JButton>();

        // Create buttons and separators
        create( "NewSequence.gif",
                "New Sequence",
                new WindowDisplayer(
                    new NewSequenceWindow.Factory( "New Sequence" ) ) );

        create( "OpenSequence.gif",
                "Open Sequence",
                new SequenceOpener() );

        create( "Save.gif",
                "Save",
                new SequenceSaver( true ) );

        addSeparator();

        create( "Print.gif",
                "Print",
                new PrintCurrentResult() );

        addSeparator();

        create( "Draw.gif",
                "Draw",
                new WindowDisplayer( new Imager.StructureFactory() ) );

        create( "FoldRNASingle.gif",
                "Fold RNA Single Strand",
                new WindowDisplayer(
                    new FoldSingleWindow.Factory( "RNA" ) ) );

        create( "OligoWalk.gif",
                "RNA OligoWalk",
                new WindowDisplayer(
                    new OligoWalkWindow.Factory( "RNA" ) ) );

        create( "Dynalign.gif",
                "RNA Dynalign",
                new WindowDisplayer(
                    new DynalignFoldWindow.Factory( "RNA" ) ) );
    }

    /**
     * Create a button for the toolbar.
     * @param path  the path to the image on the button
     * @param text  the tooltip text of the button
     * @param act  the action listener added to the button
     * @throws IOException  if resource acquisition for the toolbar failed
     */
    private void create( String path, String text, ActionBase act )
                 throws IOException {
        // Create the image for the button
        BufferedImage image = ImageIO.read(
            QuickAccessBar.class.getResourceAsStream( "/images/" + path ) );

        // Create the button
        JButton button = new JButton();
        button.addActionListener( act );
        button.setName( text );
        button.setIcon( new ImageIcon( image ) );
        button.setToolTipText( text );

        // Add the button to the toolbar and to the vector
        add( button );
        buttons.add( button );
    }

    /**
     * Set the print and save buttons either enabled or disabled. These are the
     * only two buttons whose enabled statuses fluctuate
     */
    public void setButtonsEnabled( boolean one, boolean two ) {
        buttons.get( 2 ).setEnabled( one );
        buttons.get( 3 ).setEnabled( two );
    }
}