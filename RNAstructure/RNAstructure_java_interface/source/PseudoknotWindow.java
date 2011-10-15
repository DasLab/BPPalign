/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JButton;
import javax.swing.JTextField;

/**
 * A class that opens a window allowing for pseudoknot breakage in an RNA
 * sequence prior to folding or other analysis.
 * @author Jessica Reuter
 */
public class PseudoknotWindow extends InputWindow {
    /**
     * The file name input panel
     */
    private TextInputPanel panel;

    /**
     * A serialized ID (required when extending GUI classes)
     */
    private static final long serialVersionUID = 7432098675411123751L;

    /**
     * Constructor
     */
    public PseudoknotWindow() {
        // Set titles
        setTitles( "Break Pseudoknots" );

        // Set the menu bar
        setSections( 1, 3, 4, 5, 7 );
        SmartMenuBar bar = RNAstructure.getBar();
        setMenus( bar.getMenuMaker().createTemperatureMenu() );
        bar.refresh( sections, menus );

        // Create the main panel and its actions
        panel = new TextInputPanel();
        panel.create( 2, 1, 30, "Input CT File", "Output CT File" );
        setBorder( panel, 1, null );

        // Create and set actions for buttons
        // First: Get the file to remove pseudoknots from and set its name
        // Second: Create a strand from the pseudoknotted file
        // Third: Fill in the default name for the output file
        // Fourth: Select an output file name other than the default and set it
        OpenSaveSetter first = OpenSaveSetter.open( 0, "CT Files", "ct" );
        StrandCreator second = StrandCreator.createCTStrand( 0, true );
        NoPseudoOutputFiller third = new NoPseudoOutputFiller();
        OpenSaveSetter fourth = OpenSaveSetter.save( 1, "CT Files", "ct" );

        ChainedAction pseudoSequence =
            new ChainedAction( first, second, third );

        panel.setAction( 0, pseudoSequence );
        panel.setAction( 1, fourth );

        // Add components in their proper places
        setGrid( 2, 1 );
        placeComponent( 0, 0, panel );

        setGrid( 1, 1 );
        setPad( 25, 25 );
        setInsets( 0, 0, 10, 0 );
        setAnchor( 2 );
        placeComponent( 1, 1, createStartButton( new RemovePseudoknots() ) );
    }

    /**
     * Check to make sure all required data has been entered into the window.
     * @return  true if all data is there, false if not
     */
    public boolean isReady() { return panel.checkValues(); }

    /**
     * Get all data required for breaking pseudoknots
     * @throws Exception  if setting of data runs into a problem
     */
    public void setData() throws Exception {
        JTextField[] fields = panel.getFields();

        removeData( 0 );
        addData( fields[1].getText() );
    }

    /**
     * An inner class which fills a specialized "no_pseudo" file name unique to
     * a PseudoknotWindow.
     * @author Jessica Reuter
     */
    class NoPseudoOutputFiller extends ActionBase {
        /**
         * Fill a file name into the window that doesn't hold pseudoknots.
         */
        public void act() {
            String original = panel.getFields()[0].getText();
            String newFile =
                original.replaceAll( "\\.ct", "\\_no\\_pseudo\\.ct" );
            panel.getFields()[1].setText( newFile );
        }
    }

    /**
     * An inner class which allows a PseudoknotWindow to be constructed from
     * the more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * Create a PseudoknotWindow.
         * @return  a new OligoScreenWindow
         */
        public PseudoknotWindow createWindow() {
            return new PseudoknotWindow();
        }
    }
}
