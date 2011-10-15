/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.GridLayout;

import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

/**
 * A class that opens the OligoScreen input window.
 * @author Jessica Reuter
 */
public class OligoScreenWindow extends OligoWindow {
    /**
     * The button group showing oligomer chemistry selection
     */
    private ButtonGroup group;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -4459432886935988738L;

    /**
     * Constructor
     */
    public OligoScreenWindow() {
        // Set titles and defaults
        progress = false;
        setTitles( "OligoScreen" );

        // Set the menu bar
        setSections( 1, 3, 4, 5, 7 );
        SmartMenuBar bar = RNAstructure.getBar();
        setMenus( bar.getMenuMaker().createTemperatureMenu() );
        bar.refresh( sections, menus );

        // Create the main input panel
        String[] labels = { "Oligomer List", "Output Filename" };
        panel = new TextInputPanel();
        panel.create( 2, 1, 30, labels );
        setBorder( panel, 1, null );

        // Create and set actions on panel
        OpenSaveSetter first =
            OpenSaveSetter.open( 0, "List Files", "lis" );
        OligoWalkObjectCreator second =
            new OligoWalkObjectCreator( true );
        OutputFiller third = OutputFiller.fillSingle( 1, "rep" );
        OpenSaveSetter fourth =
            OpenSaveSetter.save( 1, "Report Files", "rep" );

        ChainedAction screenAction =
            new ChainedAction( first, second, third );

        panel.setAction( 0, screenAction );
        panel.setAction( 1, fourth );

        // Create the oligomer chemistry panel
        JPanel panel2 = new JPanel( new GridLayout( 1, 2 ) );
        group = new ButtonGroup();
        String[] chemistry = { "DNA", "RNA" };
        for( int i = 0; i < 2; i++ ) {
            JRadioButton button = new JRadioButton( chemistry[i] );
            group.add( button );
            if( i == 1 ) { button.setSelected( true ); }

            button.addActionListener( new ActionBase() {
                public void act() {
                    String acid = event.getActionCommand();
                    boolean rna = acid.equals( "RNA" );
                    oligoObject.setIsrna( rna );
                }
            });

            panel2.add( button );
        }
        setBorder( panel2, 2, "Oligomer Chemistry" );

        // Create the start button
        JButton startButton = createStartButton( new RunOligoScreen() );

        // Set all components in their proper places
        setGrid( 2, 1 );
        setFill( 1 );
        placeComponent( 0, 0, panel );

        setGrid( 1, 1 );
        setFill( 2 );
        setPad( 40, 0 );
        setInsets( 0, 70, 10, 10 );
        placeComponent( 0, 1, panel2 );

        setPad( 25, 25 );
        setAnchor( 3 );
        setInsets( 0, 0, 10, 70 );
        placeComponent( 1, 1, startButton );
    }

    /**
     * Check to make sure all required file names and data have been entered in
     * this window.
     * @return  true if all data is there, false if not
     */
    public boolean isReady() {
        if( !panel.checkValues() ) { return false; }
        if( group.getSelection() == null ) { return false; }
        return true;
    }

    /**
     * Get the data from the OligoScreen window.
     */
    public void setData() {
        removeData( 0 );
        addData( panel.getFields()[0].getText() );
        addData( panel.getFields()[1].getText() );
    }

    /**
     * An inner class which allows an OligoScreenWindow to be constructed from
     * the more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * Create an OligoScreenWindow.
         * @return  a new OligoScreenWindow
         */
        public OligoScreenWindow createWindow() {
            return new OligoScreenWindow();
        }
    }
}
