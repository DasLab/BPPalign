/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JButton;

/**
 * A class that handles refolding of previous Dynalign calculations from
 * Dynalign save files.
 * @author Jessica Reuter
 */
public class DynalignRefoldWindow extends DynalignWindow {
    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = 8308768309194385130L;

    /**
     * Constructor
     */
    public DynalignRefoldWindow() {
        // Set titles
        super( "Refold From Dynalign Save File" );
	progress = true;
	progressIndeterminate = true;

        // Set the menu bar
        RNAstructure.getBar().refresh();

        // Create the file input panel
        createFilePanel( 4, "Save File", "CT File 1", "CT File 2",
            "Alignment File" );

        // Create the structure parameters panel
        createStructurePanel();

        // Create and set actions on buttons
        // First: Select the Dynalign save file name and set it
	// Second: Set the default folding parameters
        // Third: Select the first CT file and set its name
        // Fourth: Select the second CT file and set its name
        // Fifth: Select the alignment file and set its name.
        // "First" and "Second" execute in sequence on button 0.
        // "Third" executes alone on button 1.
        // "Fourth" executes alone on button 2.
        // "Fifth" executes alone on button 3.
        OpenSaveSetter first =
            OpenSaveSetter.open( 0, "Dynalign Save Files", "dsv" );
	FoldingParamFiller second =
	    FoldingParamFiller.setConstant( "20", "750", "0", "0" ); 
        OpenSaveSetter third =
            OpenSaveSetter.save( 1, "CT Files", "ct" );
        OpenSaveSetter fourth =
            OpenSaveSetter.save( 2, "CT Files", "ct" );
        OpenSaveSetter fifth =
            OpenSaveSetter.save( 3, "Alignment Files", "ali" );

	ChainedAction refolder = new ChainedAction( first, second );

        inputPanel.setAction( 0, refolder );
        inputPanel.setAction( 1, third );
        inputPanel.setAction( 2, fourth );
        inputPanel.setAction( 3, fifth );

        // Set components in their proper places
        setFill( 1 );
        setGrid( 2, 1 );
        placeComponent( 0, 0, inputPanel );

        setGrid( 1, 1 );
        placeComponent( 0, 1, inputPanel2 );
        setGrid( 1, 1 );

        setFill( 2 );
        setPad( 25, 25 );
        setInsets( 0, 0, 10, 0 );
        placeComponent( 1, 1,
            createStartButton( new RefoldDynalignRunner() ) );
    }

    /**
     * An inner class which allows a DynalignRefoldWindow to be constructed
     * from the more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * Create a DynalignRefoldWindow.
         * @return  a new DynalignRefoldWindow
         */
        public DynalignRefoldWindow createWindow() {
            return new DynalignRefoldWindow();
        }
    }
}
