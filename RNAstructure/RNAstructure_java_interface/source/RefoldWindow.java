/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that extends the basic FoldWindow, to refold a strand or strands of
 * nucleic acids generated from a folding save file.
 * @author Jessica Reuter
 */
public class RefoldWindow extends FoldWindow {
    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -5104141132034854086L;

    /**
     * Constructor
     */
    public RefoldWindow() {
        // Set window defaults
        super( "Refold From Save File" );

        // Set the menu bar
        setSections( 1, 3, 4, 5, 7 );
        SmartMenuBar bar = RNAstructure.getBar();
        setMenus( bar.getMenuMaker().createForceRefoldMenu() );
        bar.refresh( sections, menus );

        // Create the main panel and parameter panel
        createMainPanel( 2, "Save File", "CT File" );
        inputPanel.disableField( 1 );
        createFoldingPanel();

        // Create the array of cutoff values for suboptimal structure parameter
        // initialization prior to calculation
        Integer[][] cutoffs = {
            { 1200, 20 },
            { 800, 15 },
            { 500, 11 },
            { 300, 7 },
            { 120, 5 },
            { 50, 3 },
            { 0, 2 }
        };

        // Create and set actions for buttons
        // First: Select the save file and set its name.
        // Second: Create a strand from the save file. RNA can be hardcoded
        //         because the save file takes care of the nucleic acid type.
        // Third: Set the name of the default output file.
        // Fourth: Set the default variable structure parameters.
        // Fifth: Select an output file other than default and set its name.
        // "First", "Second", "Third", "Fourth" runs in sequence on button 0.
        // "Fifth" executes alone on button 1.
        OpenSaveSetter first =
            OpenSaveSetter.open( 0, "Folding Save Files", "sav" );
        StrandCreator second = StrandCreator.createSAVStrand( 0, true );
        OutputFiller third = OutputFiller.fillSingle( 1, "ct" );
        FoldingParamFiller fourth =
            FoldingParamFiller.setVariable( "10", "20", cutoffs );
        OpenSaveSetter fifth = OpenSaveSetter.save( 1, "CT Files", "ct" );
        ChainedAction foldSequence =
            new ChainedAction( first, second, third, fourth );

        inputPanel.setAction( 0, foldSequence );
        inputPanel.setAction( 1, fifth );

        // Set components in place and designate the start button action
        setComponents( new RefoldNucleicAcids() );
    }

    /**
     * An inner class which allows a RefoldWindow to be constructed from the
     * more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * Create a RefoldWindow.
         * @return  a new RefoldWindow
         */
        public RefoldWindow createWindow() {
            return new RefoldWindow();
        }
    }
}