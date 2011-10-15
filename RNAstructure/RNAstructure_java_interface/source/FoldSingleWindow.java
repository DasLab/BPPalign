/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that extends a basic FoldWindow, to fold a single strand of nucleic
 * acids into its optimal conformation.
 * @author Jessica Reuter
 */
public class FoldSingleWindow extends FoldWindow {
    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -8160909075883113480L;

    /**
     * Constructor
     * @param acid  the type of nucleic acid this window is folding
     */
    public FoldSingleWindow( String acid ) {
        // Set window defaults
        super( acid + " Single Strand Fold" );
        boolean isRNA = acid.equals( "RNA" );

        // Set the menu bar
        setMenuBar(
            RNAstructure.getBar().getMenuMaker().createForceSingleMenu() );

        // Create the main panel
        createMainPanel( 2, "Sequence File", "CT File" );
        inputPanel.disableField( 1 );

        // Create the folding option box and parameter panel
        createBox();
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
        // First: Select the sequence file and set its name.
        // Second: Create a strand from the sequence file.
        // Third: Set the name of the default output file.
        // Fourth: Set the default variable structure parameters.
        // Fifth: Select an output file other than default and set its name.
        // "First", "Second", "Third", "Fourth" run in sequence on button 0.
        // "Fifth" executes alone on button 1.
        OpenSaveSetter first =
            OpenSaveSetter.open( 0, "Sequence Files", "seq" );
        StrandCreator second = StrandCreator.createSEQStrand( 0, isRNA );
        OutputFiller third = OutputFiller.fillSingle( 1, "ct" );
        FoldingParamFiller fourth =
            FoldingParamFiller.setVariable( "10", "20", cutoffs );
        OpenSaveSetter fifth = OpenSaveSetter.save( 1, "CT Files", "ct" );
        ChainedAction foldSequence =
            new ChainedAction( first, second, third, fourth );

        inputPanel.setAction( 0, foldSequence );
        inputPanel.setAction( 1, fifth );

        // Set components in place and designate the start button action
        setComponents( new FoldNucleicAcids( "Single Stranded Folding" ) );
    }

    /**
     * Set the files in the main input panel.
     * @param sequence  the sequence file name to set
     * @param ct  the ct file name to set
     */
    public void setFileNames( String sequence, String ct ) {
        inputPanel.getFields()[0].setText( sequence );
        inputPanel.getFields()[1].setText( ct );
    }

    /**
     * An inner class which allows a FoldSingleWindow to be constructed from
     * the more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * The type of nucleic acid to fold
         */
        private String acid;

        /**
         * Constructor
         * @param type  the type of nucleic acid to fold
         */
        public Factory( String type ) { acid = type; }

        /**
         * Create a FoldSingleWindow.
         * @return  a new FoldSingleWindow
         */
        public FoldSingleWindow createWindow() {
            return new FoldSingleWindow( acid );
        }
    }
}