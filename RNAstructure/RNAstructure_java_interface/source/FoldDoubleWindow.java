/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that extends the basic FoldWindow, to fold two strands of nucleic
 * acids into their optimal conformation.
 * @author Jessica Reuter
 */
public class FoldDoubleWindow extends FoldWindow {
    /**
     * A serial ID, required for extending GUI classes
     */
    private static final long serialVersionUID = 845745453839396944L;

    /**
     * Constructor
     * @param acid  the type of nucleic acid this window is folding
     */
    public FoldDoubleWindow( String acid ) {
        // Set titles
        super( acid + " Bimolecular Fold" );
        boolean isRNA = acid.equals( "RNA" );

        // Set the menu bar
        setMenuBar(
            RNAstructure.getBar().getMenuMaker().createForceDoubleMenu() );

        // Create the main panel
        createMainPanel( 3, "Sequence File 1", "Sequence File 2", "CT File" );
        inputPanel.disableField( 1 );
        inputPanel.disableField( 2 );

        // Create the folding option box and parameter panel
        createBox();
        createFoldingPanel();

        // Create and set actions for buttons
        // First: Select the first sequence file and set its name.
        // Second: Select the second sequence file and set its name.
        // Third: Create a strand from both sequence files.
        // Fourth: Set the name of the default output file.
        // Fifth: Set the default constant structure parameters.
        // Sixth: Select an output file other than default and set its name.
        // "First" executes alone on button 0.
        // "Second", "Third", "Fourth", "Fifth" run in sequence on button 1.
        // "Sixth" executes alone on button 2.
        OpenSaveSetter first =
            OpenSaveSetter.open( 0, "Sequence Files", "seq" );
        OpenSaveSetter second =
            OpenSaveSetter.open( 1, "Sequence Files", "seq" );
        StrandCreator third = StrandCreator.createHybridSEQStrand( 1, isRNA );
        OutputFiller fourth = OutputFiller.fillMultiple( 2, "ct" );
        FoldingParamFiller fifth =
            FoldingParamFiller.setConstant( "50", "20", "0" );
        OpenSaveSetter sixth = OpenSaveSetter.save( 2, "CT Files", "ct" );
        ChainedAction bifoldSequence =
            new ChainedAction( second, third, fourth, fifth );

        inputPanel.setAction( 0, first );
        inputPanel.setAction( 1, bifoldSequence );
        inputPanel.setAction( 2, sixth );

        // Set components in place and designate the start button action
        setComponents( new FoldNucleicAcids( "Bimolecular Folding" ) );
    }

    /**
     * An inner class which allows a FoldDoubleWindow to be constructed from
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
         * Create a FoldDoubleWindow.
         * @return  a new FoldDoubleWindow
         */
        public FoldDoubleWindow createWindow() {
            return new FoldDoubleWindow( acid );
        }
    }
}