/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that extends the basic FoldWindow, to fold one strand of nucleic
 * acids into many suboptimal conformations.
 * @author Jessica Reuter
 */
public class FoldSuboptimalWindow extends FoldWindow {
    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = 3035270354322667404L;

    /**
     * Constructor
     * @param acid  the type of nucleic acid this window is folding
     */
    public FoldSuboptimalWindow( String acid ) {
        // Set titles
        super( "Generate All Suboptimal " + acid + " Structures" );
        boolean isRNA = acid.equals( "RNA" );

        // Set the menu bar
        setMenuBar( RNAstructure.getBar().getMenuMaker()
            .createForceSuboptimalMenu() );

        // Create the main panel and the parameter panel
        createMainPanel( 2, "Sequence File", "CT File" );
        inputPanel.disableField( 1 );
        createSuboptimalPanel();

        // Create arrays of cutoff values for suboptimal structure parameter
        // initialization prior to calculation
        Integer[][] percentCutoffs = {
            { 1200, 5 },
            { 800, 8 },
            { 500, 10 },
            { 300, 15 },
            { 120, 20 },
            { 50, 25 },
            { 0, 50 }
        };

        Number[][] absoluteCutoffs = {
            { 1200, 0.25 },
            { 800, 0.5 },
            { 500, 0.75 },
            { 300, 1 },
            { 120, 1.5 },
            { 50, 3 },
            { 0, 10 }
        };

        // Create and set actions for buttons
        // First: Select the sequence file and set its name.
        // Second: Create a strand from the sequence file.
        // Third: Set the name of the default output file.
        // Fourth: Select an output file other than default and set its name.
        // Fifth: Set the default variable structure parameters.
        // "First", "Second", "Third", "Fourth" run in sequence on button 0.
        // "Fifth" executes alone on button 2.
        OpenSaveSetter first =
            OpenSaveSetter.open( 0, "Sequence Files", "seq" );
        StrandCreator second = StrandCreator.createSEQStrand( 0, isRNA );
        OutputFiller third = OutputFiller.fillSingle( 1, "ct" );
        FoldingParamFiller fourth =
            FoldingParamFiller.setVariable( percentCutoffs, absoluteCutoffs );
        OpenSaveSetter fifth = OpenSaveSetter.save( 1, "CT Files", "ct" );
        ChainedAction foldSequence =
            new ChainedAction( first, second, third, fourth );

        inputPanel.setAction( 0, foldSequence );
        inputPanel.setAction( 1, fifth );

        // Set components in place and designate the start button action
        setComponents( new GenerateSuboptimalStructures() );
    }

    /**
     * Get the data for generation of suboptimal structures.
     * @return  an array of suboptimal structure object data
     */
    public Object[] getData() {
        Float percent =
            Float.parseFloat( inputPanel2.getFields()[0].getText() );
        Double g = Double.parseDouble( inputPanel2.getFields()[1].getText() );

        Object[] data = { inputPanel.getFields()[1].getText(), percent, g };
        return data;
    }

    /**
     * An inner class which allows a FoldSuboptimalWindow to be constructed
     * from the more generic context of the RNAWindowFactory.
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
         * Create a FoldSuboptimalWindow.
         * @return  a new FoldSuboptimalWindow
         */
        public FoldSuboptimalWindow createWindow() {
            return new FoldSuboptimalWindow( acid );
        }
    }
}