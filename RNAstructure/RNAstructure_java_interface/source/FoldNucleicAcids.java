/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.util.ArrayList;

/**
 * A class that handles the action of folding one or more strands of nucleic
 * acids into their optimal conformations. The class can be used for one or two
 * strands and shouldn't be used for suboptimal structure folding or refolding.
 * @author Jessica Reuter
 */
public class FoldNucleicAcids extends StrandStarter {
    /**
     * The type of folding being done
     */
    private String type;

    /**
     * Constructor
     * @param type  the type of folding being done
     */
    public FoldNucleicAcids( String type ) { this.type = type; }

    /**
     * Fold one or more strands of nucleic acids.
     */
    public void execute() {
        try {
            boolean error = false;
            int code = 0;

	    StrandDataHolder holder = 
		((InputWindow)RNAstructure.getCurrentWindow())
		.getDataHolder();
            ArrayList<Object> data = holder.getData();
            RNA strand = holder.getStrand();

            // Get the name of the output file
            String outputFile = data.get( 0 ).toString();

            // Read in SHAPE constraints, if they exist
            int length = data.size();
            if( length > 6 ) {
		FoldWindow copy = (FoldWindow)RNAstructure.getCurrentWindow();
		error = !copy.readSHAPEConstraints( 6, length );
	    }

            // Fold the strand(s)
            if( !error ) {
                // Get folding parameters and the choice to create a save file
                float percent = ((Double)data.get( 1 ) ).floatValue();
                int maxStructures = ((Double)data.get( 2 ) ).intValue();
                int windowSize = ((Double)data.get( 3 ) ).intValue();
                int maxLoopSize = (Integer)data.get( 4 );

                // Get the name of a putative save file
                Boolean createSave = (Boolean)data.get( 5 );
                String saveFile = ( createSave ) ? 
                    outputFile.replaceAll( "\\.ct", "\\.sav" ) : "";

                // Fold the strand(s) with appropriate parameters
                if( type.equals( "Single Stranded Folding" ) ) {
                    code = strand.FoldSingleStrand( percent, maxStructures,
                        windowSize, saveFile, maxLoopSize );
                } else {
                    HybridRNA hybrid = (HybridRNA)strand;
                    hybrid.SetForbidIntramolecular(
                        !hybrid.GetForbidIntramolecular() );

                    code = hybrid.FoldBimolecular( percent, maxStructures,
                        windowSize, saveFile, maxLoopSize );
                }
                error = holder.isErrorStatus( code );
            }

            // Write the CT file
            if( !error ) {
                code = strand.WriteCt( outputFile );
                error = holder.isErrorStatus( code );
            }

            // Draw structures if the user asks to
            if( !error ) {
                RNAstructure.getCurrentWindow().dispose();

                ChainedAction confirmed =
                    new ChainedAction(
                        new CloseAction(),
                        new WindowDisplayer(
                            new Imager.StructureFactory(
                                new DrawStructure( outputFile ) ) ) );

                RNAstructureInfoDialog.confirm(
                    "Do you want to draw structures?", confirmed );
            }
        } catch( Exception ex ) {
            RNAstructureInfoDialog.error( "Undefined Error During " +
                type + "." );
        }
    }
}
