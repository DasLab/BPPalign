/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.util.ArrayList;

/**
 * A class that refolds nucleic acids from folding save files. This class can
 * be used for both single-stranded and bimolecular refolding.
 * @author Jessica Reuter
 */
public class RefoldNucleicAcids extends StrandStarter {
    /**
     * A serialized ID, required for extending an action
     */
    private static final long serialVersionUID = 8473399038737617831L;

    /**
     * Refold nucleic acids from an existing save file.
     */
    public void execute() {
        try {
            boolean error = false;

	    StrandDataHolder holder =
		((InputWindow)RNAstructure.getCurrentWindow())
		.getDataHolder();
            ArrayList<Object> data = holder.getData();
            RNA strand = holder.getStrand();

            // Get output file and folding parameters
            String outputFile = data.get( 0 ).toString();
            float percent = ((Double)data.get( 1 ) ).floatValue();
            int maxStructures = ((Double)data.get( 2 ) ).intValue();
            int windowSize = ((Double)data.get( 3 ) ).intValue();

            // Refold the strand
            int code = strand.ReFoldSingleStrand( percent, maxStructures,
                windowSize );
            error = holder.isErrorStatus( code );

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
            RNAstructureInfoDialog.error(
                "Undefined Error During Refolding." );
        }
    }
}
