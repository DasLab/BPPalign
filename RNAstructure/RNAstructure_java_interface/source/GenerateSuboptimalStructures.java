/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.util.ArrayList;

/**
 * A class that handles the action of generating suboptimal structures from a
 * strand of nucleic acids.
 * @author Jessica Reuter
 */
public class GenerateSuboptimalStructures extends StrandStarter {
    /**
     * Generate suboptimal structures.
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
            String outputFile = data.get( 0 ).toString();

            // Generate structures
            Float percent = ((Double)data.get( 1 )).floatValue();
            Double deltaG = (Double)data.get( 2 );
            code = strand.GenerateAllSuboptimalStructures( percent, deltaG );
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
        } catch( Exception e ) {
            String suboptimalUndefined = 
	        "Undefined Error During Generation of Suboptimal Structures.";
            RNAstructureInfoDialog.error( suboptimalUndefined );
        }
    }
}
