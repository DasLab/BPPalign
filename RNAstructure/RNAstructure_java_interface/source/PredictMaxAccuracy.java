/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.util.ArrayList;

/**
 * A class that handles prediction of maximum expected accuracy structures.
 * @author Jessica Reuter
 */
public class PredictMaxAccuracy extends StrandStarter {
    /**
     * Predict maximum expected accuracy structures.
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

            // Get folding parameters
	    float percent = ((Double)data.get( 1 )).floatValue();
	    int maxStructures = (Integer)data.get( 2 );
	    int windowSize = (Integer)data.get( 3 );
	    double gamma = ((Float)data.get( 4 )).doubleValue();

	    // Fold the strand(s) with appropriate parameters
	    code = strand.MaximizeExpectedAccuracy( percent, maxStructures,
						    windowSize, gamma );
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
            RNAstructureInfoDialog.error( "<html><center>" +
	        "Undefined Error During Maximum<br/>" +
		"Expected Accuracy Structure Generation." );
        }
    }
}
