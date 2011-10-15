/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.util.ArrayList;

/**
 * A class that handles the calculation of stochastic probability sampling.
 * @author Jessica Reuter
 */
public class SampleStochastic extends StrandStarter {
    /**
     * Calculate stochastic probabilities
     */
    public void execute() {
        try {
	    StrandDataHolder holder =
		((InputWindow)RNAstructure.getCurrentWindow())
		.getDataHolder();
	    ArrayList<Object> data = holder.getData();
	    RNA strand = holder.getStrand();
	    boolean error = false;
		
	    // Sample structures
	    int structures = (Integer)data.get( 1 );
	    int seed = (Integer)data.get( 2 );
	    String outputFile = data.get( 0 ).toString();

	    int code = strand.Stochastic( structures, seed );
	    error = holder.isErrorStatus( code );

	    // Write the CT File
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
            RNAstructureInfoDialog.error(
                "Undefined Error During Stochastic Sampling." );
        }
    }
}
