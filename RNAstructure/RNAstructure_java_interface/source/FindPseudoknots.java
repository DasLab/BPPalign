/*
 * (c) 2010  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

/**
 * A class responsible for finding pseudoknots in a strand of nucleic acids.
 * @author Jessica Reuter
 */
public class FindPseudoknots extends StrandStarter {
    /**
     * Identify pseudoknots.
     */
    public void execute() {
        try {
            // Get variables
            int error = 0;
	    StrandDataHolder holder = 
		((InputWindow)RNAstructure.getCurrentWindow())
		.getDataHolder();
            ArrayList<Object> data = holder.getData();
            RNA strand = holder.getStrand();

	    Integer iterations = (Integer)data.get( 1 );
	    Integer minHelixLength = (Integer)data.get( 2 );

            // Break pseudoknots
            error = strand.ProbKnot( iterations, minHelixLength );

            // Write the resultant CT file
            String ctFile = null;
            if( !holder.isErrorStatus( error ) ) {
                ctFile = data.get( 0 ).toString();
                error = strand.WriteCt( ctFile );
                holder.isErrorStatus( error );
            }

            // Draw structures if the user asks to
            if( error == 0 ) {
                RNAstructure.getCurrentWindow().dispose();

                ChainedAction confirmed =
                    new ChainedAction(
                        new CloseAction(),
			new WindowDisplayer(
                            new Imager.StructureFactory(
                                new DrawStructureCircular( ctFile ) ) ) );

                RNAstructureInfoDialog.confirm(
                    "Do you want to draw structures?", confirmed );
            }
        } catch( Exception ex ) {
            RNAstructureInfoDialog.error( 
                "Undefined error during pseudoknot breakage." );
        }
    }
}
