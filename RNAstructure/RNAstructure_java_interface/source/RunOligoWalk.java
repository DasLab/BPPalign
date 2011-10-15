/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.util.ArrayList;

/**
 * A class that runs an OligoWalk calculation.
 * @author Jessica Reuter
 */
public class RunOligoWalk extends StrandStarter {
    /**
     * Run an OligoWalk calculation.
     */
    public void execute() {
        try {
	    // Create variables
            int error = 0;
	    StrandDataHolder holder =
		((InputWindow)RNAstructure.getCurrentWindow())
		.getDataHolder();
            ArrayList<Object> data = holder.getData();

            Oligowalk_object walk =
		((OligoWalkWindow)RNAstructure.getCurrentWindow())
		.getOligoObject();

	    int length = (Integer)data.get( 0 );
	    boolean isDNA = (Boolean)data.get( 1 );
	    int option = (Integer)data.get( 2 );
	    double concentration = (Double)data.get( 3 );
	    int sub = (Integer)data.get( 4 );
	    int start = (Integer)data.get( 5 );
	    int stop = (Integer)data.get( 6 );
	    String reportFile = data.get( 7 ).toString();
	    String concString = data.get( 8 ).toString();

	    // Run OligoWalk
	    error = walk.Oligowalk( length, isDNA, option, concentration,
				    sub, start, stop );

	    // Write the report file
	    if( error == 0 ) {
		error = walk.WriteReport( reportFile, length, isDNA, option,
					  concentration, sub, start, stop );
	    }

	    if( error == 0 ) {
		RNAstructure.getCurrentWindow().dispose();

		// Open the results window. Do sizing and positioning manually
		// since this is the only place this window is ever used.
		OligoWalkResultsWindow results =
		    new OligoWalkResultsWindow( walk, length, concString,
						"OligoWalk Results: " +
						reportFile );
		results.pack();
		results.setLocationRelativeTo( RNAstructure.getFrame() );
		results.setVisible( true );
	    } else {
		RNAstructureInfoDialog.error( walk.GetErrorMessage( error ) );
	    }
        } catch( Exception ex ) {
	    ex.printStackTrace();
            RNAstructureInfoDialog.error( 
                "Undefined error during OligoWalk calculation." );
	}
    }
}
