/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.util.ArrayList;

/**
 * A class that runs an OligoScreen calculation.
 * @author Jessica Reuter
 */
public class RunOligoScreen extends StrandStarter {
    /**
     * Run an OligoScreen calculation.
     */
    public void execute() {
        try {
            int errorCode = 0;
	    StrandDataHolder holder =
		((InputWindow)RNAstructure.getCurrentWindow())
		.getDataHolder();
            ArrayList<Object> data = holder.getData();

            Oligowalk_object screen =
                ((OligoScreenWindow)RNAstructure.getCurrentWindow())
		.getOligoObject();
            String inFile = data.get( 0 ).toString();
            String outFile = data.get( 1 ).toString();

            if( !holder.isErrorStatus( errorCode ) ) {
                errorCode = screen.OligoScreen( inFile, outFile );
                holder.isErrorStatus( errorCode );
            }
        } catch( Exception ex ) {
            RNAstructureInfoDialog.error( 
                "Undefined error during OligoScreen calculation." );
            ex.printStackTrace();
        } finally {
            RNAstructureInfoDialog.message( "OligoScreen Complete." );
            RNAstructure.getCurrentWindow().dispose();
        }
    }
}
