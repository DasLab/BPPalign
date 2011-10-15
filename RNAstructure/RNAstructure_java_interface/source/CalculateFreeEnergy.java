/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
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
 * A class which calculates the free energy of a strand of nucleic acids.
 * @author Jessica Reuter
 */
public class CalculateFreeEnergy extends StrandStarter {
    /**
     * Calculate free energy.
     */
    public void execute() {
        try {
            boolean error = false;

	    StrandDataHolder holder =
		((InputWindow)RNAstructure.getCurrentWindow()).getDataHolder();
            ArrayList<Object> data = holder.getData();
            RNA strand = holder.getStrand();

            // If the user wants to write a thermodynamic details file, write
            // the complex file.
            if( (Boolean)data.get( 2 ) == true ) {
                String ctFile = data.get( 0 ).toString();

                // Write the thermodynamic details file
                String outfile = ctFile.replace( ".ct", ".out" );
                int code = strand.WriteThermodynamicDetails( outfile );
                error = holder.isErrorStatus( code );

                // Show an error if a problem occurred
                if( error ) {
                    RNAstructureInfoDialog.error( 
                        "Error writing thermodynamic details file." );
                }
            }

            // If the user doesn't want a thermodynamic details file, write
            // simple output instead.
            else {
                try {
                    // Get the number of structures to write data for
                    int structures = strand.GetStructureNumber();

                    // Initialize the file writer
                    String out = data.get( 1 ).toString();
                    BufferedWriter writer =
                        new BufferedWriter( new FileWriter( out ) );

                    // For each structure, calculate free energy
                    for( int i = 1; i <= structures; i++ ) {
                        double energy = strand.CalculateFreeEnergy( i );
                        error = holder.isErrorStatus();

                        // If no error occurred, write result to output file
                        // Otherwise, stop writing output.
                        if( !error ) {
                            writer.write( "Structure: " + 
                                Integer.toString( i ) + "   Energy = " +
                                Double.toString( energy ) );

                            writer.newLine();
                            writer.flush();
                        } else break;
                    }
                    writer.close();
                } catch( IOException e ) {
                    RNAstructureInfoDialog.error(
                        "Error writing free energy file." );
                }

            }
        } catch( Exception ex ) {
            RNAstructureInfoDialog.error( 
                "Undefined error during free energy calculation." );
        } finally {
	    RNAstructure.getCurrentWindow().dispose();
            RNAstructureInfoDialog.message( "Efn2 Complete." );
        }
    }
}
