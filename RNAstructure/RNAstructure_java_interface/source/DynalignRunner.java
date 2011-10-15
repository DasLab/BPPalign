/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JTextField;

/**
 * A class that runs a Dynalign calculation.
 * @author Jessica Reuter
 */
public class DynalignRunner extends StrandStarter {
    /**
     * Run the Dynalign calculation.
     */
    protected void execute() {
        int error = 0;
        Dynalign_object dynalign = null;

        try {
            // Draw in variables
            DynalignFoldWindow input =
		(DynalignFoldWindow)RNAstructure.getCurrentWindow();
            dynalign = input.getDynalign();

            JTextField[] fields1 = input.getFilePanel().getFields();
            JTextField[] fields2 = input.getParameterPanel().getFields();
            boolean[] checkedOptions = input.getBoxSelections();

            String ctFile1 = fields1[1].getText();
            String ctFile2 = fields1[3].getText();
            String alignFile = fields1[4].getText();

            short percent = Short.parseShort( fields2[0].getText() );
            short structures = Short.parseShort( fields2[1].getText() );
            short structWindow = Short.parseShort( fields2[2].getText() );
            short alignWindow = Short.parseShort( fields2[3].getText() );

            float gap = input.getGapPenalty();
            boolean inserts = checkedOptions[0];
            boolean createSave = checkedOptions[1];

            String saveFile = "";
	    if( createSave ) {
		saveFile = alignFile.substring( 0, alignFile.indexOf( "." ) );
		saveFile = saveFile.concat( ".dsv" );
	    }

            // Run Dynalign and write output
            error = dynalign.Dynalign( structures, structWindow, alignWindow,
                percent, (short)-99, gap, inserts, saveFile );
            if( error != 0 ) { throw new Exception(); }

            error = dynalign.GetRNA1().WriteCt( ctFile1 );
            if( error != 0 ) { throw new Exception(); }

            error = dynalign.GetRNA2().WriteCt( ctFile2 );
            if( error != 0 ) { throw new Exception(); }

            dynalign.WriteAlignment( alignFile );
            if( ( error = dynalign.GetErrorCode() ) != 0 ) {
                throw new Exception();
            }

	    // Close the calculation window and open a sketcher if the user
            // wants to
	    RNAstructure.getCurrentWindow().dispose();

	    ChainedAction confirmed =
		new ChainedAction(
		    new CloseAction(),
		    new Imager.TwoImagersFactory( ctFile1, ctFile2 ) );

	    RNAstructureInfoDialog.confirm(
		"Do you want to draw structures?", confirmed );
        } catch( Exception ex ) {
            String message = ( error != 0 ) ?
                dynalign.GetErrorMessage( error ) :
                "Undefined error during Dynalign calculation.";
            RNAstructureInfoDialog.error( message );
        }
    }
}
