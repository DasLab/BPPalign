/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JTextField;

/**
 * A class that runs a Dynalign refolding calculation.
 * @author Jessica Reuter
 */
public class RefoldDynalignRunner extends StrandStarter {
    /**
     * Refold a Dynalign calculation.
     */
    protected void execute() {
        int error = 0;
	Dynalign_object dynalign = null;
	DynalignRefoldWindow refold =
	    (DynalignRefoldWindow)RNAstructure.getCurrentWindow();

        try {
            // Draw in variables
            JTextField[] fields1 = refold.getFilePanel().getFields();
            JTextField[] fields2 = refold.getParameterPanel().getFields();

            String saveFile = fields1[0].getText();
            String ctFile1 = fields1[1].getText();
            String ctFile2 = fields1[2].getText();
            String alignFile = fields1[3].getText();

            short percent = Short.parseShort( fields2[0].getText() );
            short structures = Short.parseShort( fields2[1].getText() );
            short structWin = Short.parseShort( fields2[2].getText() );
            short alignWin = Short.parseShort( fields2[3].getText() );

            // Refold an existing save file.
	    dynalign = new Dynalign_object( saveFile, percent, structures,
                                            structWin, alignWin );
	    error = dynalign.GetErrorCode();
	    if( error != 0 ) { throw new Exception(); }

            // Write output
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
	    confirmed.act();
        } catch( Exception ex ) {
            String message = ( error != 0 ) ?
                dynalign.GetErrorMessage( error ) :
                "Undefined error during Dynalign calculation.";
            RNAstructureInfoDialog.error( message );
	    ex.printStackTrace();
        }
    }
}
