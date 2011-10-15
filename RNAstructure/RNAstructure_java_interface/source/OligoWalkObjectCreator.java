/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that creates an Oligowalk_object from various contexts.
 * @author Jessica Reuter
 */
public class OligoWalkObjectCreator extends ActionBase {
    /**
     * Boolean telling whether an object should be RNA (true) or DNA (false)
     */
    private boolean isRNA;

    /**
     * Constructor
     * @param type  true = RNA, false = DNA
     */
    public OligoWalkObjectCreator( boolean type ) { isRNA = type; }

    /**
     * Create an Oligowalk_object
     */
    protected void act() {
	// Get file name
	OligoWindow window = (OligoWindow)RNAstructure.getCurrentWindow();
	String file = window.getFileInputPanel().getFields()[0].getText();
    
	// Create the Oligowalk_object if the file is OK
	if( !file.equals( "" ) ) {
	    Oligowalk_object object =
		( window instanceof OligoScreenWindow ) ?
		    new Oligowalk_object( isRNA ) :
		    new Oligowalk_object( file, 1 );

	    int code = object.GetErrorCode();
	    if( code == 0 ) {
		object.setIsrna( isRNA );
		window.setOligoObject( object );
		RNAstructure.getBar().enableMenus();
	    } else {
		RNAstructureInfoDialog.error( object.GetErrorMessage( code ) );
	    }
	}
    }
}
