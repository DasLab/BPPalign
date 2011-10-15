/* 
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JButton;

/**
 * A class that creates a Dynalign_object in multiple contexts.
 * @author Jessica Reuter
 */
public class DynalignObjectCreator extends ActionBase {
    /**
     * The array of data used to create the Dynalign_object (may or may not be
     * used for a specific creation)
     */
    private Object[] data;

    /**
     * The file used to create the Dynalign_object (may or may not be used for
     * a specific creation)
     */
    private String file;
    
    /**
     * The index used to create the Dynalign_object (may or may not be used for
     * a specific creation)
     */
    private int index;

    /**
     * File constructor
     * @param file  the file used to create the Dynalign_object
     */
    public DynalignObjectCreator( String file ) {
        this.file = file;
        index = -1;
    }

    /**
     * Index constructor
     * @param index  the index used to create the Dynalign_object
     */
    public DynalignObjectCreator( int index ) {
        file = "";
        this.index = index;
    }

    /**
     * Multiple argument constructor
     * @param idx1  the index of the first file
     * @param idx2  the index of the second file
     * @param type  the type of objects to use in creation
     * @param acid  the backbone type
     */
    public DynalignObjectCreator( int idx1, int idx2, int type, String acid ) {
        data = new Object[4];
        data[0] = idx1;
        data[1] = idx2;
        data[2] = type;
        data[3] = acid.equals( "RNA" );
    }

    /**
     * Create a Dynalign_object.
     */
    protected void act() {
        Dynalign_object obj = null;
        TextInputPanel panel =
            (TextInputPanel)(((JButton)event.getSource()).getParent());
    
        // If a Dynalign object for folding is created, create a new object
        if( file == null ) {
            int index1 = (Integer)data[0];
            int index2 = (Integer)data[1];
            int type = (Integer)data[2];
            boolean isRNA = (Boolean)data[3];

            String file1 = panel.getFields()[index1].getText();
            String file2 = panel.getFields()[index2].getText();
            obj = new Dynalign_object( file1, type, file2, type, isRNA );
        }

        // Otherwise, if a save file is being read, create the "old" object
        else if( !file.equals( "" ) ) {
            if( index != -1 ) { file = panel.getFields()[index].getText(); }
            obj = new Dynalign_object( file );
        }

        // Check and set the object
        if( obj != null ) {
            int code = obj.GetErrorCode();

            if( code == 0 ) {
		DynalignWindow window =
		    (DynalignWindow)RNAstructure.getCurrentWindow();
                window.setDynalign( obj );
                window.setWorkingStrand( obj.GetRNA1() );
                RNAstructure.getBar().enableMenus();
            } else {
                RNAstructureInfoDialog.error( obj.GetErrorMessage( code ) );
            }
        }
    }
}
