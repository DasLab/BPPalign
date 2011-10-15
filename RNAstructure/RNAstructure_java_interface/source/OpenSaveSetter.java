/* 
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JButton;

/**
 * A class that handles entry of a file to open or save in a TextInputPanel
 * @author Jessica Reuter
 */
public class OpenSaveSetter extends DataSelector {
    /**
     * The data needed to create a setter
     */
    private Object[] data;

    /**
     * Constructor
     * @param data  the data needed to create a setter 
     */
    private OpenSaveSetter( Object... data ) { this.data = data; }

    /**
     * Create an action that allows for selecting a file name to open
     * @param index  the index of the button in the panel
     * @param desc  the file filter description
     * @param ext  the file filter extension
     * @return  the open action listener
     */
    public static OpenSaveSetter open( int index, String desc, String ext ) {
        return new OpenSaveSetter( index, "Open", desc, ext, false );
    }

    /**
     * Create an action that allows for selecting a file name to save
     * @param index  the index of the button in the panel
     * @param desc  the file filter description
     * @param ext  the file filter extension
     * @return  the save action listener
     */
    public static OpenSaveSetter save( int index, String desc, String ext ) {
        return new OpenSaveSetter( index, "Save", desc, ext, true );
    }

    /**
     * Select a file, either to open or save
     */
    protected void selectInput() {
        // Get data
        int index = (Integer)data[0];
        boolean allowed = (Boolean)data[4];
        String type = data[1].toString();
        String description = data[2].toString();
        String extension = data[3].toString();

        boolean ok = true;
        TextInputPanel panel =
            (TextInputPanel)(((JButton)event.getSource()).getParent());

        // Check to make sure all previous selections are present
        for( int i = 0; i < index; i++ ) {
            String next = panel.getFields()[i].getText();
            if( next.trim().equals( "" ) ) {
                String message = "This file cannot be selected " +
                    "at this time.";
                RNAstructureInfoDialog.warning( message );
                ok = false;
                break;
            }
        }

        // Select a file if necessary
        if( ok ) {
            FileSelector selector = new FileSelector( allowed, type,
                description + "," + extension );
            String newFile = selector.showSelector();

            if( !newFile.equals( "" ) ) {
                panel.getFields()[index].setText( newFile );
            }
        }
    }
}