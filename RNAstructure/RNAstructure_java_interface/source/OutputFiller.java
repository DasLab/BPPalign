/* 
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.io.File;

import javax.swing.JButton;

/**
 * A class that fills a box in a TextInputPanel with a file name, based on
 * parameters from a different box.
 * @author Jessica Reuter
 */
public class OutputFiller extends DataSelector {
    /**
     * The data needed to fill a box.
     */
    private Object[] data;

    /**
     * Constructor
     * @param data  the data needed to fill a box
     */
    private OutputFiller( Object... data ) { this.data = data; }

    /**
     * Create an action that fills in the output file for two strands. This
     * method assumes that both files which create the strands have the same
     * extension and type. This action fills in an output file based on the two
     * previous TextInputPanel indices to the one specified as a parameter for
     * this method.
     * @param index  the index of the button on the panel
     * @param extension  the extension of the output file
     * @return  the multiple strand output filler
     */
    public static OutputFiller fillMultiple( int index, String extension ) {
        return new OutputFiller( index, extension, true );
    }

    /**
     * Create an action that fills in the output file for two strands. This
     * method assumes that both files which create the strands have the same
     * extension and type. This action fills in an output file based on input
     * from two specified indices of a TextInputPanel.
     * @param button  the index of the button on the panel
     * @param idx  the index of the first field on the panel
     * @param idx2  the index of the second field on the panel
     * @param extension  the extension of the output file
     * @return  the multiple strand output filler
     */
    public static OutputFiller fillMultiple( int button, int idx, int idx2,
                                             String extension ) {
        return new OutputFiller( button, idx, idx2, extension );
    }

    /**
     * Create an action that fills in the output file for a single strand
     * @param index  the index of the button on the panel
     * @param extension  the extension of the output file
     * @return  the single strand output filler
     */
    public static OutputFiller fillSingle( int index, String extension ) {
        return new OutputFiller( index, extension, false );
    }

    /**
     * Fill a box with a file name.
     */
    protected void selectInput() {
        // If one file or two is set using one specified index
        if( data.length == 3 ) {
            int index = (Integer)data[0];
            if( !(Boolean)data[2]) {
                fill( index, index - 1, -1, (String)data[1] );
            } else {
                fill( index, index - 2, index - 1, (String)data[1] );
            }
        }

        // If a file is set by two distinct specified indices
        else {
            fill( (Integer)data[0], (Integer)data[1], (Integer)data[2],
                  (String)data[3] );
        }
    }

    /**
     * Fill a file name based on two input names.
     * @param idx  the index to fill
     * @param idx2  the index of the first input name
     * @param idx3  the index of the second input name
     * @param ext  the extension of the new file name
     */
    private void fill( int idx, int idx2, int idx3, String ext ) {
        TextInputPanel panel =
            (TextInputPanel)(((JButton)event.getSource()).getParent());
        String file1 = panel.getFields()[idx2].getText().trim();

        if( !file1.equals( "" ) ) {
            String newFile = null;
            String extOrig =
                file1.substring( file1.lastIndexOf( "." ) );

            if( idx3 != -1 ) {
                String file2 = panel.getFields()[idx3].getText().trim();
                String name1 = new File( file1 ).getName();
                String name2 = new File( file2 ).getName();
                String dir = file1.replace( name1, "" );

                newFile = dir + ( name1.replace( extOrig, "" ) ) + "_" +
                    name2.replace( extOrig, "" ) + "." + ext;
            } else { newFile = file1.replace( extOrig, "." + ext ); }

            // Set the new file name in the proper index
            panel.getFields()[idx].setText( newFile );
        }
    }
}