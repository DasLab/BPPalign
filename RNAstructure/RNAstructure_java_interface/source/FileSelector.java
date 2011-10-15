/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.io.File;

import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;

/**
 * A class that handles file selectors for the RNAstructure GUI.
 * @author Jessica Reuter
 */
public class FileSelector extends JFileChooser {
    /**
     * The last directory a file chooser selected a file from
     */
    private static String directory = System.getProperty( "user.home" );

    /**
     * The file selected from this selector.
     */
    private String file;

    /**
     * A serialized ID, required when extending GUI classes
     */
    private static final long serialVersionUID = -6893507606169775466L;

    /**
     * The type of file selector this is
     */
    private String type;

    /**
     * Constructor
     * 
     * Each element in the list of file name extensions has two parts, the
     * description of the extension type, followed by the extension itself,
     * with no dot preceeding the extension. The two parts of the extension are
     * separated by a comma; neither the description nor the extension can
     * contain a comma.
     * 
     * @param all  true if all files are allowed in this chooser, false if not
     * @param type  the type of file chooser to open
     * @param list  the list to be used for file name extensions
     */
    public FileSelector( boolean all, String type, String... list ) {
        this.type = type;
        file = "";

        // Set up the file chooser and its general file filters
        setAcceptAllFileFilterUsed( all );
        setCurrentDirectory( new File( directory ) );
        setFileSelectionMode( JFileChooser.FILES_ONLY );

        // Add specific file filters by extensions
        for( String element: list ) {
            String[] data = element.split( "," );
            String description = data[0] + " (*." + data[1] + ")";
            addChoosableFileFilter(
                new FileNameExtensionFilter( description, data[1] ) );
        }
    }
    
    /**
     * Show the file selector.
     * @return  the selected file.
     */
    public String showSelector() {
        int value;
        String safeType = type.toLowerCase();

        if( safeType.equals( "open" ) ) { value = showOpenDialog( null ); }
        else if( safeType.equals( "save" ) ) {
           value = showSaveDialog( null );
        } else {
            safeType = safeType.replaceFirst(
                Character.toString( safeType.charAt( 0 ) ),
                Character.toString( safeType.charAt( 0 ) ).toUpperCase() );
            value = showDialog( null, safeType );
        }

        // Get the selected file, if there is one
        if( value == JFileChooser.APPROVE_OPTION ) {
            file = getSelectedFile().getAbsolutePath();
            String name = getSelectedFile().getName();
            directory = file.replace( name, "" );

            String description = getFileFilter().getDescription();
            int start = description.lastIndexOf( "." );
            int end = description.lastIndexOf( ")" );
            String extension = description.substring( start, end );

            if( !file.endsWith( extension ) ) {
                file = file.concat( extension );
            }
        }

        return file;
    }
}