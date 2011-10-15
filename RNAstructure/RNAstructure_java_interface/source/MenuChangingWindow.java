/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A base for all RNAstructure windows that require their own menus.
 * @author Jessica Reuter
 */
public class MenuChangingWindow extends BaseWindow {   
    /**
     * A boolean array that tells whether save (index 0) and print (index 2)
     * buttons should be enabled for this window.
     */
    protected boolean[] enabledArray;

    /**
     * Array of menus that are specific to this particular input window, and
     * disappear when this window does.
     */
    protected SmartMenu[] menus;

    /**
     * Array of integers that tells which sections of the File menu are visible
     * with this window
     */
    protected Integer[] sections;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -278456025264872720L;

    /**
     * Constructor
     */
    protected MenuChangingWindow() {
        super( RNAstructure.getFrame() );

        enabledArray = new boolean[2];
        enabledArray[0] = false;
        enabledArray[1] = false;

        addWindowListener( new WindowManager() );
        RNAstructure.getWindowList().add( this );
    }

    /**
     * Get the permissions of the buttons which are variably enabled.
     * @return  the array of permissions
     */
    public boolean[] getEnabledPermissions() { return enabledArray; }

    /**
     * Set the unique menus that are visible for this window.
     * @param menus  the unique menus that are visible
     */
    public void setMenus( SmartMenu... menus ) { this.menus = menus; }

    /**
     * Set the sections of the file menu that are visible for this window.
     * @param indices  the indices of the file menu that are visible
     */
    public void setSections( Integer... indices ) { sections = indices; }
}
