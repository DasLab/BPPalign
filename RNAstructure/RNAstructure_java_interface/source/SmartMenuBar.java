/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JMenuBar;

/**
 * A menu bar holding hardcoded menus which can enable, disable, or disappear
 * based on the user's needs. Since there is only one instance of this in the
 * RNAstructure GUI, three menus (RNA, DNA, and Help) never change, while the
 * File menu is always in position 0, although changeable.
 * @author Jessica Reuter
 */
public class SmartMenuBar extends JMenuBar {
    /**
     * The menu maker for this bar
     */
    private MenuMaker menuMaker;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -213657577331397211L;

    /**
     * Constructor
     */
    public SmartMenuBar() { menuMaker = new MenuMaker(); }

    /**
     * Disable the variable menus in the bar. Note that the menus File, RNA,
     * DNA, and Help are not affected by this.
     */
    public void disableVariableMenus() {
	int size = getMenuCount() - 1;
	for( int i = 3; i < size; i++ ) { getMenu( i ).setEnabled( false ); }
    }

    /**
     * Enable all menus in this bar.
     */
    public void enableMenus() {
        int size = getMenuCount();
        for( int i = 0; i < size; i++ ) { getMenu( i ).setEnabled( true ); }
    }

    /**
     * Get the menu maker for this bar.
     * @return  the menu maker
     */
    public MenuMaker getMenuMaker() { return menuMaker; }

    /**
     * Set the menu bar back to its initial state, when there are no windows
     * open except the main frame.
     */
    public void refresh() {
        if( getMenuCount() > 0 ) {
            while( !getMenu( 3 ).getText().equals( "Help" ) ) { remove( 3 ); }
            ((SmartMenu)getMenu( 0 )).setSectionsVisible( 1, 3, 4, 5, 7 );
        } else {
            SmartMenu fileMenu = menuMaker.createFileMenu();
            fileMenu.setSectionsVisible( 1, 3, 4, 5, 7 );

            add( fileMenu );
            add( menuMaker.createRNAMenu() );
            add( menuMaker.createDNAMenu() );
            add( menuMaker.createHelpMenu() );
        }

        revalidate();
        repaint();
    }

    /**
     * Remove changeable menus and add new ones.
     * This method is used for menus that are not always present. The following
     * menus should not be included in this menu list: File, RNA, DNA, Help.
     * @param sections  list of sections of the file menu that should be set as
     *                  visible
     * @param menus  the list of changeable menus to add
     */
    public void refresh( Integer[] sections, SmartMenu... menus ) {
        // Remove all changeable menus
        while( !getMenu( 3 ).getText().equals( "Help" ) ) { remove( 3 ); }

        // Add all new menus
        for( SmartMenu menu: menus ) {
            SmartMenu temp = menu;
            temp.setEnabled( false );
            add( temp, getMenuCount() - 1 );
        }

        // Set specific sections of file menu visible
        ((SmartMenu)getMenu( 0 )).setSectionsVisible( sections );

        revalidate();
        repaint();
    }
}
