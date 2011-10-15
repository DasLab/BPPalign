/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.Font;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JCheckBoxMenuItem;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JSeparator;
import javax.swing.KeyStroke;

/**
 * A menu which is hardcoded for all of its actions and items, but can show
 * different combinations of its items when the user requests
 * @author Jessica Reuter
 */
public class SmartMenu extends JMenu {
    /**
     * A serialized ID, required when extending GUI classes
     */
    private static final long serialVersionUID = -5973688271178011972L;

    /**
     * A list of the starting points for the different sections within a menu
     */
    private List<Integer> sections = new ArrayList<Integer>();

    /**
     * Constructor
     * @param title  the name of the menu
     */
    public SmartMenu( String title ) {
        super( title );
        setFont( new Font( getFont().getFontName(), 0, 16 ) );
        sections.add( 0 );
    }

    /**
     * Add a section start to the menu.
     * @param index  the index of the section start
     */
    public void addSectionStart( int index ) { sections.add( index ); }

    /**
     * Create a check box menu item.
     * @param txt  the text of the menu item
     * @param act  the action this menu item does
     * @param lbl  the text this item displays in the main label on rollover
     */
    public void createCheckItem( String txt, ActionBase act, String lbl ) {
        add( initializeItem( txt, act, lbl, 2 ) );
    }

    /**
     * Create a menu item.
     * @param text  the text of the menu item
     * @param act  the action this menu item does
     * @param label  the text this item displays in the main label on rollover
     */
    public void createItem( String text, ActionBase act, String label ) {
        add( initializeItem( text, act, label, 1 ) );
    }

    /**
     * Create a menu item with an associated keyboard shortcut.
     * @param text  the text of the menu item
     * @param act  the action this menu item does
     * @param label  the text this item displays in the main label on rollover
     * @param key  the keystroke this item responds to
     */
    public void createItem( String text, ActionBase act, String label,
                            String key ) {
        JMenuItem item = initializeItem( text, act, label, 1 );
        item.setAccelerator( KeyStroke.getKeyStroke( key ) );
        add( item );
    }

    /**
     * Get the number of sections in this menu.
     * @return  the number of sections in this menu.
     */
    public int getNumberOfSections() { return sections.size(); }

    /**
     * Initialize a menu item with text and particular actions.
     * @param text  the text of the menu item
     * @param act  the action this menu item does
     * @param label  the text this item displays in the main label on rollover
     * @param type  the type of menu item this is
     *        1 = simple item, only text
     *        2 = text and a check box
     * @return  the new menu item
     */
    private JMenuItem initializeItem( String text, ActionBase act,
                                      final String label, int type ) {
        // Create the item with its action listener
        JMenuItem item = ( type == 1 ) ?
            new JMenuItem( text ) : new JCheckBoxMenuItem( text );

        item.setFont( new Font( item.getFont().getFontName(), 0, 16 ) );
        item.addActionListener( act );

        // Add a mouse listener for rollovers
        MouseActionBase rollover = new MouseActionBase() {
            public void enterMouse() { RNAstructure.setLabel( label ); }

            public void exitMouse() {
                RNAstructure.setLabel( "For help, press F1." );
            }
        };

        item.addMouseListener( rollover );
        return item;
    }

    /**
     * Set multiple sections visible.
     * For this method to work correctly, the sections must be listed in
     * numerical order.
     * @param list  the variable list of sections whose visibility to set
     *                 The list is one-indexed.
     */
    public void setSectionsVisible( Integer... list ) {
        int length = sections.size();

        // First, set all sections invisible, then set only the required
        // sections visible.
        for( int i = 1; i <= length; i++ ) {
            setSectionVisibility( i, false );
        }

        for( Integer part: list ) { setSectionVisibility( part, true ); }

        // If the last visible component of the menu is a separator, set the
        // separator invisible
        for( int i = getMenuComponentCount() - 1; i >= 0; i-- ) {
            if( getMenuComponent( i ).isVisible() ) {
                if( getMenuComponent( i ) instanceof JSeparator ) {
                    getMenuComponent( i ).setVisible( false );
                }
                break;
            }
        }
    }

    /**
     * Set a particular section of a menu visible or invisible. By default, all
     * menu sections start out visible.
     * @param section  the section whose visibility to set, one-indexed
     * @param view  true if the section is visible, false if not
     */
    public void setSectionVisibility( int section, boolean view ) {
        // Get the bounds of the section whose visiblity to change
        int start = ( section != 1 ) ? sections.get( section - 1 ) : 0;

        int end = ( section >= sections.size() ) ?
            getItemCount() : sections.get( section );

        // Change the visiblity of the entire section, including separators
        for( int i = start; i < end; i++ ) {
            getMenuComponent( i ).setVisible( view );
        }
    }
}