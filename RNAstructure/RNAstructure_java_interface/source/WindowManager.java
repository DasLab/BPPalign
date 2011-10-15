/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.Window;

/**
 * A class that removes a window from the memory list, refreshing the main
 * screen with the details of the window before it. This class also can
 * refresh details of an existing window freshly selected.
 * @author Jessica Reuter
 */
public class WindowManager extends WindowActionBase {
    /**
     * Set a particular window's menus and focus on activation.
     */
    public void activateWindow() {
        enableWindow( (MenuChangingWindow)event.getSource() );
    }

    /**
     * Remove a window from the list of main windows and focus the next.
     */
    public void closeWindowDoAfter() {
        RNAstructure.getWindowList()
            .remove( (MenuChangingWindow)event.getSource() );
        Window window = ( RNAstructure.getWindowList().size() != 0 ) ?
	    RNAstructure.getWindowList().getLast() : null;

        // If other windows are open, focus the most recent
        // If not, focus the main frame
        if( window != null ) {
            enableWindow( (MenuChangingWindow)window );
        } else {
            RNAstructure.getBar().refresh();
            RNAstructure.getFrame().setTitle( "RNAstructure" );
            RNAstructure.getToolBar().setButtonsEnabled( false, false );
            RNAstructure.getFrame().requestFocus();
        }
    }

    /**
     * Enable a window and its menus.
     * @param window  the MenuChangingWindow to enable
     */
    private void enableWindow( MenuChangingWindow top ) {
	if( top.sections != null ) {
	    RNAstructure.getBar().refresh( top.sections, top.menus );
        } else {
	    RNAstructure.getBar().refresh();
	}

        // Enable menus if the current strand is not null. There are three
	// exceptions to this rule where menus must be active:
        //    1. NewSequenceWindow
        //    2. Imager
	//    3. OligoWalkResultsWindow
        boolean autoEnable =
            top instanceof NewSequenceWindow ||
            top instanceof Imager ||
	    top instanceof OligoWalkResultsWindow;

	boolean enable = autoEnable;
	if( enable == false ) {
	    if( top instanceof DynalignWindow ) {
		enable = ((DynalignWindow)top).getWorkingStrand() != null;
	    } else if( top instanceof OligoWindow ) {
		enable = ((OligoWindow)top).getOligoObject() != null;
	    } else {
		enable =
		    ((InputWindow)top).getDataHolder().getStrand() != null;
	    }
	}

        if( enable ) { RNAstructure.getBar().enableMenus(); }

        top.setTitles( top.getTitle() );
        boolean[] enabledArray = top.getEnabledPermissions();
        RNAstructure.getToolBar()
            .setButtonsEnabled( enabledArray[0], enabledArray[1] );

	RNAstructure.setCurrentWindow( top );
    }
}
