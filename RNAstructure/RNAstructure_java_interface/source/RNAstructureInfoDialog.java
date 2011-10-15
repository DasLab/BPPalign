/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class which creates dialogs for the RNAstructure GUI. This class uses
 * static methods of the class PopupDialog, which does the major work of
 * creating the dialogs.
 * @author Jessica Reuter
 */
public class RNAstructureInfoDialog {
    /**
     * Create a message dialog that shows a confirmation message, with options
     * Yes, No, or Cancel.
     * @param yesAction  the action for when the user presses "Yes"
     * @param noAction  the action for when the user presses "No"
     */
    public static void choose( String message, ActionBase yesAction,
                               ActionBase noAction ) {
        ButtonPanel panel = PopupDialog.yesNoCancel( yesAction, noAction );
        PopupDialog.msg( message, panel, false, false );
    }

    /**
     * Create a message dialog.that shows a confirmation message, with options
     * OK or Cancel
     * @param message  the message to be shown
     * @param action  the action listener for the OK button
     */
    public static void confirm( String message, ActionBase action ) {
        ButtonPanel buttonPanel = PopupDialog.OKCancel( action );
        PopupDialog.msg( message, buttonPanel, false, false );
    }

    /**
     * Create a message dialog.that shows a simple error message.
     * @param message  the message to be shown
     */
    public static void error( String message ) {
        PopupDialog.msg( message, PopupDialog.OK(), true, true );
    }

    /**
     * Create a message dialog.that allows the user to input data
     * @param message  the message to be shown
     * @param inputString  the default string for the input dialog
     * @return  the data the user inputs, as a string
     */
    public static String input( String message, String inputString ) {
        return PopupDialog.input( message, inputString );
    }

    /**
     * Create a message dialog.that shows a simple message.
     * @param message  the message to be shown
     */
    public static void message( String message ) {
        PopupDialog.msg( message, PopupDialog.OK(), false, false );
    }

    /**
     * Create a message dialog.that shows a simple warning message.
     * @param message  the message to be shown
     */
    public static void warning( String message ) {
        PopupDialog.msg( message, PopupDialog.OK(), false, true );
    }
}
