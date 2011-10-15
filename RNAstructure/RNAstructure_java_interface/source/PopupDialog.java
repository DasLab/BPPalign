/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.BorderLayout;
import java.awt.Dialog;

import javax.swing.BorderFactory;
import javax.swing.Icon;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.UIManager;

/**
 * A class that makes confirmation and input dialogs for the RNAstructure GUI.
 * This class and its methods should not be called directly.
 * @author Jessica Reuter
 */
public class PopupDialog extends BaseWindow {
    /**
     * The text field for the input dialog
     */
    private static JTextField field;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -3112803059649834936L;

    /**
     * Private constructor, used to set defaults
     * @param error  true if this dialog will show an error, false if not
     */
    private PopupDialog( Boolean error ) {
        // Create window defaults
        super( RNAstructure.getCurrentWindow() );
        setLayout( new BorderLayout() );
        setModalityType( Dialog.ModalityType.APPLICATION_MODAL );
        setResizable( false ); 

        // Set the title based on the error state
        if( error == null ) { setTitle( "RNAstructure Input" ); }
        else if( error ) { setTitle( "RNAstructure Error" ); }
        else { setTitle( "RNAstructure" ); }
    }

    /**
     * Create a message panel on the dialog.
     * @param iconRequired  true if an icon is required, false if not
     * @param message  the message on the panel
     * @return  the message panel
     */
    private JPanel createMessagePanel( boolean iconRequired, String message ) {
        JPanel mainPanel = new JPanel( new BorderLayout() );

        // Create the message label
        JLabel label = ( iconRequired ) ?
            new JLabel( "<html><br/>".concat( message ) ) :
            new JLabel( "<html><center><br/>".concat( message ) );
        label.setHorizontalAlignment( ( iconRequired ) ?
            JLabel.LEFT :
            JLabel.CENTER );
        label.setBorder( ( iconRequired ) ?
            BorderFactory.createEmptyBorder( 0, 0, 0, 10 ) :
            BorderFactory.createEmptyBorder( 0, 10, 0, 10 ) );

        // If an icon is required, a more complex panel must be created
        if( iconRequired ) {
            // Create the appropriate icon
            String title = getTitle();

            Icon icon = ( title.equals( "RNAstructure" ) ) ?
                UIManager.getIcon( "OptionPane.warningIcon" ) :
                ( title.equals( "RNAstructure Error" ) ) ?
                    UIManager.getIcon( "OptionPane.errorIcon" ) :
                    UIManager.getIcon( "OptionPane.informationIcon" );

            // Create the icon label
            JLabel image = new JLabel();
            image.setIcon( icon );
            image.setBorder(
                BorderFactory.createEmptyBorder( 10, 10, 0, 10 ) );

            // If this is an input window, add a text field for input
            JPanel textPanel = null;
            if( title.equals( "RNAstructure Input" ) ) {
                textPanel = new JPanel( new BorderLayout() );
                textPanel.add( field, BorderLayout.CENTER );
                textPanel.add( label, BorderLayout.NORTH );
                textPanel.setBorder( label.getBorder() );
            }

            // Create the encompassing panel and add it to the main panel
            JPanel iconPanel = new JPanel( new BorderLayout() );
            iconPanel.add( image, BorderLayout.WEST );

            if( textPanel == null ) { iconPanel.add( label ); }
            else { iconPanel.add( textPanel ); }

            mainPanel.add( iconPanel );
        }

        // If no icon required, simply add the message label to the main panel
        else { mainPanel.add( label ); }

        return mainPanel;
    }

    /**
     * Finalize dialog parameters
     */
    private void finish() {
        pack();
        setLocationRelativeTo( null );
        setVisible( true );
    }

    /**
     * Create a dialog which prompts the user for input
     * @param message  the message on the dialog
     * @param input  the field's initial input
     */
    public static String input( String message, String input ) {
        // Initialize the input field
        field = new JTextField();
        if( input != null ) { field.setText( input ); }

        // Create the action to cancel input
        ChainedAction cancelAction = new ChainedAction(
            new SetInputField( field ), new CloseAction() );

        // Create the input button panel (OK, Cancel)
        ButtonPanel panel = new ButtonPanel( 2 );
        panel.setAction( new CloseAction(), 0, "OK" );
        panel.setAction( cancelAction, 1, "Cancel" );
        panel.addButtons();

        // Create the dialog and return input entered into it
        PopupDialog.msg( message, panel, null, true );
        return field.getText();
    }

    /**
     * Create a dialog, with the proper message, ButtonPanel, error state, and
     * icon state. The error Boolean must be the boolean object wrapper rather
     * than the primitive type because a value of null for the error signifies
     * an input dialog (which also holds a message).
     * @param message  the message to display on the dialog
     * @param panel  the button panel needed for this dialog
     * @param error  true if this is an error dialog, false if not
     * @param iconRequired  true if an icon is required, false if not
     */
    public static void msg( String message, ButtonPanel panel, Boolean error,
                            boolean iconRequired ) {
        PopupDialog dialog = new PopupDialog( error );
        dialog.add( dialog.createMessagePanel( iconRequired, message ) );
        dialog.add( panel, BorderLayout.SOUTH );
        dialog.finish();
    }

    /**
     * Create a button panel that holds a single "OK" button
     * @return the panel
     */
    public static ButtonPanel OK() {
        ButtonPanel panel = new ButtonPanel( 1 );
        panel.setAction( new CloseAction(), 0, "OK" );
        panel.addButtons();
        return panel;
    }

    /**
     * Create a button panel with both an "OK" button and a "Cancel" button
     * @param ok  the action for the "OK" button
     * @return the panel
     */
    public static ButtonPanel OKCancel( ActionBase ok ) {
        ButtonPanel panel = new ButtonPanel( 2 );
        panel.setAction( ok, 0, "OK" );
        panel.setAction( new CloseAction(), 1, "Cancel" );
        panel.addButtons();
        return panel;
    }

    /**
     * Create a button panel with three buttons: "Yes", "No", and "Cancel".
     * @param yesAction  the action for when the user presses "Yes"
     * @param noAction  the action for when the user presses "No"
     * @return  the panel
     */
    public static ButtonPanel yesNoCancel( ActionBase yesAction,
                                           ActionBase noAction ) {
        ButtonPanel panel = new ButtonPanel( 3 );
        ActionBase[] acts = { yesAction, noAction, new CloseAction() };
        String[] labels = { "Yes", "No", "Cancel" };
        panel.setActions( acts, labels );
        panel.addButtons();
        return panel;
    }
}
