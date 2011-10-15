/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridLayout;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextPane;

/**
 * A class responsible for creating and displaying a window that enables a user
 * to input data for creating a new sequence of either RNA or DNA, or to view
 * an existing sequence of nucleic acids.
 * Since this class can hold existing sequences as well as create new ones, its
 * name is somewhat of a misnomer.
 * @author Jessica Reuter
 */
public class NewSequenceWindow extends MenuChangingWindow {
    /**
     * Whether any part of the information in this window has been edited
     */
    private boolean edited;

    /**
     * The input areas in this window.
     */
    private JTextPane[] inputAreas;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -1467391257056258674L;

    /**
     * The button which turns out loud reading of sequence on and off
     */
    private JButton talkButton;

    /**
     * Constructor
     * In general, the title of this window is "New Sequence", but a sequence
     * that already exists can also be opened here, in which case its path name
     * will be the window's title.
     * @param title  the title on the window
     */
    public NewSequenceWindow( String title ) {
        // Set titles and defaults
        setTitles( title );
        edited = false;

        // Set menu bar and toolbar
        setSections( 1, 2, 3, 4, 5, 7 );
        SmartMenuBar bar = RNAstructure.getBar();
        setMenus( bar.getMenuMaker().createEditMenu(),
                  bar.getMenuMaker().createReadMenu() );
        RNAstructure.getBar().refresh( sections, menus );
        RNAstructure.getBar().enableMenus();
        RNAstructure.getToolBar().setButtonsEnabled( true, false );

        // Make the three major input areas and the button panel
        JPanel[] inputPanels = makeInputAreas();
        JPanel buttonPanel = makeOptionPanel();

        // Add an extra window listener so the user doesn't accidentally lose a
        // file that has been edited
        addWindowListener( new EditChecker() );

        // Add components in their proper places
        setFill( 1 );
        setInsets( 10, 10, 10, 10 );
        for( int i = 0; i < 3; i++ ) {
            placeComponent( 0, i, inputPanels[i] );
         }

        placeComponent( 0, 3, buttonPanel );

        setDefaultCloseOperation( JDialog.DO_NOTHING_ON_CLOSE );
    }

    /**
     * Get a particular input area.
     * @param index  the index of the area, from top to bottom (one-indexed)
     * @return  the selected area
     */
    public JTextPane getArea( int index ) { return inputAreas[index-1]; }

    /**
     * Get the button which switches the talk pane on and off.
     * @return  the talker button
     */
    public JButton getTalkButton() { return talkButton; }

    /**
     * Check if the NewSequenceWindow has been edited.
     * @return  true if the window has been edited, false if not
     */
    public boolean isEdited() { return edited; }

   /**
    * Create the three main input areas for the window.
    * @return  the array of input panels
    */
   private JPanel[] makeInputAreas() {
       String[] labels = { "Title:", "Comment:", 
                           "Sequence: (Nucleotides in lower case are forced " +
                               "single stranded in structure predictions.)" };
       JPanel[] inputPanels = new JPanel[3];
       inputAreas = new JTextPane[3];

       for( int i = 0; i <= 2; i++ ) {
           // Create the input label and text area
           JPanel panel = new JPanel( new BorderLayout() );
           panel.add( new JLabel( labels[i] ), BorderLayout.NORTH );

           JTextPane area = ( i != 2 ) ? new JTextPane() : new TalkPane();
           area.setFont( new Font( "Monospaced", 0, 14 ) );
           area.setBackground( Color.WHITE );
           area.addKeyListener( new EditListener() );

           // Add the panel to a scroll pane
           JScrollPane pane = new JScrollPane( area );
           pane.setHorizontalScrollBarPolicy( ( i == 0 ) ?
               JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS :
               JScrollPane.HORIZONTAL_SCROLLBAR_NEVER );
           pane.setVerticalScrollBarPolicy( ( i == 0 ) ?
               JScrollPane.VERTICAL_SCROLLBAR_NEVER :
               JScrollPane.VERTICAL_SCROLLBAR_ALWAYS );

           // Set the size of the scroll pane
           int height = ( i == 0 ) ? 40 : ( i == 1 ) ? 50 : 200;
           pane.setPreferredSize( new Dimension( 100, height ) );

           // Set the input panel and text pane in their data structures
           panel.add( pane );
           inputPanels[i] = panel;
           inputAreas[i] = area;
       }

       return inputPanels;
   }

   /**
    * Create the option panel at the bottom of the window.
    * @return  the button panel
    */
   private JPanel makeOptionPanel() {
       JPanel buttonPanel = new JPanel( new GridLayout( 1, 0 ) );

       // Create button that formats a sequence
       JButton formatButton = new JButton( "Format Sequence" );
       formatButton.addActionListener( new Formatter() );
       buttonPanel.add( formatButton );

       // Create button that turns out loud reading of sequence on and off
       talkButton = new JButton( "Read Sequence" );
       talkButton.addActionListener( ReadAloud.readAutomatic() );
       buttonPanel.add( talkButton );

       // Create buttons that fold DNA and RNA
       String[] acids = { "DNA", "RNA" };
       for( int i = 0; i < 2; i++ ) {
           JButton acidButton = new JButton( "Fold as " + acids[i] );
           acidButton.addActionListener( new PrepareSequenceFold( acids[i] ) );
           buttonPanel.add( acidButton );
       }

       return buttonPanel;
   }

   /**
    * Set the edited property of the window.
    * @param edited  true if the window was edited, false if not
    */
   public void setEdited( boolean edited ) { this.edited = edited; }

   /**
    * Switch the text of the button which turns talking on and off to reflect
    * its curent state.
    */
   public void switchTalkButton() {
       boolean idle = talkButton.getText().equals( "Read Sequence" );
       if( idle ) { talkButton.setText( "Stop Reading" ); }
       else { talkButton.setText( "Read Sequence" ); }
   }

    /**
     * An inner class which allows a NewSequenceWindow to be constructed from
     * the more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * The title of the NewSequenceWindow
         */
        private String title;

        /**
         * Constructor
         * @param title  the title of the new window
         */
        public Factory( String title ) { this.title = title; }

        /**
         * Create a NewSequenceWindow.
         * @return  a new NewSequenceWindow
         */
        public NewSequenceWindow createWindow() {
            return new NewSequenceWindow( title );
        }
    }
}