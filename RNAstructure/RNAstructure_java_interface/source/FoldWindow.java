/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JTextField;

/**
 * A window enabling the user to input data allowing folding of a single strand
 * or multiple strands of nucleic acids
 * @author Jessica Reuter
 */
public class FoldWindow extends InputWindow {
    /**
     * A check box that holds the option of creating a save file.
     * This box is used in all subclasses of FoldWindow except RefoldWindow.
     */
    protected JCheckBox box;

    /**
     * The input panel holding file name input
     */
    protected TextInputPanel inputPanel;

    /**
     * The input panel holding folding parameters
     */
    protected TextInputPanel inputPanel2;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -4920389665645273479L;

    /**
     * Constructor
     * @param title  the title of the window
     */
    protected FoldWindow( String title ) { setTitles( title ); }

    /**
     * Create the box that gives the option to generate a save file.
     */
    protected void createBox() { box = new JCheckBox( "Generate Save File" ); }

    /**
     * Create the main input panel.
     * @param amount  the number of rows in the main input panel
     * @param names  the array of button names
     */
    protected void createMainPanel( int amount, String... names ) {
        inputPanel = new TextInputPanel();
        inputPanel.create( amount, 1, 30, names );
        setBorder( inputPanel, 1, null );
    }

    /**
     * Create the folding parameters input panel which is not involved in
     * suboptimal structure folding.
     */
    protected void createFoldingPanel() {
        String[] names = {
            "Max % Energy Difference",
            "Max Number of Structures",
            "Window Size"
        };

        createParamPanel( 3, names );
    }

    /**
     * Create one of two suboptimal structure parameter panels.
     * @param amount  the number of rows in the parameter panel
     * @param names  the array of button names
     */
    private void createParamPanel( int amount, String[] names ) {
        inputPanel2 = new TextInputPanel();
        inputPanel2.create( amount, 2, 15, names );
        setBorder( inputPanel2, 2, "Suboptimal Structure Parameters" );
    }

    /**
     * Create folding parameters input panel which is involved in suboptimal
     * structure folding
     */
    protected void createSuboptimalPanel() {
        String[] names = { 
            "Max % Energy Difference",
            "Max Absolute Energy Difference"
        };

        createParamPanel( 2, names );
    }

    /**
     * Get the parameter panel.
     * @return  the parameter panel
     */
    public TextInputPanel getParameterPanel() { return inputPanel2; }

    /**
     * Check the input panels to see if all required files and parameter values
     * are present. This method is a wrapper to call the TextInputPanel method
     * checkValues on both input panels; see TextInputPanel.java for more
     * information.
     * @return  true if all values are defined, false if not
     */
    protected boolean isReady() {
        return ( inputPanel.checkValues() && inputPanel2.checkValues() );
    }

    /**
     * Place the components of this FoldWindow into their proper positions.
     * @param action  the action on the start button
     */
    protected void setComponents( StrandStarter action ) {
        setGrid( 2, 1 );
        setFill( 2 );
        placeComponent( 0, 0, inputPanel );

        if( box != null ) { placeComponent( 0, 1, box ); }

        setGrid( 1, 1 );
        placeComponent( 0, 2, inputPanel2 );

        JButton startButton = createStartButton( action );

        setPad( 25, 25 );
        setInsets( 0, 0, 0, 10 );
        placeComponent( 1, 2, startButton );
    }

    /**
     * Set all data from this FoldWindow into its data holder.
     * @throws Exception  if an error occurred during data setting
     */
    protected void setData() throws Exception {
        // Get the arrays of text fields and their lengths
        JTextField[] one = inputPanel.getFields();
        JTextField[] two = inputPanel2.getFields();

        int len1 = one.length;
        int len2 = two.length;

        Integer size = (Integer)holder.getData().get( 0 );
        removeData( 0 );

        addData( one[len1 - 1].getText() );

        for( int i = 0; i < len2; i++ ) {
           addData( Double.parseDouble( two[i].getText() ) );
        }
        addData( size );

        if( box != null ) { addData( box.isSelected() ); }
    }

    /**
     * Set the menu bar and its menus.
     * @param forceMenu  the unique force menu for a window type
     */
    protected void setMenuBar( SmartMenu forceMenu ) {
        setSections( 1, 3, 4, 5, 7 );
        setMenus( RNAstructure.getBar().getMenuMaker().createTemperatureMenu(),
                  forceMenu,
                  RNAstructure.getBar().getMenuMaker().createMaxLoopMenu() );
        RNAstructure.getBar().refresh( sections, menus );
    }
}
