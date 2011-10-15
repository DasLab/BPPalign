/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.GridLayout;
import java.awt.Insets;

import java.util.ArrayList;

import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;
import javax.swing.border.EmptyBorder;

public class PostscriptDialogWindow extends DialogWindow {
    /**
     * The array of check boxes in this window that selects
     * structure drawing types.
     */
    private final JCheckBox[] boxesA;

    /**
     * The array of check boxes in this window that selects annotation
     * types.
     */
    private final JCheckBox[] boxesB;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = 5207689097905284764L;

    /**
     * Constructor
     */
    public PostscriptDialogWindow() {
        // Set the title
        // Create the main data label
        setTitle( "Postscript Output" );
        JLabel label = new JLabel( "Select type(s) of Postscript below." );
	label.setBorder( new EmptyBorder( 10, 0, 10, 0 ) );

        // Create the main check box panel.
        JPanel panel = new JPanel( new GridLayout( 2, 1 ) );

	// Create the structure drawing type panel.
	JPanel panelA = new JPanel( new GridLayout( 1, 2 ) );

	String[] labelsA = {
	    "<html><center>Radial Structures<br/>" +
	        "(No Pseudoknots Only)",
	    "<html><center>Circular Structures"
	};

	boxesA = new JCheckBox[2];
	for( int i = 1; i <= 2; i++ ) {
	    JCheckBox box = new JCheckBox( labelsA[i-1] );
	    box.setHorizontalAlignment( SwingConstants.CENTER );
	    box.setVerticalAlignment( SwingConstants.TOP );
	    box.setHorizontalTextPosition( SwingConstants.CENTER );
	    box.setVerticalTextPosition( SwingConstants.BOTTOM );
	    box.setMargin( new Insets( 0, 0, 5, 0 ) );
	    boxesA[i-1] = box;
	    panelA.add( boxesA[i-1] );
	}
	panel.add( panelA );

	// Create the annotation type panel.
	JPanel panelB = new JPanel( new GridLayout( 1, 3 ) );

        String[] labelsB = {
            "<html><center>Unannotated",
            "<html><center>Probability<br/>Annotated",
            "<html><center>SHAPE<br/>Annotated",
        };

        boxesB = new JCheckBox[3];
        for( int i = 1; i <= 3; i++ ) {
            JCheckBox box = new JCheckBox( labelsB[i-1] );
	    box.setHorizontalAlignment( SwingConstants.CENTER );
	    box.setVerticalAlignment( SwingConstants.TOP );
	    box.setHorizontalTextPosition( SwingConstants.CENTER );
	    box.setVerticalTextPosition( SwingConstants.BOTTOM );
	    box.setMargin( new Insets( 5, 0, 0, 0 ) );
	    boxesB[i-1] = box;
	    panelB.add( boxesB[i-1] );
        }
	panel.add( panelB );

	// Set the check boxes checked that apply to the current structure
	// being shown.
	DrawStructure structure =
	    (DrawStructure)(((Imager)RNAstructure.getCurrentWindow())
	    .getSketcher());

	int typeIndex1 = 0;
	if( structure instanceof DrawStructureCircular ) { typeIndex1 = 1; }
	boxesA[typeIndex1].setSelected( true );

	int typeIndex2 = 0;
	boolean isProbability =
	    !structure.getProbabilityAnnotationFile().equals( "" );
	boolean isSHAPE = !structure.getSHAPEAnnotationFile().equals( "" );
	if( isProbability ) { typeIndex2 = 1; }
	else if( isSHAPE ) { typeIndex2 = 2; }
	boxesB[typeIndex2].setSelected( true );

        // Create the button panel
        ActionBase select = new ActionBase() {
            public void act() {
		ArrayList<String> types = new ArrayList<String>();
		if( boxesA[0].isSelected() ) { types.add( "Radial" ); }
		if( boxesA[1].isSelected() ) { types.add( "Circular" ); }
		if( boxesB[0].isSelected() ) { types.add( "Plain" ); }
		if( boxesB[1].isSelected() ) { types.add( "Probability" ); }
		if( boxesB[2].isSelected() ) { types.add( "SHAPE" ); }

                PostscriptStructureFileWriter.exportToPostscript( types )
		    .act();
                ((PostscriptDialogWindow)boxesA[0].getTopLevelAncestor())
                    .dispose();
            }
        };

        ButtonPanel buttonPanel = new ButtonPanel( 2 );
        buttonPanel.setAction( select, 0, "OK" );
        buttonPanel.setAction( new CloseAction(), 1, "Cancel" );
        buttonPanel.addButtons();

        // Place components in the right places
        placeComponent( 0, 0, label );
        placeComponent( 0, 1, panel );
        placeComponent( 0, 2, buttonPanel );
    }

    /**
     * An inner class which allows a PostscriptDialogWindow to be constructed
     * from the more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * Create a PostscriptDialogWindow for structure type selection.
         * @return  a new SpinnerWindow
         */
        public PostscriptDialogWindow createWindow() {
            return new PostscriptDialogWindow();
        }
    }
}