/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.util.ArrayList;

import javax.swing.BorderFactory;
import javax.swing.JCheckBox;
import javax.swing.JMenu;

/**
 * A class that handles traditional Dynalign calculations.
 * @author Jessica Reuter
 */
public class DynalignFoldWindow extends DynalignWindow {
    /**
     * A list of alignment constraints
     * Column 0: Index in sequence 1
     * Column 1: Index in sequence 2
     */
    private ArrayList<Integer[]> alignments;

    /**
     * The check box saying whether single base pair inserts are allowed
     */
    private JCheckBox box;

    /**
     * A check box saying whether a save file should be generated
     */
    private JCheckBox box2;

    /**
     * The input panel holding gap penalty data
     */
    private TextInputPanel inputPanel3;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -3770136130269724217L;

    /**
     * Constructor
     * @param acid  the type of nucleic acid this window is using
     */
    public DynalignFoldWindow( String acid ) {
        // Set titles
        super( acid + " Dynalign" );

        // Set the menu bar
        setSections( 1, 3, 4, 5, 7 );
        SmartMenuBar bar = RNAstructure.getBar();

	alignments = new ArrayList<Integer[]>();

        SmartMenu constraint1 =
            bar.getMenuMaker().createSequenceConstraintsMenu( 1 );
        constraint1.addMouseListener( new DynalignFoldWindowMenuGrabber() );

        SmartMenu constraint2 =
            bar.getMenuMaker().createSequenceConstraintsMenu( 2 );
        constraint2.addMouseListener( new DynalignFoldWindowMenuGrabber() );

        setMenus(
            bar.getMenuMaker().createTemperatureMenu(),
            constraint1, constraint2,
            bar.getMenuMaker().createDynalignAlignmentMenu()
        );
        RNAstructure.getBar().refresh( sections, menus );

        // Create the file input panel
        createFilePanel( 5, "Sequence File 1", "CT File 1", "Sequence File 2",
            "CT File 2", "Alignment File" );
        inputPanel.disableField( 1 );
        inputPanel.disableField( 3 );

        // Create the arrays of cutoffs for structure and alignment windows
        Integer[][] structureCutoffs = {
            { 1200, 20 },
            { 800, 15 },
            { 500, 11 },
            { 300, 7 },
            { 120, 5 },
            { 50, 3 },
            { 0, 2 }
        };

        Integer[][] alignmentCutoffs = {
            { 500, 3 },
            { 300, 2 },
            { 50, 1 },
            { 0, 0 }
        };

        // Create and set actions for buttons
        // First: Select the first sequence file and set its name.
        // Second: Set the name of the default output file for the first input
        //         sequence file.
        // Third: Select an output file other than the default, for the first
        //        input file, and set its name.
        // Fourth: Select the second sequence file and set its name.
        // Fifth: Set the name of the default output file for the second input
        //        sequence file.
        // Sixth: Select an output file other than the default, for the second
        //        input file, and set its name.
        // Seventh: Create a Dynalign object.
        // Eighth: Fill in the name of the default output alignment file, using
        //         the two previous output file names.
        // Ninth: Fill the folding parameters in.
        // Tenth: Select an alignment file other than default and set its name.
        // "First" and "Second" execute in sequence on button 0.
        // "Third" executes alone on button 1.
        // "Fourth", "Fifth", "Seventh", "Eighth", and "Ninth" execute in
        // sequence on button 2.
        // "Sixth" executes alone on button 3.
        // "Tenth" executes alone on button 4.
        OpenSaveSetter first =
            OpenSaveSetter.open( 0, "Sequence Files", "seq" );
        OutputFiller second = OutputFiller.fillSingle( 1, "ct" );
        OpenSaveSetter third =
            OpenSaveSetter.save( 1, "CT Files", "ct" );
        OpenSaveSetter fourth =
            OpenSaveSetter.open( 2, "Sequence Files", "seq" );
        OutputFiller fifth = OutputFiller.fillSingle( 3, "ct" );
        OpenSaveSetter sixth =
            OpenSaveSetter.save( 3, "CT Files", "ct" );
        DynalignObjectCreator seventh =
            new DynalignObjectCreator( 0, 2, 2, acid );
        OutputFiller eighth = OutputFiller.fillMultiple( 4, 1, 3, "ali" );
        FoldingParamFiller ninth = FoldingParamFiller.setVariable( "20", "20",
            structureCutoffs, alignmentCutoffs );
        OpenSaveSetter tenth =
            OpenSaveSetter.save( 4, "Alignment Files", "ali" );

        ChainedAction firstFileFill = new ChainedAction( first, second );
        ChainedAction dynalignSequence =
            new ChainedAction( fourth, fifth, seventh, eighth, ninth );

        inputPanel.setAction( 0, firstFileFill );
        inputPanel.setAction( 1, third );
        inputPanel.setAction( 2, dynalignSequence );
        inputPanel.setAction( 3, sixth );
        inputPanel.setAction( 4, tenth );

        // Create complex input with gap penalties and base pair inserts
        String[] gapLabel = { "Gap Penalty" };
        inputPanel3 = new TextInputPanel();
        inputPanel3.create( 1, 2, 5, gapLabel );
        inputPanel3.setBorder(
            BorderFactory.createEmptyBorder( 0, 10, 0, 0 ) );
        inputPanel3.getFields()[0].setText( "0.4" );

        box = new JCheckBox( "<html>Single BP Inserts<br>Allowed" );
        box.setSelected( true );
        box.setBorder( BorderFactory.createEmptyBorder( 0, 10, 0, 0 ) );

        box2 = new JCheckBox( "Generate Save File" );
        box2.setBorder( BorderFactory.createEmptyBorder( 0, 10, 0, 0 ) );

        // Create the structure parameters panel
        createStructurePanel();

        // Set components in their proper places
        setFill( 1 );
        setGrid( 2, 1 );
        placeComponent( 0, 0, inputPanel );

        setGrid( 1, 3 );
        placeComponent( 0, 1, inputPanel2 );

        setGrid( 1, 1 );
        setPad( 0, 20 );
        placeComponent( 1, 1, inputPanel3 );

        setInsets( 0, 10, 0, 0 );
        placeComponent( 1, 2, box );
        placeComponent( 1, 3, box2 );

        setGrid( 2, 1 );
        setFill( 2 );
        setPad( 25, 25 );
        setInsets( 0, 0, 10, 0 );
        placeComponent( 0, 4, createStartButton( new DynalignRunner() ) );
    }

    /**
     * Add a forced alignment.
     * @param seq1index  the index in sequence 1 to align
     * @param seq2Index  the index in sequence 2 to align
     */
    public void addForcedAlignment( int seq1Index, int seq2Index ) {
	Integer[] align = { seq1Index, seq2Index };
	alignments.add( align );
    }

    /**
     * Clear all forced alignments.
     */
    public void clearForcedAlignments() { alignments.clear(); }

    /**
     * Get the forced alignments.
     * Convert the list into a 2D array for easier handling.
     * @return  the forced alignments, as a 2D array
     */
    public Integer[][] getForcedAlignments() {
	int size = alignments.size();

	Integer[][] alignArray = new Integer[size][];
	for( int i = 0; i < size; i++ ) {
	    alignArray[i] = alignments.get( i );
	}

	return alignArray;
    }

    /**
     * Get an array which holds the box selections.
     * @return  the array of boolean selections.
     */
    public boolean[] getBoxSelections() {
        boolean[] selections = { box.isSelected(), box2.isSelected() };
        return selections;
    }

    /**
     * Get the gap penalty.
     * @return  the gap penalty
     */
    public float getGapPenalty() {
        return Float.parseFloat( inputPanel3.getFields()[0].getText() );
    }

    /**
     * Check if all required data fields have been filled in. This method is a
     * wrapper to call the TextInputPanel method checkValues; see
     * TextInputPanel.java for more information.
     * @return  true if all required data is filled, false if not
     */
    public boolean isReady() {
        return super.isReady() && inputPanel3.checkValues();
    }

    /**
     * An inner class that gets the last menu clicked on while this window
     * was active.
     * @author Jessica Reuter
     */
    class DynalignFoldWindowMenuGrabber extends MouseActionBase {
        /**
         * Grab the name of the most recent clicked menu
         */
        public void enterMouse() {
            JMenu menu = (JMenu)event.getSource();

            if( menu.isEnabled() ) {
                String lastMenuText = menu.getText();

                int num = Integer.parseInt( lastMenuText.substring(
                    lastMenuText.lastIndexOf( " " ) + 1 ) );

                currentStrand =
                    ( num == 1 ) ? dynalign.GetRNA1() : dynalign.GetRNA2();
            }
        }
    }

    /**
     * An inner class which allows a DynalignFoldWindow to be constructed from
     * the more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * The type of nucleic acid being used in Dynalign
         */
        private String acid;

        /**
         * Constructor
         * @param acid  the type of nucleic acid
         */
        public Factory( String acid ) { this.acid = acid; }

        /**
         * Create a DynalignFoldWindow.
         * @return  a new DynalignFoldWindow
         */
        public DynalignFoldWindow createWindow() {
            return new DynalignFoldWindow( acid );
        }
    }
}
