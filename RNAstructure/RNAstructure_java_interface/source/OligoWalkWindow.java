/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.GridLayout;

import java.util.Enumeration;

import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/**
 * A class that handles creation of an input window for OligoWalk.
 * @author Jessica Reuter
 */
public class OligoWalkWindow extends OligoWindow {
    /**
     * Check box that tells whether to include target suboptimal structures
     */
    private JCheckBox box;

    /**
     * The combobox that allows a user to select the oligo concentration
     */
    private JComboBox comboBox;

    /**
     * The text field that holds the oligo concentration
     */
    private JTextField field;

    /**
     * The button group that tells which mode to use
     */
    private ButtonGroup group;

    /**
     * Boolean, true if this window handles DNA, false if not
     */
    private boolean isDNA;

    /**
     * The oligo length input spinner
     */
    private JSpinner spinner;

    /**
     * The array of spinners used in target structure limits
     */
    private JSpinner[] spinners;

    /**
     * The input panel that holds siRNA data (may or may not be present)
     */
    private TextInputPanel siRNAPanel;

    /**
     * A serialized ID (required when extending GUI classes)
     */
    private static final long serialVersionUID = 8231065364850273876L;

    /**
     * Constructor
     * @param type  the type of nucleic acid being analyzed
     */
    public OligoWalkWindow( String type ) {
        setTitles( type + " OligoWalk" );
	isDNA = false;
        
        // Set sections and menus
        setSections( 1, 3, 4, 5, 7 );
        SmartMenuBar bar = RNAstructure.getBar();
        setMenus( bar.getMenuMaker().createTemperatureMenu() );
        bar.refresh( sections, menus );

        // Create the file name input panel
        createFilePanel();

        // Create the panel of mode buttons
        JPanel modePanel = createModePanel();

        // Create the option to include suboptimal structures
        box = new JCheckBox( "Include Target Suboptimal Structures in Free " +
            "Energy Calculation" );
        
	// Create the oligomer chemistry panel
	JPanel oligomerPanel = createOligomerPanel();

        // Create the panel which holds oligo length and concentration
        JPanel oligoPanel = createOligoPanel();

        // Create target structure limits panel
        JPanel limitPanel = createLimitsPanel();

        // Create the siRNA panel, if applicable
        if( type.equals( "siRNA" ) ) {
            siRNAPanel = new TextInputPanel();
            siRNAPanel.create( 3, 2, 10,
              "Antisense Strand Unimolecular Folding",
              "Target Opening Free Energy",
              "First Nearest Neighbor Free Energy");
            setBorder( siRNAPanel, 2,
                "Screening Parameters for siRNA Selection" );
        }

        // Place all components in their proper places
        int gridy = 0;
        setPad( 0, 20 );
        setGrid( 2, 1 );

        placeComponent( 0, gridy++, panel );
        placeComponent( 0, gridy++, modePanel );

        if( type.equals( "siRNA" ) ) {
            placeComponent( 0, gridy++, siRNAPanel );
        }

	placeComponent( 0, gridy++, oligomerPanel );
        placeComponent( 0, gridy++, box );
        placeComponent( 0, gridy++, oligoPanel );

        setGrid( 1, 1 );
        setPad( 80, 0 );
        setInsets( 0, 10, 0, 0 );
        placeComponent( 0, gridy, limitPanel );

        setPad( 25, 25 );
        setInsets( 5, 0, 0, 10 );
        placeComponent( 1, gridy, createStartButton( new RunOligoWalk() ) );
    }
    
    /**
     * Create the top input panel that allows for input of file names.
     */
    private void createFilePanel() {
        // Create the panel
        panel = new TextInputPanel();
        panel.create( 2, 1, 30, "Input File", "Report File" );
        setBorder( panel, 1, null );

	// Add actions on buttons
	// First: Determine the input file and set its name.
	// Second: Create the OligoWalkObject.
	// Third: Set the default output file.
	// Fourth: Edit the GUI models and put in default information.
	// Fifth: Determine an output file other than the default and set it.
	// "First", "Second". "Third", "Fourth" run in sequence on button 0.
	// "Fifth" executes alone on button 1.
	OpenSaveSetter first = OpenSaveSetter.open( 0, "CT Files", "ct" );
	OligoWalkObjectCreator second = new OligoWalkObjectCreator( !isDNA );
	OutputFiller third = OutputFiller.fillSingle( 1, "rep" );
	ActionBase fourth = new ActionBase() {
	    public void act() {
		if( !panel.getFields()[0].getText().trim().equals( "" ) ) {
		    ((OligoWalkWindow)RNAstructure.getCurrentWindow())
			.editModels();
		}
	    }
	};
	OpenSaveSetter fifth =
	    OpenSaveSetter.save( 1, "Report Files", "rep" );

	ChainedAction sequence =
	    new ChainedAction( first, second, third, fourth );

	panel.setAction( 0, sequence );
	panel.setAction( 1, fifth );
    }
    
    /**
     * Create the panel that holds constraints for target structure limits.
     * @return  the limits panel
     */
    private JPanel createLimitsPanel() {
        // Create the panel
        JPanel limitPanel = new JPanel( new GridLayout( 1, 2 ) );
        setBorder( limitPanel, 2, "Target Structure Limits for Walk" );
        
        // Create each of the two spinners and add them to the panel.
        spinners = new JSpinner[2];
        for( int i = 0; i <= 1; i++ ) {
            JPanel panelC = new JPanel();

            if( i == 0 ) { panelC.add( new JLabel( "Start: " ) ); }
            else { panelC.add( new JLabel( "Stop: " ) ); }

            JSpinner targetSpinner =
                new JSpinner( new SpinnerNumberModel( 0, 0, 0, 0 ) );
            ((JSpinner.DefaultEditor)targetSpinner.getEditor())
                .getTextField().setColumns( 5 );
            
            spinners[i] = targetSpinner;
            panelC.add( targetSpinner );
            limitPanel.add( panelC );
        }
        
        return limitPanel;
    }
    
    /**
     * Create the input panel which sets the calculation mode.
     * @return  the input panel
     */
    private JPanel createModePanel() {
        // Create the panel and the button group on it
        JPanel modePanel = new JPanel( new GridLayout( 3, 1 ) );
        setBorder( modePanel, 2, "Mode" );
        group = new ButtonGroup();
        
        // For each mode button, place it both in the group and on the panel
        String[] modeData = {
            "Break Local Structure",
            "Refold Whole RNA for Each Oligomer (Slowest)",
            "Do Not Consider Target Structure (Fastest)" };

        for( int i = 0; i <= 2; i++ ) {
            JRadioButton button = new JRadioButton( modeData[i] );
	    if( i == 0 ) { button.setSelected( true ); }
            group.add( button );
            modePanel.add( button );
        }
        
        return modePanel;
    }
    
    /**
     * Create the panel that holds oligo length and concentration information.
     * @return  the oligo panel
     */
    private JPanel createOligoPanel() {
        // Create the oligo length panel
        JPanel panelA = new JPanel();
        JLabel label1 = new JLabel( "Oligo Length: " );
        label1.setHorizontalAlignment( JLabel.RIGHT );
        panelA.add( label1 );

        // Create the oligo length spinner. with a change listener that allows
        // target structure limits to change with the length
        spinner = new JSpinner( new SpinnerNumberModel( 18, 1, 10000, 1 ) );
        spinner.addChangeListener( new ChangeListener() {
            public void stateChanged( ChangeEvent e ) { editModels(); }
        });
        panelA.add( spinner );

        // Create the oligo concentration panel
        JPanel panelB = new JPanel();
        JLabel label2 = new JLabel( "Oligo Concentration: " );
        label2.setHorizontalAlignment( JLabel.RIGHT );
        panelB.add( label2 );

        field = new JTextField( "1" );
        field.setColumns( 4 );
        panelB.add( field );

        String[] amounts = { "mM", "uM", "nM", "pM" };
        comboBox = new JComboBox( amounts );
        comboBox.setSelectedItem( amounts[1] );
        panelB.add( comboBox );

        // Combine oligo length panel and oligo concentration panel
        JPanel oligoPanel = new JPanel( new GridLayout( 2, 1 ) );
        oligoPanel.add( panelA );
        oligoPanel.add( panelB );
        return oligoPanel;
    }

    /**
     * Create the panel which dictates oligomer chemistry.
     * @return  the oligomer panel
     */
    private JPanel createOligomerPanel() {
	JPanel panel = new JPanel( new GridLayout( 1, 2 ) );
	String[] chemistry = { "DNA", "RNA" };
	for( int i = 1; i <= 2; i++ ) {
	    JRadioButton button = new JRadioButton( chemistry[i-1] );
	    if( i == 2 ) { button.setSelected( true ); }

	    button.addActionListener( new ActionBase() {
		public void act() {
		    String acid = event.getActionCommand();
		    isDNA = acid.equals( "DNA" );
		}
	    });

	    panel.add( button );
	}

        setBorder( panel, 2, "Oligomer Chemistry" );
	return panel;
    }
    
    /**
     * Change the spinner number models to coincide with changes in either the
     * nucleic acid strand or other spinners.
     */
    protected void editModels() {
        int length = RNAstructure.getCurrentStrand().GetSequenceLength();
        
        int spinLength = (Integer)spinner.getModel().getValue();
        int max = length - spinLength + 1;
        
        spinner.setModel(
	    new SpinnerNumberModel( spinLength, 1, length, 1 ) );
        ((JSpinner.DefaultEditor)spinner.getEditor())
            .getTextField().setColumns( 5 );
        
        spinners[0].setModel( new SpinnerNumberModel( 1, 1, max, 1 ) );
        ((JSpinner.DefaultEditor)spinners[0].getEditor())
            .getTextField().setColumns( 5 );
        
        spinners[1].setModel( new SpinnerNumberModel( max, 1, max, 1 ) );
        ((JSpinner.DefaultEditor)spinners[1].getEditor())
            .getTextField().setColumns( 5 );
    }
    

    /**
     * Check to see if all required data has been added into this window.
     * @return  true if all required data is present, false if not
     */
    public boolean isReady() {
	boolean text = panel.checkValues();
	boolean conc = !field.getText().equals( "" );
	return text && conc;
    }

    /**
     * Set the required data for this window.
     * @throws Exception  if an error occurs during data setting
     */
    public void setData() throws Exception {
	removeData( 0 );

	addData( (Integer)spinner.getModel().getValue() );
	addData( isDNA );

	int index = 0;
	for( Enumeration e = group.getElements(); e.hasMoreElements(); ) {
	    index++;

	    JRadioButton b = (JRadioButton)e.nextElement();
	    if( b.getModel() == group.getSelection() ) { break; }
	}
	addData( index );

	Double value = Double.parseDouble( field.getText() );
	String adjustType = comboBox.getSelectedItem().toString();
	double adjust =
	    ( adjustType.equals( "mM" ) ) ? 0.001 :
	    ( adjustType.equals( "uM" ) ) ? 0.000001 :
	    ( adjustType.equals( "nM" ) ) ? 0.000000001 :
	    0.000000000001;
	addData( value * adjust );

	addData( ( !box.isSelected() ) ? 0 : 3 );
	addData( (Integer)spinners[0].getModel().getValue() );
	addData( (Integer)spinners[1].getModel().getValue() );
	addData( panel.getFields()[1].getText() );

	String concString = Double.toString( value ) + " " + adjustType;
	addData( concString );
    }
    
    /**
     * An inner class which allows an OligoWalkWindow to be constructed from 
     * the more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * The type of nucleic acid the created window uses
         */
        private String acid;

        /**
         * Constructor
         * @param acid  the type of nucleic acid the created window uses
         */
        public Factory( String acid ) { this.acid = acid; }

        /**
         * Create an OligoWalkWindow.
         * @return  a new OligoWalkWindow
         */
        public OligoWalkWindow createWindow() {
            return new OligoWalkWindow( acid );
        }
    }
}
