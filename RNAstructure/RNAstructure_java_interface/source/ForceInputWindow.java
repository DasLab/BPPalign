/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.GridLayout;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

/**
 * A class that holds a window which forces constraints on either a sequence or
 * an alignment, as the case may be.
 * @author Jessica Reuter
 */
public class ForceInputWindow extends DialogWindow {
    /**
     * A TextInputPanel for Dynalign alignment constraints
     */
    private static TextInputPanel dynalignPanel;

    /**
     * An array of text fields for the complex constraints
     */
    private static JTextField[] fields;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -9141513320939555143L;

    /**
     * A single text field for the simple constraints
     */
    private static JTextField simpleField;

    /**
     * Constructor
     * @param constraint  the constraint this window is getting input for
     */
    public ForceInputWindow( String constraint ) {
        setTitle( constraint );

        // Declare variables that may be needed, depending on the constraint
        TextInputPanel panel = null;
        JPanel gridPanel = null;
        JPanel simplePanel = null;

        fields = null;
        simpleField = null;

        boolean complex = false;

        // If constraint deals with forcing an alignment, create a top panel
        if( constraint.equals( "Force Alignment" ) ) {
            panel = new TextInputPanel();
            String[] alignLabels = {
                "Nucleotide in Sequence 1:",
                "Nucleotide in Sequence 2:"
            };

            panel.create( 2, 2, 10, alignLabels );
            dynalignPanel = panel;
        }

        // Otherwise, if constraint deals with sequence, input is different
        else {
            complex = ( constraint.equals( "Force Pair" ) || 
                constraint.equals( "Prohibit Base Pairs" ) );

            // If multiple pieces of data are needed, the input is complex
            if( complex ) {
                gridPanel = new JPanel( new GridLayout() );
                fields = new JTextField[3];

                JPanel miniPanel;
                JLabel label;
                JTextField field;

                String[] labels = { "Base 1", "Base 2", "Helix Length" };
                for( int i = 0; i <= 2; i++ ) {
                    miniPanel = new JPanel( new BorderLayout() );

                    label = new JLabel( labels[i] );
                    label.setHorizontalAlignment( JLabel.CENTER );
                    miniPanel.add( label, BorderLayout.NORTH );

                    field = new JTextField( "0" );
                    field.setColumns( 5 );
                    fields[i] = field;
                    miniPanel.add( field );

                    miniPanel.setBorder(
                        BorderFactory.createEmptyBorder( 5, 20, 5, 20 ) );

                    gridPanel.add( miniPanel );
                }
            }

            // If only one piece of data is needed, the input is simple
            else {
                simplePanel = new JPanel( new BorderLayout() );
                JLabel simpleLabel = new JLabel( "Base Number" );
                simpleLabel.setHorizontalAlignment( JLabel.CENTER );
                simplePanel.add( simpleLabel, BorderLayout.NORTH );

                simpleField = new JTextField( "0" );
                simpleField.setColumns( 5 );
                simplePanel.add( simpleField );
            }
        }

        // The bottom buttons of all force windows are identical
        // Create the top button panel
        ButtonPanel buttonPanel = new ButtonPanel( 2 );
        buttonPanel.setAction( new SetForcedInput( constraint ), 0, "OK" );
        buttonPanel.setAction( new CloseAction(), 1, "Cancel" );
        buttonPanel.addButtons();

        // Create the lowermost button
        JButton closeButton = new JButton( "OK and Close" );
        closeButton.addActionListener(
            new ChainedAction(
                new SetForcedInput( constraint ),
                new CloseAction() ) );
        closeButton.setPreferredSize( new Dimension( 125, 30 ) );

        JPanel closePanel = new JPanel();
        closePanel.add( closeButton );
        closePanel.setBorder( BorderFactory.createEmptyBorder( 5, 5, 5, 5 ) );

        // Add all components in their proper places
        setFill( 2 );
        if( complex ) { placeComponent( 0, 0, gridPanel ); }
        else {
            if( panel != null ) { placeComponent( 0, 0, panel ); }
            else { placeComponent( 0, 0, simplePanel ); }
        }

        placeComponent( 0, 1, buttonPanel );
        placeComponent( 0, 2, closePanel );
    }

    /**
     * Get the input for complex constraints.
     * @return an array of the numerical constraints
     */
    public static Integer[] getComplexConstraints() {
        Integer[] constraints = new Integer[3];
        constraints[0] = Integer.parseInt( fields[0].getText() );
        constraints[1] = Integer.parseInt( fields[1].getText() );
        constraints[2] = Integer.parseInt( fields[2].getText() );
        return constraints;
    }

    /**
     * Get the input for the Dynalign alignment constraints.
     * @return an array of the numerical constraints
     */
    public static Integer[] getDynalignConstraints() {
        Integer[] constraints = new Integer[2];
        JTextField[] fields = dynalignPanel.getFields();
        constraints[0] = Integer.parseInt( fields[0].getText() );
        constraints[1] = Integer.parseInt( fields[1].getText() );
        return constraints;
    }

    /**
     * Get the input for simple constraints.
     * @return the numerical constraint
     */
    public static Integer getSimpleConstraint() {
        return Integer.parseInt( simpleField.getText() );
    }

    /**
     * Reset the text field or fields to their initial value.
     */
    public static void reset() {
        if( simpleField != null ) { simpleField.setText( "0" ); }

        if( fields != null ) {
            fields[0].setText( "0" );
            fields[1].setText( "0" );
            fields[2].setText( "0" );
        }
    }

    /**
     * An inner class which allows a ForceInputWindow to be constructed from
     * the more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * The constraint type that is being set in the created window
         */
        private String constraint;

        /**
         * Constructor
         * @param type  the type of constraint
         */
        public Factory( String type ) { constraint = type; }

        /**
         * Create a ForceInputWindow.
         * @return  a new ForceInputWindow
         */
        public ForceInputWindow createWindow() {
            return new ForceInputWindow( constraint );
        }
    }
}