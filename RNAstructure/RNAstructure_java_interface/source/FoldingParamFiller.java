/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JComponent;

/**
 * A class that fills a particular TextInputPanel with folding parameters.
 * @author Jessica Reuter
 */
public class FoldingParamFiller extends DataSelector {
    /**
     * The data needed to fill parameters.
     */
    private Object[] data;

    /**
     * Constructor
     * @param input  the data needed to fill parameters
     */
    private FoldingParamFiller( boolean constant, Object[] input ) {
        int length = input.length;
        data = new Object[length + 1];
        data[0] = constant;
        for( int i = 0; i < length; i++ ) { data[i+1] = input[i]; }
    }

    /**
     * Set parameters in a FoldWindow.
     */
    protected void selectInput() {
        MenuChangingWindow top =
            (MenuChangingWindow)(((JComponent)event.getSource())
            .getTopLevelAncestor());
	
	boolean fill =
	    top instanceof DynalignWindow ||
	    top instanceof FoldWindow;

        if( fill ) {
            boolean constant = (Boolean)data[0];            
            int length = data.length;

	    TextInputPanel panel = null;
	    if( top instanceof DynalignWindow ) {
		panel = ((DynalignWindow)top).getParameterPanel();
	    } else if( top instanceof FoldWindow ) {
		panel = ((FoldWindow)top).getParameterPanel();
	    }

            // If parameters are constant, use a simple for loop to set
            if( constant ) {
                for( int i = 1; i < length; i++ ) {
                    panel.getFields()[i-1].setText( (String)data[i] );
                }
            }

            // If parameters are variable, use more complex logic.
            else {
                for( int i = 1; i < length; i++ ) {
                    // There may be some constant variables in the
                    // variable group, so if the next object is a
                    // constant value, set it as such
                    if( data[i] instanceof String ) {
                        panel.getFields()[i-1].setText( (String)data[i] );
                    }

                    // If object is array of lengths and values, use 
                    // them to set a field
                    else {
                        // Get the nucleic acid strand length.
                        int acids = RNAstructure.getCurrentStrand()
			    .GetSequenceLength();

                        // All lengths and values are dealt with as
                        // thresholds, a query must be greater than a
                        // specific value to meet a category.
                        // For each cutoff value, check how it compares
                        // to the sequence length, and if a cutoff
                        // matches, set the appropriate value
                        Object[][] values = (Object[][])data[i];
                        int choices = values.length;

                        for( int j = 0; j < choices; j++ ) {
                            double next =
                                ((Number)values[j][0]).doubleValue();

                            if( acids > next ) {
                                panel.getFields()[i-1]
                                    .setText( values[j][1].toString() );
                                j = choices + 1;
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Create an action that allows for setting constant structure parameters.
     * @param data  the data itself, unformatted
     * @return  the constant parameters listener
     */
    public static FoldingParamFiller setConstant( Object... data ) {
        return new FoldingParamFiller( true, data );
    }

    /**
     * Create an action that allows for setting variable structure parameters.
     * @param data  the data itself, unformatted
     * @return the variable parameters listener
     */
    public static FoldingParamFiller setVariable( Object... data ) {
        return new FoldingParamFiller( false, data );
    }
}
