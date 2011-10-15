/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JTextField;

/**
 * A class holding a window that modifies the range of a dot plot.
 * @author Jessica Reuter
 */
public class PlotRangeWindow extends DialogWindow {
    /**
     * The main input panel
     */
    private TextInputPanel panel;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -8616494515994365054L;

    /**
     * Constructor
     */
    public PlotRangeWindow() {
        // Get the plot bounds and currently visible range
        DrawDotPlot plotter =
	    (DrawDotPlot)((Imager)RNAstructure.getCurrentWindow())
            .getSketcher();
        float[] bounds = plotter.getPlotBounds();
        float[] range = plotter.getVisibleRange();

        // Create the main input panel
        panel = new TextInputPanel();
        panel.create( 2, 2, 20, "Min", "Max" );

        JTextField[] fields = panel.getFields();
        fields[0].setText( Float.toString( range[0] ) );
        fields[1].setText( Float.toString( range[1] ) );

        // Create the middle reset button panel
        ChainedAction reset =
            new ChainedAction(
                new SetInputField( fields[0], Float.toString( bounds[0] ) ),
                new SetInputField( fields[1], Float.toString( bounds[1] ) ) );

        ButtonPanel buttonPanel = new ButtonPanel( 1 );
        buttonPanel.setAction( reset, 0, "Reset" );
        buttonPanel.addButtons();

        // Create the bottom button panel
        ChainedAction changer =
            new ChainedAction(
                PlotChanger.setRangeAction(),
                new CloseAction() );

        ButtonPanel buttonPanel2 = new ButtonPanel( 2 );
        buttonPanel2.setAction( changer, 0, "OK" );
        buttonPanel2.setAction( new CloseAction(), 1, "Cancel" );
        buttonPanel2.addButtons();

        // Place components in their proper places
        placeComponent( 0, 0, panel );
        placeComponent( 0, 1, buttonPanel );
        placeComponent( 0, 2, buttonPanel2 );
    }

    /**
     * Get the input fields in this window.
     * @return  the input fields
     */
    public JTextField[] getFields() { return panel.getFields(); }

    /**
     * An inner class that allows a PlotRangeWindow to be constructed from the
     * more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * Create a PlotRangeWindow.
         */
        public Object createWindow() {
            return new PlotRangeWindow();
        }
    }
}
