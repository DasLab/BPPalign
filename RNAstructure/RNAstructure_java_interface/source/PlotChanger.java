/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JComponent;
import javax.swing.JTextField;

/**
 * A class which handles actions that modify a dot plot.
 * @author Jessica Reuter
 */
public class PlotChanger extends ActionBase {
    /**
     * The type of action being done
     */
    private int type;

    /**
     * Private constructor
     * @param type  the type of action being done
     */
    private PlotChanger( int type ) { this.type = type; }

    /**
     * Do a plot modifying action.
     */
    protected void act() {
        DrawDotPlot plotter =
            (DrawDotPlot)(((Imager)RNAstructure.getCurrentWindow())
	    .getSketcher());

        // Change colors in a plot
        if( type == 1 ) {
            String value = ((SpinnerWindow)((JComponent)event.getSource())
                .getTopLevelAncestor()).getSpinner().getValue().toString();
            double number = Double.parseDouble( value );

            plotter.changeColors( (int)number );
            ((Imager)RNAstructure.getCurrentWindow())
                .repaintLegend( plotter.getLegend() );
        }

        // Change the range of a plot
        else if( type == 2 ) {
            JTextField[] fields =
                ((PlotRangeWindow)((JComponent)event.getSource())
                .getTopLevelAncestor()).getFields();
            float newMin = Float.parseFloat( fields[0].getText() );
            float newMax = Float.parseFloat( fields[1].getText() );

            plotter.setVisibleRange( newMin, newMax );
            plotter.getImagePanel().repaint();
            plotter.createLegend( plotter.getEntries() );
            ((Imager)RNAstructure.getCurrentWindow())
                .repaintLegend( plotter.getLegend() );
        }

        // Validate the plot
        RNAstructure.getCurrentWindow().validate();
        RNAstructure.getCurrentWindow().repaint();
    }

    /**
     * Create an action that modifies the number of colors in a plot
     * @return  the color changing PlotChanger
     */
    public static PlotChanger makeColorAction() {
        return new PlotChanger( 1 );
    }

    /**
     * Create an action that modifies the range of a plot
     * @return  the range modifying PlotChanger
     */
    public static PlotChanger setRangeAction() { return new PlotChanger( 2 ); }
}
