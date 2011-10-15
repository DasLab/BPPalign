/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.Color;
import java.awt.Point;
import java.awt.Robot;

/**
 * A class that displays data about a dot in a plot when the mouse is clicked.
 * @author Jessica Reuter
 */
public class DotAnalyzer extends MouseActionBase {
    /**
     * A constant that holds the RGB value of the color black
     */
    private final int BLACK = Color.black.getRGB();

    /**
     * A constant that holds the RGB value of the color white
     */
    private final int WHITE = Color.white.getRGB();

    /**
     * The sketcher that originally drew the dot plot being analyzed
     */
    private DrawDotPlot plotter;

    /**
     * Constructor
     * @param plotter  the sketcher that originally drew the dot plot
     */
    public DotAnalyzer( DrawDotPlot plotter ) { this.plotter = plotter; }

    /**
     * Display information about the plot when a dot is clicked.
     */
    public void clickMouse() {
        try {
            Point click = event.getLocationOnScreen();
            int color = new Robot().getPixelColor( click.x, click.y ).getRGB();

            if( color != BLACK && color != WHITE ) {
                Point point = event.getPoint();
                int ptX = point.x, ptY = point.y;

                // Determine which point was clicked
                short[][] cutoffs = plotter.getCutoffs();
                int length = cutoffs.length;
                int x = 1, y = 1;

                for( int i = 1; i <= length; i++ ) { 
                    if( ptX >= cutoffs[i-1][0] ) { x = i; }
                }

                for( int i = 1; i <= length; i++ ) { 
                    if( ptY >= cutoffs[i-1][1] ) { y = i; }
                }

                // Get the plot type values and dot value
                int type = plotter.getPlotType();
                double value = plotter.getDotValue( x, y );

                String valueType =
                    ( type == 0 ) ? "Free Energy (kcal/mol)" :
                    ( type == 1 ) ? "-log10(BP Probability)" :
                    null;

                char nuc1 = '?', nuc2 = '?';
                if( plotter.getSequenceArray() != null ) {
                    nuc1 = plotter.getSequenceArray()[y-1];
                    nuc2 = plotter.getSequenceArray()[x-1];
                }

                String dataString = Integer.toString( y ) + "(" + nuc1 +
                    ") - " + Integer.toString( x ) + "(" + nuc2 + "); " +
                    valueType + " = " + Double.toString( value );
                ((Imager)RNAstructure.getCurrentWindow())
		    .setInfoLabel( dataString );
            } else {
                ((Imager)RNAstructure.getCurrentWindow())
		    .setInfoLabel( "DOT PLOT" );
            }
        } catch( Exception ex ) {
            RNAstructureInfoDialog.error( "Error retrieving dot data." );
        }
    }
}
