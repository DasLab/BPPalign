/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.text.DecimalFormat;

/**
 * A class that draws a dot plot.
 * @author Jessica Reuter
 */
public class DrawDotPlot extends SingleSketcher {
    /**
     * The array of cutoffs where dots can be.
     * Column 0 is the start of each column in the plot.
     * Column 1 is the end of each column in the plot.
     */
    private short[][] cutoffs;

    /**
     * The divider in the legend that tells what type of values are described.
     */
    private String divider;

    /**
     * The entries on the scale of the dot plot.
     */
    private final int division = 10;

    /**
     * The array of dot values.
     */
    protected float[][] dots;

    /**
     * The number of categories in the dot plot
     */
    private int entries;

    /**
     * The decimal format used to format given values
     */
    private final DecimalFormat formatter = new DecimalFormat( "0.00" );

    /**
     * The minimum and maximum values in the dot plot.
     */
    private float[] plotBounds;

    /**
     * The array of thresholds at which dot values separate.
     */
    private double[] thresholds;

    /**
     * The translation point of this plot.
     */
    private final Point translatePoint =
        new Point( 50, 50 + ( division * 5 ) );

    /**
     * The type of dot plot being drawn
     * 0 = energy plot
     * 1 = probability plot
     */
    private int type;

    /**
     * The visible range of the dots on the plot. This range may or may not be
     * the same as the maximum and minimum values in the dot plot
     */
    private float[] visibleRange;

    /**
     * Constructor
     * @param file  the name of the file this plot comes from
     * @param type  the type of dot plot to be drawn
     * @param entries  the number of entries in the dot plot 
     */
    protected DrawDotPlot( String file, int type, int entries ) {
        super( file );
        this.type = type;
        this.entries = entries;

        imagePanel = new PrintablePanel() {
            private static final long serialVersionUID = 1960856843396790364L;

            public void paintComponent( Graphics g ) {
                // Create the graphics context
                Graphics2D g2 = (Graphics2D)g;
                g2.setRenderingHint( RenderingHints.KEY_INTERPOLATION,
                    RenderingHints.VALUE_INTERPOLATION_BILINEAR );

                // Paint the canvas white
                g2.setColor( Color.white );
                g2.fillRect( 0, 0, getWidth(), getHeight() );
                g2.scale( scale, scale );

                // Translate and draw the plot
                g2.translate( translatePoint.x, translatePoint.y );
                drawFramework( g2 );
                createDots( g2 );
            }
        };
    }

    /**
     * Redraw the plot, with a different number of colors.
     * @param colors  the new number of colors to use
     */
    public void changeColors( int colors ) {
        entries = colors;
        createLegend( entries );
    }

    /**
     * Create the dots. The method, for efficiency, only draws the dots visible
     * in the current viewing area.
     * @param g2  the graphics context being drawn with
     */
    private void createDots( Graphics2D g2 ) {
        Rectangle visible = imagePanel.getVisibleRect();

        // Adjust the translation for the scale
        Point adjusted = new Point( (int)(translatePoint.x * scale),
            (int)(translatePoint.y * scale ) );

        int rows = dots.length;
        for( int i = 1; i <= rows; i++ ) {
            int rowLength = dots[i-1].length;

            for( int j = 1; j <= rowLength; j++ ) {
                // Get the raw x,y value and position of the point
                float value = dots[i-1][j-1];

                int x = ( j * 5 ) + ( i * 5 ) - 7;
                int y = ( i * 5 ) - 2;

                // Adjust for scaling if necessary
                int adjustedX = (int)( x * scale ) + adjusted.x;
                int adjustedY = (int)( y * scale ) + adjusted.y;

                double adjustedXend = adjustedX + 4*scale;
                double adjustedYend = adjustedY + 4*scale;

                boolean drawPoint =
                    ( ( visible.contains( adjustedX, adjustedY ) ) ||
                      ( visible.contains( adjustedXend, adjustedYend ) ) ) &&
                    ( value != Float.POSITIVE_INFINITY ) &&
                    ( getDotIsVisible( value ) );

                // Draw a point if it's visible
                if( drawPoint ) {
                    for( int k = 0; k < entries; k++ ) {
                        if( value < thresholds[k] ) {
                            g2.setColor( (Color)legend[k][1] );
                            g2.fillRect( x, y, 4, 4 );
                            break;
                        }
                    }
                }
            }
        }

        // Initialize the cutoffs array if it hasn't already been created
        // (create it on the first paint for use later)
        if( cutoffs == null ) {
            cutoffs = new short[rows][];
            for( int i = 0; i < rows; i++ ) { cutoffs[i] = new short[2]; }
            refreshCutoffs();
        }
    }

    /**
     * Create a legend.
     * @param entries  the number of entries in the legend
     */
    public void createLegend( int entries ) {
        this.entries = entries;
        legend = new Object[entries][2];

        double min = visibleRange[0];
        double max = visibleRange[1];

        // Set colors for the legend sections
        boolean odd = ( entries % 2 ) != 0;

        int last = entries - 1;
        int spaces = ( odd ) ? ( last / 2 ) - 1 : ( entries - 2 ) / 2;
        double split = 1 / (double)( spaces + 1 );

        legend[0][1] = Color.red;
        legend[entries-1][1] = Color.blue;
        if( odd ) { legend[entries/2][1] = Color.green; }

        for( int i = 1; i <= spaces; i++ ) {
            double next = split * i;

            int hue1 = (int)( ( 1 - next ) * 255 );
            int hue2 = (int)( next * 255 );

            Color shade = new Color( hue1, hue2, 0 );
            Color inverse = new Color( 0, hue2, hue1 );

            legend[i][1] = shade;
            legend[last - i][1] = inverse;
        }

        // Set text for legend sections
        double increment = ( max - min ) / entries;

        for( int i = 0; i < entries; i++ ) {
            double multi = i * increment;
            double current = min + multi;
            double next = current + increment;

            String varies = ( i != entries - 1 ) ? " < " : " <= "; 
            legend[i][0] =
                format( current ) + " <= " + divider + varies + format( next );
        }
    }

    /**
     * Draw an energy dot plot.
     * @param file  the file name of the plot to draw
     * @param divs  the number of entries in the dot plot
     */
    public static DrawDotPlot drawEnergyPlot( String file, int divs ) {
        DrawDotPlot action = new DrawDotPlot( file, 0, divs );
        action.setStrandType( 4 );
        action.setFilters( "Save Files,sav" );
        action.setDivider( "Free Energy (kcal/mol)" );
        return action;
    }

    /**
     * Draw the grid and scale for the plot.
     * @param size  the maximum size of the plot
     * @param g2  the graphics context being drawn on
     */
    private void drawFramework( Graphics2D g2 ) {
        int newSize = getSequenceSize() * 5;
        int newDivision = division * 5;
        int leftover = newSize % newDivision;

        g2.setColor( Color.black );
        FontMetrics metrics = g2.getFontMetrics();

        g2.drawLine( 0, 0, newSize, 0 );
        g2.drawLine( newSize, 0, newSize, newSize );
        g2.drawLine( 0, 0, newSize, newSize );

        int counter = 0;
        while( counter <= newSize ) {
            int point = newSize - counter;
            if( leftover != 0 ) { point -= leftover; }

            g2.drawLine( point, point, newSize + newDivision, point );
            g2.drawString( Integer.toString( point / 5 ),
                newSize + 5, point + 15 );

            g2.drawLine( point, point, point, 0 - newDivision );
            g2.rotate( Math.toRadians( 90 ) );

            int number2 = newSize - point;
            if( leftover != 0 ) { number2 -= leftover; }
            String numberString = Integer.toString( number2 / 5 );
            int width = metrics.stringWidth( numberString );

            g2.drawString( numberString, -width - 5, -counter - 5 );
            g2.rotate( Math.toRadians( -90 ) );

            counter += newDivision;
        }

        if( leftover != 0 ) {
            String sizeString = Integer.toString( getSequenceSize() );
            g2.drawLine( newSize, newSize, newSize + newDivision, newSize );
            g2.drawString( sizeString, newSize + 5, newSize + 15 );

            g2.drawLine( newSize, newSize, newSize, 0 - newDivision );
            g2.rotate( Math.toRadians( 90 ) );
            int width = metrics.stringWidth( sizeString );

            g2.drawString( sizeString, -width - 5, -newSize - 5 );
            g2.rotate( Math.toRadians( -90 ) );
        }

        thresholds = new double[entries];
        for( int i = 0; i < entries; i++ ) {
            String legendString = legend[i][0].toString();
            int end = legendString.lastIndexOf( ' ' ) + 1;
            String threshold = legendString.substring( end );
            thresholds[i] = Double.parseDouble( threshold );
        }
    }

    /**
     * Draw a probability dot plot.
     * @param file  the file name of the plot to draw
     * @param divs  the number of entries in the dot plot
     */
    public static DrawDotPlot drawProbabilityPlot( String file, int divs ) {
        DrawDotPlot action = new DrawDotPlot( file, 1, divs );
        action.setStrandType( 3 );
        action.setFilters( "Partition Function Save Files,pfs" );
        action.setDivider( "-log10(BP Probability)" );
        return action;
    }

    /**
     * Format an increment value to two decimal places only.
     * @param value  the value to format
     * @return  the value formatted as a string
     */
    private String format( double value ) {
        String formatted = formatter.format( value );
        if( formatted.equals( "-0.00" ) ) { formatted = "0.00"; }
        return formatted;
    }

    /**
     * Get the array of drawn dot row and column cutoffs.
     * @return  the array of cutoffs
     */
    public short[][] getCutoffs() { return cutoffs; }

    /**
     * Get whether a particular dot is visible.
     * @param value  the value of the dot
     * @return  true if the dot is visible, false if not
     */
    public boolean getDotIsVisible( float value ) {
        return ( value >= visibleRange[0] && value <= visibleRange[1] );
    }

    /**
     * Get the value for a particular dot.
     * @param x  the x value of the dot, one-indexed
     * @param y  the y value of the dot, one-indexed
     * @return  the value of the dot
     */
    public float getDotValue( int x, int y ) { return dots[y-1][x-y]; }

    /**
     * Get the number of entries in the plot.
     * @return  the number of entries in the plot
     */
    public int getEntries() { return entries; }

    /**
     * Get the bounds of the plot.
     * @return  the array of the bound
     */
    public float[] getPlotBounds() { return plotBounds; }

    /**
     * Get the dot plot type.
     * @return  the dot plot type
     */
    public int getPlotType() { return type; }

    /**
     * Get the preferred image panel size at 100% zoom.
     * @return  the preferred image panel size
     */
    protected Dimension getPreferredPanelSize() {
        int adjustedDivision = division * 5;
        int adjust = (getSequenceSize() * 5) + 100 + adjustedDivision;
        return new Dimension( adjust, adjust );
    }

    /**
     * Get the translation point.
     * @return  the translation point
     */
    public Point getTranslationPoint() { return translatePoint; }

    /**
     * Get the visible range of this plot.
     * @return  the array of the visible range
     */
    public float[] getVisibleRange() { return visibleRange; }

    /**
     * Gather all data for drawing a dot plot.
     */
    public void prepare() {
        // Get strand
        RNA strand = holder.getStrand();
        preparationError = false;

        // Get the sequence length
        int length = strand.GetSequenceLength();
        if( holder.isErrorStatus() ) {
            preparationError = true;
            return;
        }

        // Initialize some graphics and stuff that needs to be refreshed with
        // each plot
        dots = new float[length][];
        scale = (double)595 / (double)(getPreferredPanelSize().width);

        // Calculate the values for each dot
        float min = Float.MAX_VALUE, max = Float.MAX_VALUE * -1;

        for( int i = 1; i <= length; i++ ) {
            float[] row = new float[length - (i-1)];

            for( int j = i; j <= length; j++ ) {
                double raw =
                    ( type == 0 ) ? strand.GetPairEnergy( i, j ) :
                    ( type == 1 ) ? strand.GetPairProbability( i, j ) :
                    0.0;

                if( holder.isErrorStatus() ) {
                    preparationError = true;
                    return;
                }

                if( !( preparationError = holder.isErrorStatus() ) ) {
                    if( type == 1 ) { raw = -Math.log10( raw ); }
                    float value = Float.valueOf( Double.toString( raw ) );

                    boolean infinite =
                        value == Float.POSITIVE_INFINITY ||
                        value == Float.NEGATIVE_INFINITY;

                    if( value > max ) {;
                        if( type == 0 && value > 0.0 ) {
                            value = Float.POSITIVE_INFINITY;
                        } else if( infinite ) {
                            value = Float.POSITIVE_INFINITY;
                        } else { max = value; }
                    } else if( value < min ) { min = value; }

                    row[j-i] = value;
                } else break;
            }

            dots[i-1] = row;
        }

        setLegendAndBounds( min, max );
    }

    /**
     * Create the array of cutoffs for the first time, or refresh its values
     * when the zooming scale changes.
     */
    public void refreshCutoffs() {
        Point adjusted = new Point( (int)(translatePoint.x * scale),
            (int)(translatePoint.y * scale ) );
        int rows = dots.length;

        for( int i = 1; i <= rows; i++ ) {
            int rowLength = dots[i-1].length;

            for( int j = 1; j <= rowLength; j++ ) {
                int x = ( j * 5 ) + ( i * 5 ) - 7;
                int y = ( i * 5 ) - 2;

                int adjustedX = (int)( x * scale ) + adjusted.x;
                int adjustedY = (int)( y * scale ) + adjusted.y;

                if( i == 1 ) { cutoffs[j-1][0] = (short)(adjustedX - 1); }
                cutoffs[i-1][1] = (short)(adjustedY - 1);
            }
        }
    }

    /**
     * Set the divider for this plot.
     * @param divider  the divider to set
     */
    public void setDivider( String divider ) { this.divider = divider; }

    /**
     * Set the legend and plot bounds.
     * @param min  the minimum value of the plot
     * @param max  the maximum value of the plot
     */
    protected void setLegendAndBounds( float min, float max ) {
        // Save the plot bounds and visible range
        plotBounds = new float[2];
        plotBounds[0] = min;
        plotBounds[1] = max;

        visibleRange = new float[2];
        visibleRange[0] = min;
        visibleRange[1] = max;

        // Create the legend and save the thresholds for dots
        createLegend( entries );
    }

    /**
     * Set the plot type.
     * @param type  the plot type
     */
    protected void setPlotType( int type ) { this.type = type; }

    /**
     * Set the visible range.
     * @param min  the new minimum
     * @param max  the new maximum
     */
    public void setVisibleRange( float min, float max ) {
        visibleRange[0] = min;
        visibleRange[1] = max;
    }
}