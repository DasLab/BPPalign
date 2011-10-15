/*
 * (c) 2010  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.RenderingHints;
import java.awt.geom.AffineTransform;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;

import java.util.LinkedList;

/**
 * A class that draws a circular structure.
 * @author Jessica Reuter
 */
public class DrawStructureCircular extends DrawStructure {
    /**
     * The angle between nucleotides on the circle.
     */
    private double angle;

    /**
     * The X coordinate of the center of the circle.
     */
    private int centerX;

    /**
     * The Y coordinate of the center of the circle.
     */
    private int centerY;

    /**
     * The points at which nucleotides are found on the circle.
     */
    private Object[][] circleStructureData;

    /**
     * The radius of the circle.
     */
    private int radius;

    /**
     * Constructor
     * @param file  the name of the file this structure comes from 
     */
    public DrawStructureCircular( String file ) { super( file ); }

    /**
     * Create the image panel.
     */
    protected void createImagePanel() {
        imagePanel = new PrintablePanel() {
            private static final long serialVersionUID = -6100421234532478124L;

            public void paintComponent( Graphics g ) {
                // Get the graphics context and prepare the canvas
                Graphics2D g2 = (Graphics2D)g;
                g2.setRenderingHint( RenderingHints.KEY_INTERPOLATION,
                    RenderingHints.VALUE_INTERPOLATION_BILINEAR );

                g2.setColor( Color.white );
                g2.fillRect( 0, 0, getWidth(), getHeight() );
                g2.setColor( Color.black );

                // Draw structure.
                int length = circleStructureData.length;

                for( int i = 0; i < length; i++ ) {
                    // Create backbone path and draw it with nucleotide labels
                    Double drawAngle = (Double)circleStructureData[i][2];
                    String nucleotide = Character.toString(
                        (Character)circleStructureData[i][1] );

                    Point2D place = (Point2D)circleStructureData[i][0];
		    AffineTransform transform = new AffineTransform();
		    transform.scale( scale, scale );
		    transform.translate( place.getX() + 60,
					 place.getY() + 60 );
		    transform.rotate( drawAngle );
		    g2.setTransform( transform );
                    g2.drawString( nucleotide, 0, 0 );
		    g2.setTransform( new AffineTransform() );
		    g2.scale( scale, scale );

                    // Draw pair, if necessary.
                    if( circleStructureData[i].length == 5 ) {
			// Determine the paired point and angle.
                        Point2D pairPt = (Point2D)circleStructureData[i][3];
			Double pairAngle = (Double)circleStructureData
			    [(Integer)circleStructureData[i][4] - 1][2];

			// Determine the end points of the arc.
			double x1 =
			    place.getX() + 60 + ( Math.cos( drawAngle ) * 6 );
			double y1 =
			    place.getY() + 60 + ( Math.sin( drawAngle ) * 6 );
			double x2 =
			    pairPt.getX() + 60 + ( Math.cos( pairAngle ) * 6 );
			double y2 =
			    pairPt.getY() + 60 + ( Math.sin( pairAngle ) * 6 );

			// Calculate midpoints.
			double midX = ( x1 + x2 ) / 2.0;
			double midY = ( y1 + y2 ) / 2.0;

			// Intialize the default arc center point.
			Point2D centerPt =
			    new Point2D.Double( centerX, centerY );

			// Determine the number of bases between the pair.
			int between =
			    (Integer)circleStructureData[i][4] - ( i + 1 );
			if( between > ( length / 2 ) ) {
			    between = length - between;
			}

			// Determine gamma, the center threshold, and the
			// distance between the arc ends.
			double gamma = 0.9;
			double centerThresh =
			    ( ( (double)radius * 2.0 ) / 8.0 ) * 5.0;
			double ends = Math.sqrt(
			    Math.pow( x2 - x1, 2 ) + Math.pow( y2 - y1, 2 ) );

			// Calculate the appropriate X and Y bend for curves.
			double distance =
			    ( ( 2.0 * (double)between ) / (double)length ) *
			    (double)radius * gamma;
			double lineAngle = Math.atan2(
			    (double)centerY - midY, (double)centerX - midX );

			// Calculate the distance along the bend for curves.
			double distX = distance * Math.cos( lineAngle ) * 2.0;
			double distY = distance * Math.sin( lineAngle ) * 2.0;

			// Edit the default control point if necessary.
			if( ends < centerThresh ) {
			    centerPt = new Point2D.Double(
			        midX + distX, midY + distY );
			}

			// Draw the arc.
                        Path2D arc = new Path2D.Double();
			arc.moveTo( x1, y1 );
			arc.curveTo( x1, y1,
				     centerPt.getX(), centerPt.getY(),
				     x2, y2 );
                        g2.draw( arc );
		    }

                    // Add number labels, where necessary.
                    int index = i + 1;
                    boolean doLabel =
                        ( index % 10 == 0 ) ||
                        ( index == 1 ) ||
                        ( index == length );

		    if( doLabel ) {
			g2.setFont( g2.getFont().deriveFont( Font.BOLD ) );
                        g2.translate( place.getX() + 60, place.getY() + 60 );
			g2.rotate( drawAngle );
			g2.drawString( nucleotide, 0, 0 );
                        g2.drawString( Integer.toString( i + 1 ), 0, -12 );
			g2.setFont( g2.getFont().deriveFont( Font.PLAIN ) );
                    }
                }

                // Display the structure information on the label
                displayStructureInformation();
            }
        };
    }

    /**
     * Determine the position of a nucleotide.
     * @param nuc  the nucleotide to find
     * @return  the nucleotide position as a point
     */
    protected Point2D findNucleotidePosition( int nuc ) {
        double baseAngle = ( angle * (double)nuc ) - ( 3.1415 / 2.0 );
	double radiusD = (double)radius;

        Double baseX = (double)centerX + ( radiusD * Math.cos( baseAngle ) );
        Double baseY = (double)centerY + ( radiusD * Math.sin( baseAngle ) );

	Integer x = baseX.intValue();
	Integer y = baseY.intValue();

	if( x < structureBounds[0] ) { structureBounds[0] = x; }
	if( x > structureBounds[1] ) { structureBounds[1] = x; }
	if( y < structureBounds[2] ) { structureBounds[2] = y; }
	if( y > structureBounds[3] ) { structureBounds[3] = y; }

        return new Point2D.Double( baseX, baseY );
    }

    /**
     * Prepare for structure drawing by getting proper structure information.
     */
    public void prepare() {
        RNA strand = holder.getStrand();
        offset = new Point( 0, 0 );
        int fontSize = 12;

	int length = strand.GetSequenceLength();
	radius = ( 3 * length );
	centerX = radius;
	centerY = radius;

	structureBounds = new int[4];
	structureBounds[0] = Integer.MAX_VALUE;
	structureBounds[1] = Integer.MIN_VALUE;
	structureBounds[2] = Integer.MAX_VALUE;
	structureBounds[3] = Integer.MIN_VALUE;

	if( circleStructureData == null ) {
	    int extraBases =
		( length >= 10000 ) ? 4 :
		( length >= 1000 ) ? 3 :
		( length >= 100 ) ? 2 :
		1;

	    angle = ( 2.0 * 3.1415 ) / (double)( length + extraBases );

	    // Read structure data
	    if( !holder.isErrorStatus() ) {
		circleStructureData = new Object[length][];

		// For each nucleotide, get its data
		for( int i = 1; i <= length; i++ ) {
		    Object[] row;
		    int rowLength;

		    // Check if this nucleotide is paired
		    int pair = strand.GetPair( i, drawingData[0] );
		    if( !holder.isErrorStatus() ) {
			rowLength = ( pair > i ) ? 5 : 3;
			row = new Object[rowLength];
		    } else {
			preparationError = true;
			return;
		    }

		    // Get position of the nucleotide and the nucleotide itself
		    row[0] = findNucleotidePosition( i );
		    if( row[0] == null ) { break; }

		    char nuc = strand.GetNucleotide( i );
		    if( !holder.isErrorStatus() ) { row[1] = nuc; }
		    else {
			preparationError = true;
			break;
		    }

		    // Determine the angle of the nucleotide on the circle.
		    row[2] = angle * i;

		    // If the nucleotide is paired, find position of the pair
		    if( rowLength == 5 ) {
			hasPairs = true;
			row[3] = findNucleotidePosition( pair );
			if( preparationError ) { break; }
			row[4] = pair;
		    }

		    circleStructureData[i-1] = row;
		}
	    }

	    // Translate as necessary to ensure positive structure bounds.
	    if( structureBounds[0] < 0 ) {
		structureBounds[0] -= structureBounds[0];
		structureBounds[1] -= structureBounds[0];

		for( int i = 0; i < length; i++ ) {
		    Object[] row = circleStructureData[i];

		    Point2D.Double p1X = (Point2D.Double)row[0];
		    double p1X_xVal = p1X.getX();
		    double p1X_yVal = p1X.getY();
		    p1X.setLocation( p1X_xVal - structureBounds[0], p1X_yVal );
		    circleStructureData[i][0] = p1X;

		    if( row.length == 5 ) {
			Point2D.Double p2X = (Point2D.Double)row[3];
			double p2X_xVal = p2X.getX();
			double p2X_yVal = p2X.getY();
			p2X.setLocation(
			    p2X_xVal - structureBounds[2], p2X_yVal );
			circleStructureData[i][3] = p2X;
		    }
		}
	    }

            if( structureBounds[2] < 0 ) {
                structureBounds[2] -= structureBounds[2];
                structureBounds[3] -= structureBounds[2];

                for( int i = 0; i < length; i++ ) {
		    Object[] row = circleStructureData[i];

		    Point2D.Double p1Y = (Point2D.Double)row[0];
                    double p1Y_xVal = p1Y.getX();
                    double p1Y_yVal = p1Y.getY();
                    p1Y.setLocation( p1Y_xVal, p1Y_yVal - structureBounds[2] );
                    circleStructureData[i][0] = p1Y;

                    if( row.length == 5 ) {
			Point2D.Double p2Y = (Point2D.Double)row[3];
			double p2Y_xVal = p2Y.getX();
                        double p2Y_yVal = p2Y.getY();
                        p2Y.setLocation(
			    p2Y_xVal, p2Y_yVal - structureBounds[2] );
                        circleStructureData[i][3] = p2Y;
                    }
		}
            }

	    // Set the initial canvas size
	    createImagePanel();
	    Dimension preferred = getPreferredPanelSize();
	    double xScale = (double)595 / (double)preferred.width;
	    double yScale = (double)575 / (double)preferred.height;
	    scale = Math.min( xScale, yScale );

	    Dimension scaled = new Dimension(
                (int)(preferred.width * scale),
		(int)(preferred.height * scale) );
	    imagePanel.setPreferredSize( scaled );
	}
    }
}
