/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
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
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;

import java.util.LinkedList;

/**
 * A class that draws a structure.
 * @author Jessica Reuter
 */
public class DrawStructure extends SingleSketcher {
    /**
     * The array of probability color annotation values. These values are read
     * once for each structure, and do not change.
     */
    private Color[] annotations1;

    /**
     * The array of SHAPE color annotation values. These values are read from a
     * file, and can be changed at any time depending on the file being read.
     */
    private Color[] annotations2;

    /**
     * The annotation file for probability annotations.
     */
    private String annotationFile1;

    /**
     * The annotation file for SHAPE annotations.
     */
    private String annotationFile2;

    /**
     * Render the structure clockwise (true) or counterclockwise (false)
     */
    private boolean clockwise;

    /**
     * Array of drawing data: which structure and which type
     * Structures are 1-indexed, types are as follows:
     * 0 = plain, 1 = energy, 2 = probability
     */
    protected int[] drawingData;

    /**
     * Boolean of whether this structure contains pairs.
     */
    protected boolean hasPairs;

    /**
     * The array of structure label data
     */
    private Point[] labelData;

    /**
     * The initial offset for drawing, used to create a border
     */
    protected Point offset;

    /**
     * The array of structure bounds (min/max x and min/max y)
     */
    protected int[] structureBounds;

    /**
     * The array of structure data
     */
    private Object[][] structureData;

    /**
     * Constructor
     * @param file  the name of the file this structure comes from 
     */
    public DrawStructure( String file ) {
        super( file );
        setStrandType( 1 );
        setFilters( "CT Files,ct" );

        hasPairs = false;
        clockwise = true;
        scale = 1;

        drawingData = new int[2];
        drawingData[0] = 1;
        drawingData[1] = 0;

	annotationFile1 = "";
	annotationFile2 = "";
    }

    /**
     * Get whether this structure contains pairs.
     * @return  true if this structure contains pairs, false if not.
     */
    public boolean containsPairs() { return hasPairs; }

    /**
     * Create the particular sequence of drawings on the input panel.
     */
    protected void createImagePanel() {
        imagePanel = new PrintablePanel() {
            private static final long serialVersionUID = -6100421208832478124L;

            public void paintComponent( Graphics g ) {
                // Get the graphics context and prepare the canvas
                Graphics2D g2 = (Graphics2D)g;
                g2.setRenderingHint( RenderingHints.KEY_INTERPOLATION,
                    RenderingHints.VALUE_INTERPOLATION_BILINEAR );

                g2.setColor( Color.white );
                g2.fillRect( 0, 0, getWidth(), getHeight() );
                g2.setColor( Color.black );
                g2.scale( scale, scale );

                // Create the backbone path and initialize its starting point,
                // then draw it
                GeneralPath path = new GeneralPath( Path2D.WIND_EVEN_ODD );
                Point start = (Point)structureData[0][0];
                path.moveTo( start.x, start.y );

                int length = structureData.length;
                for( int i = 0; i < length; i++ ) {
                    Point next = (Point)structureData[i][0];
                    path.lineTo( next.x, next.y );

		    if( ( i + 1 ) % 10 == 0 ) {
			int labelX = labelData[i/10].x;
			int labelY = labelData[i/10].y;
			if( ( labelX != 0 ) && ( labelY != 0 ) ) {
			    path.lineTo( labelX, labelY );
			    path.lineTo( next.x, next.y );
			}
		    } 
                }
                g2.draw( path );

                // Draw pairs
                g2.setStroke( new BasicStroke( (float)5 ) );
                for( int i = 0; i < length; i++ ) {
                    if( structureData[i].length == 3 ) {
                        Point current = (Point)structureData[i][0];
                        Point pair = (Point)structureData[i][2];

                        g2.draw( new Line2D.Double(
                            current.x, current.y, pair.x, pair.y ) );
                    }
                }
                g2.setStroke( new BasicStroke( (float)1 ) );

                // Add nucleotide and number labels
                g2.setFont( new Font( g2.getFont().getName(), 1, 16 ) );
                FontMetrics metrics = g2.getFontMetrics();
                int tall = metrics.getAscent();

                for( int i = 0; i < length; i++ ) {
                    Character nucleotide = (Character)structureData[i][1];
                    Point next = (Point)structureData[i][0];

                    int wide = metrics.charWidth( nucleotide );
                    int placeX = next.x - ( wide / 2 );
                    int placeY = next.y + ( tall / 2 );

                    Color annotationColor =
                        ( drawingData[1] == 0 ) ? Color.black :
                        ( drawingData[1] == 1 ) ? annotations1[i] :
                        annotations2[i];
                    
                    g2.setColor( Color.white );
                    g2.fillRect(
                        placeX - 1, placeY - tall - 1, wide + 2, tall + 2 );
                    g2.setColor( annotationColor );
                    g2.drawString( nucleotide.toString(), placeX, placeY );

                    if( (i+1) % 10 == 0 ) {
			int index = ( ( i + 1 ) / 10 ) - 1;
                        int labelPlaceX = labelData[index].x;
                        int labelPlaceY = labelData[index].y;

			if( !( labelPlaceX == 0 && labelPlaceY == 0 ) ) {
			    String number = Integer.toString( i + 1 );
			    int stringWidth = metrics.stringWidth( number );
			    g2.setColor( Color.white );
			    g2.fillRect(
				labelPlaceX - 1, labelPlaceY - tall - 1,
				stringWidth + 2, tall + 2 );
			    g2.setColor( Color.black );
			    g2.drawString( number, labelPlaceX, labelPlaceY );
			}
                    }
                }

                // Display the structure information on the label
                displayStructureInformation();
            }
        };
    }

    /**
     * Display the structure's information on its Imager's label.
     */
    public void displayStructureInformation() {
        RNA strand = holder.getStrand();
        String text = "STRUCTURE " + Integer.toString( drawingData[0] ) +
            " of " + Integer.toString( strand.GetStructureNumber() ) +
            "; " + strand.GetCommentString( drawingData[0] );

        LinkedList<MenuChangingWindow> list = RNAstructure.getWindowList();
        for( MenuChangingWindow window: list ) {
            boolean target =
                window instanceof Imager &&
                ((Imager)window).getSketcher() == this;

            if( target ) {
                ((Imager)window).setInfoLabel( text );
                break;
            }
        }
    }

    /**
     * Determine the position of a nucleotide, and if necessary a label.
     * @param nuc  the nucleotide to find
     * @return  the nucleotide position as a point
     */
    protected Point2D findNucleotidePosition( int nuc ) {
        RNA strand = holder.getStrand();

        int x = strand.GetNucleotideXCoordinate( nuc );
        if( !holder.isErrorStatus() ) {
            int y = strand.GetNucleotideYCoordinate( nuc );
            if( !holder.isErrorStatus() ) {
                if( x < structureBounds[0] ) { structureBounds[0] = x; }
                else if( x > structureBounds[1] ) { structureBounds[1] = x; }

                if( y < structureBounds[2] ) { structureBounds[2] = y; }
                else if( y > structureBounds[3] ) { structureBounds[3] = y; }

		if( nuc % 10 == 0 ) {
		    int labelX = strand.GetLabelXCoordinate( nuc );
		    if( !holder.isErrorStatus() ) {
			if( labelX < structureBounds[0] ) {
			    structureBounds[0] = labelX;
			} else if( labelX > structureBounds[1] ) {
			    structureBounds[1] = labelX;
			}
		    } else { preparationError = true; }

		    int labelY = strand.GetLabelYCoordinate( nuc );
		    if( !holder.isErrorStatus() ) {
			if( labelY < structureBounds[2] ) {
			    structureBounds[2] = labelY;
			} else if( labelY > structureBounds[3] ) {
			    structureBounds[3] = labelY;
			}
		    } else { preparationError = true; }

		    labelData[(nuc / 10) - 1] = new Point( labelX, labelY );
		}

                return new Point2D.Double( x, y );
            } else {
                preparationError = true;
                return null;
            }
        } else {
            preparationError = true;
            return null;
        }
    }

    /**
     * Get the current structure number.
     * @return  the current structure number
     */
    public int getCurrentStructure() { return drawingData[0]; }

    /**
     * Get the preferred image panel size.
     * @return  the preferred image panel size
     */
    protected Dimension getPreferredPanelSize() {
        int width = structureBounds[1] + offset.x + 30;
        int height = structureBounds[3] + offset.y + 30;
        return new Dimension( width, height );
    }

    /**
     * Get the current probability annotation file.
     * @return  the current probability annotation file.
     */
    public String getProbabilityAnnotationFile() { return annotationFile1; }

    /**
     * Get the current SHAPE annotation file.
     * @return  the current probability annotation file.
     */
    public String getSHAPEAnnotationFile() { return annotationFile2; }

    /**
     * Return whether the structure is rendered clockwise (true) or
     * counterclockwise (false)
     * @return  the render direction
     */
    public boolean isClockwise() { return clockwise; }

    /**
     * Prepare for structure drawing by getting proper structure information.
     */
    public void prepare() {
	// Get the RNA strand and offset.
        RNA strand = holder.getStrand();
        offset = new Point( 30, 30 );

        // Initialize the structure bounds array
        if( structureBounds == null ) {
            structureBounds = new int[4];
            structureBounds[0] = Integer.MAX_VALUE;
            structureBounds[1] = Integer.MIN_VALUE;
            structureBounds[2] = Integer.MAX_VALUE;
            structureBounds[3] = Integer.MIN_VALUE;
        }

        // Read structure data
        if( structureData == null ) {
            // Determine drawing coordinates
            int coordError =
                strand.DetermineDrawingCoordinates( 30, 30, drawingData[0] );

            if( holder.isErrorStatus( coordError ) ) {
                preparationError = true;
                return;
            }

            // Get data on nucleotide positions and angles
            int length = strand.GetSequenceLength();
            if( !holder.isErrorStatus() ) {
                structureData = new Object[length][];
		labelData = new Point[length / 10];

                // For each nucleotide, get its data
                for( int i = 1; i <= length; i++ ) {
                    Object[] row;
                    int rowLength;

                    // Check if this nucleotide is paired
                    int pair = strand.GetPair( i, drawingData[0] );
                    if( !holder.isErrorStatus() ) {
                        rowLength = ( pair == 0 ) ? 2 : 3;
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

                    // If the nucleotide is paired, find position of the pair
                    if( pair != 0 ) {
                        hasPairs = true;
                        row[2] = findNucleotidePosition( pair );
                        if( preparationError ) { break; }
                    }

		    // Add the new row of structure data.
                    structureData[i-1] = row;
                }

		// If structure is probability annotated, read annotation data,
		// as this is structure dependent.
		if( !annotationFile1.equals( "" ) ) {
		    annotations1 = AnnotationMapHolder.readProbabilityColors(
			annotationFile1, drawingData[0] );
		}
            }
        }

        // Re-translate points to create positive structure bounds
        translateCoordinates( true );
        translateCoordinates( false );

	structureBounds[1] += 30;
	structureBounds[3] += 30;

        // Set the initial canvas size
	createImagePanel();
        Dimension preferred = getPreferredPanelSize();
	
	double xScale = (double)595 / (double)preferred.width;
	double yScale = (double)575 / (double)preferred.height;
	scale = Math.min( xScale, yScale );

        Dimension scaled = new Dimension(
            (int)(preferred.width * scale), (int)(preferred.height * scale) );
        imagePanel.setPreferredSize( scaled );
    }

    /**
     * Render the structure either clockwise or counterclockwise.
     * @param clockwise  true if clockwise, false if not
     */
    public void renderClockwise( boolean clockwise ) {
        this.clockwise = clockwise;
        setStructureNumber( drawingData[0] );
	imagePanel.validate();
        imagePanel.repaint();
    }

    /**
     * Set the structure as unannotated, and repaint it.
     */
    public void setUnannotated() {
        drawingData[1] = 0;
        imagePanel.repaint();
    }

    /**
     * Set the structure as probability annotated, and repaint it.
     */
    public void setProbabilityAnnotated() {
        drawingData[1] = 1;
        imagePanel.repaint();
    }

    /**
     * Set the probability annotation file.
     * @param file  the probability annotation file to set.
     */
    public void setProbabilityAnnotationFile( String file ) {
	annotationFile1 = file;
    }

    /**
     * Set the probability annotation array.
     * @param array  the new probability annotation array
     */
    public void setProbabilityArray( Color[] array ) { annotations1 = array; }

    /**
     * Set the structure as SHAPE annotated, and repaint it.
     */
    public void setSHAPEAnnotated() {
        drawingData[1] = 2;
        imagePanel.repaint();
    }

    /**
     * Set the SHAPE annotation file.
     * @param file  the SHAPE annotation file to set.
     */
    public void setSHAPEAnnotationFile( String file ) {
	annotationFile2 = file;
    }

    /**
     * Set the SHAPE annotation array.
     * @param array  the new SHAPE annotation array
     */
    public void setSHAPEArray( Color[] array ) { annotations2 = array; }

    /**
     * Set a particular structure number to be drawn.
     * @param number  the new structure number
     */
    public void setStructureNumber( int number ) {
        double scale = getScale();
        drawingData[0] = number;
        structureBounds = null;
        structureData = null;
        prepare();
        setScale( scale );
        imagePanel.repaint();
    }

    /**
     * Translate coordinates to ensure positive structure bounds
     * @param xCoords  true if translating x coordinates, false if not
     */
    private void translateCoordinates( boolean xCoords ) {
        int index = ( xCoords ) ? 0 : 2;
	int length = structureData.length;
	int length2 = labelData.length;

        if( structureBounds[index] < 0 ) {
            if( xCoords ) { offset.x -= structureBounds[index]; }
            else { offset.y -= structureBounds[index]; }

            for( int i = 0; i < length; i++ ) {
                Point2D next2D = (Point2D)structureData[i][0];
		Point next =
		    new Point( (int)next2D.getX(), (int)next2D.getY() );
                if( xCoords ) { next.x += offset.x; }
                else { next.y += offset.y; }
                structureData[i][0] = next;

                if( structureData[i].length == 3 ) {
                    Point2D pair2D = (Point2D)structureData[i][2];
		    Point pair =
			new Point( (int)pair2D.getX(), (int)pair2D.getY() );

                    if( xCoords ) { pair.x += offset.x; }
                    else { pair.y += offset.y; }
                    structureData[i][2] = pair;
                }
            }

	    for( int i = 0; i < length2; i++ ) {
		if( !( labelData[i].x == 0 && labelData[i].y == 0 ) ) {
		    if( xCoords ) { labelData[i].x += offset.x; }
		    else { labelData[i].y += offset.y; }
		}
	    }
        }

	// If the structure is to be rendered clockwise, flip the X coordinates
	if( index == 0 && clockwise ) {
	    int adjustment = structureBounds[1] + offset.x + 30;
	    for( int i = 0; i < length; i++ ) {
		((Point)structureData[i][0]).x *= -1;
		((Point)structureData[i][0]).x += adjustment;

		if( structureData[i].length == 3 ) {
		    ((Point)structureData[i][2]).x *= -1;
		    ((Point)structureData[i][2]).x += adjustment;
		}
	    }

	    for( int i = 0; i < length2; i++ ) {
		if( !( labelData[i].x == 0&& labelData[i].y == 0 ) ) {
		    labelData[i].x *= -1;
		    labelData[i].x += adjustment;
		}
	    }
	}
    }
}
