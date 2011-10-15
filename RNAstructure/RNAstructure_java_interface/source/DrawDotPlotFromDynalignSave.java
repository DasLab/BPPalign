/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.io.File;

/**
 * A class that draws one dot plot from a Dynalign save file. Two instances of
 * this class must be created to draw both plots from a save file.
 * @author Jessica Reuter
 */
public class DrawDotPlotFromDynalignSave extends DrawDotPlot {
    /**
     * The strand to draw, first or second.
     */
    private int strand;

    /**
     * Constructor
     * @param file  the file being analyzed
     * @param entries  the number of entries in the dot plot
     * @param strand  the strand to draw
     */
    private DrawDotPlotFromDynalignSave( String file, int entries,
					 int strand ) {
        super( file, 0, entries );
	this.strand = strand;
    }

    /**
     * Check to make sure the file is valid
     */
    public void create() {
        if( file.equals( "" ) ) { preparationError = true; }
    }

    /**
     * Create an action that draws a dot plot for strand 1 from a
     * Dynalign save file.
     * @param file  the Dynalign save file
     * @return  the dot plot strand 1 drawing action
     */
    public static DrawDotPlotFromDynalignSave drawFirstPlot( String file ) {
        return new DrawDotPlotFromDynalignSave( file, 5, 1 );
    }

    /**
     * Create an action that draws a dot plot for strand 2 from a
     * Dynalign save file.
     * @param file  the Dynalign save file
     * @return  the dot plot strand 2 drawing action
     */
    public static DrawDotPlotFromDynalignSave drawSecondPlot( String file ) {
	return new DrawDotPlotFromDynalignSave( file, 5, 2 );
    }

    /**
     * Get the strand number handled in this Dynalign plot.
     * @return  the strand number.
     */
    public int getStrandNumber() { return strand; }

    /**
     * Read in the dots array from the file.
     */
    public void prepare() {
        try {
	    Dynalign_object dynalign = new Dynalign_object( file );

	    RNA acids = ( strand == 1 ) ?
		dynalign.GetRNA1() : dynalign.GetRNA2();

            sequenceLength = acids.GetSequenceLength();
	    if( acids.GetErrorCode() != 0 ) {
		preparationError = true;
		return;
	    }

            dots = new float[sequenceLength][];
            for( int i = sequenceLength; i > 0; i-- ) {
                dots[sequenceLength - i] = new float[i];
                for( int j = 0; j < i; j++ ) {
                    dots[sequenceLength - i][j] = Float.POSITIVE_INFINITY;
                }
            }

            setDivider( "Free Energy (kcal/mol)" );
            setPlotType( 0 );

            float min = Float.MAX_VALUE, max = Float.MAX_VALUE * -1;
	    for( int i = 1; i <= sequenceLength; i++ ) {
		float[] row = new float[sequenceLength - (i-1)];

		for( int j = i; j <= sequenceLength; j++ ) {
		    double raw = dynalign.GetBestPairEnergy( strand, i, j );

		    if( dynalign.GetErrorCode() != 0 ) {
			RNAstructureInfoDialog.error( dynalign.GetErrorMessage(
			    dynalign.GetErrorCode() ) );
			preparationError = true;
			return;
		    }

		    if( !preparationError ) {
			float value = Float.valueOf( Double.toString( raw ) );

			boolean infinite =
			    value == Float.POSITIVE_INFINITY ||
			    value == Float.NEGATIVE_INFINITY;
			
			if( infinite || value > 0 ) {
			    value = Float.POSITIVE_INFINITY;
			}
			else if( value > max ) { max = value; }
			else if( value < min ) { min = value; }

			row[j-i] = value;
		    }
		}

		dots[i-1] = row;
	    }	    

            scale = (double)595 / (double)(getPreferredPanelSize().width);
            setLegendAndBounds( min, max );

	    // Edit the file to distinguish between sequences
	    String identifier = "";
	    if( file.contains( "_" ) ) {
		String name = new File( file ).getName();
		if( strand == 1 ) {
		    identifier = name.substring( 0, name.indexOf( "_" ) );
		} else {
		    identifier = name.substring(
			name.indexOf( "_" ) + 1, name.indexOf( "." ) );
		}
	    } else { identifier = Integer.toString( strand ); }

	    file = file.concat( "  Sequence: " ).concat( identifier );
        } catch( Exception e ) {
	    RNAstructureInfoDialog.error( "Error preparing dot plot." );
	    preparationError = true;
	}
    }
}
