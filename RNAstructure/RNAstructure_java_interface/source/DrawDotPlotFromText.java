/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.io.BufferedReader;
import java.io.FileReader;

/**
 * A class that draws a dot plot from a text file, without the use of the RNA
 * class for holding data.
 * @author Jessica Reuter
 */
public class DrawDotPlotFromText extends DrawDotPlot {
    /**
     * Constructor
     * @param entries  the number of entries in the dot plot
     */
    private DrawDotPlotFromText( int entries ) {
        super( "", -1, entries );
    }

    /**
     * Get and set the file name.
     */
    public void create() {
        FileSelector selector =
            new FileSelector( true, "open", "Dot Plot Files,dp" );
        String file = selector.showSelector();

        if( !file.equals( "" ) ) { setFile( file ); }
        else { preparationError = true; }
    }

    /**
     * Create an action that draws a dot plot from a text file.
     * @return  the dot plot drawing action
     */
    public static DrawDotPlotFromText draw() {
        return new DrawDotPlotFromText( 5 );
    }

    /**
     * Read in the dots array from the file.
     */
    public void prepare() {
        try {
            BufferedReader reads =
                new BufferedReader( new FileReader( file ) );
            String line = null;

            sequenceLength = Integer.parseInt( reads.readLine() );
            dots = new float[sequenceLength][];
            for( int i = sequenceLength; i > 0; i-- ) {
                dots[sequenceLength - i] = new float[i];
                for( int j = 0; j < i; j++ ) {
                    dots[sequenceLength - i][j] = Float.POSITIVE_INFINITY;
                }
            }

            String label = reads.readLine().split( "\t" )[2];
            String dividerType = ( label.equals( "kcal/mol" ) ) ?
                "Free Energy (kcal/mol)" :"-log10(BP Probability)";
            setDivider( dividerType );
            setPlotType( ( label.equals( "kcal/mol" ) ) ? 0 : 1 );

            float min = Float.MAX_VALUE, max = Float.MAX_VALUE * -1;
            while( ( line = reads.readLine() ) != null ) {
                String[] data = line.split( "\t" );
                if( data.length == 3 ) {
                    int x = Integer.parseInt( data[0] );
                    int y = Integer.parseInt( data[1] );

                    float value = Float.parseFloat( data[2] );
                    if( value < min ) { min = value; }
                    else if( value > max ) { max = value; }

                    dots[x-1][y-x] = value;
                } else { showError(); }
            }

            reads.close();

            scale = (double)595 / (double)(getPreferredPanelSize().width);
            setLegendAndBounds( min, max );
        } catch( Exception e ) { showError(); }
    }

    /**
     * Show an error message and set the error flag.
     */
    private void showError() {
        RNAstructureInfoDialog.error( "Error preparing dot plot." );
        preparationError = true;
    }
}
