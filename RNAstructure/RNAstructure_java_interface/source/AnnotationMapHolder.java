/*
 *  (c) 2009  Mathews Lab, University of Rochester Medical Center
 *  
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Vector;

/**
 * A class which holds maps that hold annotation color key data.
 * @author Jessica Reuter
 */
public final class AnnotationMapHolder {
    /**
     * The array of colors used in annotation maps
     * Row 0: Probability Colors
     * Row 1: SHAPE colors
     */
    private static Color[][] colors = {
        { Color.red, new Color( 249, 169, 57 ), Color.pink,
          new Color( 0, 100, 0 ), new Color( 127, 255, 0 ),
          new Color( 127, 255, 212 ), Color.blue, Color.magenta },
        { new Color( 169, 169, 169 ), Color.red, new Color( 249, 169, 57 ),
          Color.black }
    };

    /**
     * The vector of annotation maps
     */
    private static Vector<Map<String, Color>> maps = createMaps();

    /**
     * Create the maps
     */
    private static Vector<Map<String, Color>> createMaps() {
        Vector<Map<String, Color>> maps = new Vector<Map<String, Color>>( 2 );

        LinkedHashMap<String, Color> probs =
            new LinkedHashMap<String, Color>();

        probs.put( "BP Probability >= 99%", colors[0][0] );
        probs.put( "99% > BP Probability >= 95%", colors[0][1] );
        probs.put( "95% > BP Probability >= 90%", colors[0][2] );
        probs.put( "90% > BP Probability >= 80%", colors[0][3] );
        probs.put( "80% > BP Probability >= 70%", colors[0][4] );
        probs.put( "70% > BP Probability >= 60%", colors[0][5] );
        probs.put( "60% > BP Probability >= 50%", colors[0][6] );
        probs.put( "50% > BP Probability", colors[0][7] );
        maps.add( Collections.unmodifiableMap( probs ) );

        LinkedHashMap<String, Color> shape =
            new LinkedHashMap<String, Color>();

        shape.put( "No Data", colors[1][0] );
        shape.put( "SHAPE >= 0.7", colors[1][1] );
        shape.put( "0.7 > SHAPE >= 0.3", colors[1][2] );
        shape.put( "0.3 > SHAPE", colors[1][3] );
        maps.add( Collections.unmodifiableMap( shape ) );

        return maps;
    }

    /**
     * Get the probability map.
     * @return  the probability map
     */
    public static Map<String, Color> getChanceMap() { return maps.get( 0 ); }

    /**
     * Read and return the probability annotation values for a particular file
     * and structure
     * @param file  the file holding the annotation values to read
     * @param structure  the structure to read values for
     * @return  the array of color values
     */
    public static Color[] readProbabilityColors( String file, int structure ) {
        RNA strand = null;
        int code = 0;
        Color[] annotations = null;

        try {
            strand = new RNA( file, 3, true );
            code = strand.GetErrorCode();

            if( code == 0 ) {
                int length = strand.GetSequenceLength();
                annotations = new Color[length];
                RNA strandOrig = RNAstructure.getCurrentStrand();

                for( int i = 1; i <= length; i++ ) {
                    int pair = strandOrig.GetPair( i, structure );
                    code = strandOrig.GetErrorCode();

                    if( code == 0 ) {
			if( pair == 0 ) { annotations[i-1] = Color.black; }

                        if( ( pair != 0 ) && ( pair > i ) ) {
                            double value =
                                strand.GetPairProbability( i, pair );
                            code = strand.GetErrorCode();

                            if( code == 0 ) {
                                annotations[i-1] =
                                    ( value >= 0.99 ) ? colors[0][0] :
                                    ( value >= 0.95 ) ? colors[0][1] :
                                    ( value >= 0.90 ) ? colors[0][2] :
                                    ( value >= 0.80 ) ? colors[0][3] :
                                    ( value >= 0.70 ) ? colors[0][4] :
                                    ( value >= 0.60 ) ? colors[0][5] :
                                    ( value >= 0.50 ) ? colors[0][6] :
                                    colors[0][7];
				annotations[pair-1] = annotations[i-1];
                            } else { throw new Exception(); }
                        }
                    } else { throw new Exception(); }
                }
            } else { throw new Exception();}
        } catch( Exception e ) {
            RNAstructureInfoDialog.error( strand.GetErrorMessage( code ) );
        }

        return annotations;
    }

    /**
     * Read and return the SHAPE annotation values for a particular file.
     * @param file  the file holding the annotation values to read
     * @param bases  the number of bases that should be read
     * @return  the array of color values
     */
    public static Color[] readSHAPEColors( String file, int bases ) {
        Color[] annotations = new Color[bases];

        try {
            String line = null;
            BufferedReader reader =
                new BufferedReader( new FileReader( file ) );
            int counter = 0;

            while( ( line = reader.readLine() ) != null ) {
                double value = Double.parseDouble( line.split( "\t" )[1] );

                annotations[counter++] =
                    ( value == -999 ) ? colors[1][0] :
                    ( value >= 0.7 ) ? colors[1][1] :
                    ( value >= 0.3 ) ? colors[1][2] :
                    colors[1][3];
            }
            reader.close();
        } catch( Exception e ) {
            RNAstructureInfoDialog.error( "Error reading SHAPE data file." );
        }

        return annotations;
    }

    /**
     * Get the SHAPE map.
     * @return  the probability map
     */
    public static Map<String, Color> getSHAPEMap() { return maps.get( 1 ); }
}