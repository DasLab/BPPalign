/* 
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that sets annotation on a drawn structure.
 * @author Jessica Reuter
 */
public class AnnotationWriter extends ActionBase {
    /**
     * The type of annotation (or lack thereof)
     */
    private int type;

    /**
     * Constructor
     * @param type  the type of writer to create
     */
    private AnnotationWriter( int type ) { this.type = type; }

    /**
     * Create a writer that places probability annotation
     * @return  the writer
     */
    public static AnnotationWriter probability() {
        return new AnnotationWriter( 0 );
    }

    /**
     * Create a writer that removes all annotation.
     * @return  the writer
     */
    public static AnnotationWriter remove() {
        return new AnnotationWriter( 2 );
    }

    /**
     * Create a writer that places SHAPE annotation.
     * @return  the writer
     */
    public static AnnotationWriter shape() {
        return new AnnotationWriter( 1 );
    }

    /**
     * Annotate, or remove annotation from, a structure
     */
    protected void act() {
        RNAstructure.setLabel( "For help, press F1." );

        DrawStructure sketcher = 
            (DrawStructure)((Imager)RNAstructure.getCurrentWindow())
	    .getSketcher();

        // Create probability annotation
        if( type == 0 ) {
            int structure = sketcher.getCurrentStructure();

            FileSelector selector1 = new FileSelector( false, "open",
                "Partition Function Save Files,pfs" );
            String file = selector1.showSelector();

            if( !file.equals( "" ) ) {
                sketcher.setProbabilityArray(
                    AnnotationMapHolder.readProbabilityColors(
                        file, structure ) );
		sketcher.setProbabilityAnnotationFile( file );
                sketcher.setSHAPEArray( null );
                sketcher.setProbabilityAnnotated();
            }
        }

        // Create SHAPE annotation
        else if( type == 1 ) {
            int bases = sketcher.getSequenceSize();

            FileSelector selector2 = new FileSelector( false, "open",
                "SHAPE Files,shape" );
            String file = selector2.showSelector();

            if( !file.equals( "" ) ) {
                sketcher.setSHAPEArray(
                    AnnotationMapHolder.readSHAPEColors( file, bases ) );
		sketcher.setSHAPEAnnotationFile( file );
		sketcher.setProbabilityArray( null );
                sketcher.setSHAPEAnnotated();
            }
        }

        // Remove all annotation
        else {
            sketcher.setProbabilityArray( null );
            sketcher.setSHAPEArray( null );
            sketcher.setUnannotated();
        }

	// Set the focus to the imager.
	((Imager)RNAstructure.getCurrentWindow()).requestFocus();
    }
}
