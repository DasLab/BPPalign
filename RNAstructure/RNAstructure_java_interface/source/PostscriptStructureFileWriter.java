/*
 * (c) 2010  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.io.BufferedWriter;
import java.io.FileWriter;

import java.util.ArrayList;

/**
 * A class that writes structure files in Postscript.
 * @author Jessica Reuter
 */
public class PostscriptStructureFileWriter extends StructureFileWriter {
    /**
     * Boolean, true if circular structures are exported, false if not.
     */
    private boolean circular;

    /**
     * Boolean, true if plain structures are exported, false if not.
     */
    private boolean plain;

    /**
     * Boolean, true if probability annotated structures are exported, false if
     * not.
     */
    private boolean probability;

    /**
     * Boolean, true if radial structures are exported, false if not.
     */
    private boolean radial;

    /**
     * Boolean, true if SHAPE annotated structures are exported, false if not.
     */
    private boolean shape;

    /**
     * Constructor
     * @param types  the list of selections denoting which file types to write
     */
    private PostscriptStructureFileWriter( ArrayList<String> types ) {
        type = "Postscript Files,ps";

	circular = types.contains( "Circular" );
	radial = types.contains( "Radial" );

	plain = types.contains( "Plain" );
	probability = types.contains( "Probability" );
	shape = types.contains( "SHAPE" );
    }

    /**
     * Select an annotation file.
     * @param annoType  the annotation file type
     * @param filter  the filter for the FileSelector
     * @return  the chosen file name.
     */
    private String annotate( String annoType, String filter ) {
	String file = open( annoType, filter );

	if( file.equals( "" ) ) {
	    RNAstructureInfoDialog.error(
                "<html><center>" +
                annoType + " file not selected.<br/>" +
                "An annotated file will not be written."
	    );
	}

	return file;
    }

    /**
     * Choose a file with a FileSelector.
     * @param act the type of action the file chooser does.
     * @param fileType  the type of file to choose.
     * @param filter  the filter type for the FileSelector.
     * @return the chosen file name.
     */
    private String chooseFile( String act, String fileType, String filter ) {
	String caption = act + " " + fileType + " File";
	FileSelector selector = new FileSelector( false, caption, type );
	return selector.showSelector();
    }

    /**
     * Prepare to export the structure to one or more postscript files.
     * @param types  the array of selections denoting which file types to write
     * @return  a new PostscriptStructureFileWriter
     */
    public static PostscriptStructureFileWriter exportToPostscript(
        ArrayList<String> types ) {

        return new PostscriptStructureFileWriter( types );
    }

    /**
     * Open a file with a FileSelector.
     * @param fileType  the type of file to choose.
     * @param filter  the filter type for the FileSelector.
     * @return  the chosen file name
     */
    private String open( String fileType, String filter ) {
	return chooseFile( "Open", fileType, filter );
    }

    /**
     * Open a probability annotation file with a FileSelector.
     * @return  the selected file name
     */
    private String probability() {
	return annotate(
	    "Probability Annotation",
            "Partition Function Save Files,pfs"
        );
    }

    /**
     * Save a file with a FileSelector.
     * @param fileType  the type of file to choose.
     * @param filter  the filter type for the FileSelector.
     * @return  the chosen file name
     */
    private String save( String fileType, String filter ) {
	return chooseFile( "Save", fileType, filter );
    }

    /**
     * Open a SHAPE annotation file with a FileSelector.
     * @return  the selected file name
     */
    private String shape() {
        return annotate(
	    "SHAPE Annotation",
            "SHAPE Data Files,shape"
        );
    }

    /**
     * Show dialog box describing the result of a structure writing operation.
     * @param success  true if the writing finished successfully, false if not.
     * @param structureType  the type of structure that was written.
     */
    private void showResult( boolean success, String structureType ) {
	String ok =
	    "<html><center>" + structureType + "<br/>Structure File Written.";
	String error =
	    "<html><center>" + "Error writing " + structureType + "<br/>" +
	    "Structure File.";

	if( success ) { RNAstructureInfoDialog.message( ok ); }
	else { RNAstructureInfoDialog.error( error ); }
    }

    /**
     * Write the structure to a postscript file.
     */
    protected void writeFile() {
        // Get the input file name and the structure drawing.
        String inFile = ((Imager)RNAstructure.getCurrentWindow()).getTitle();
	DrawStructure sketcher =
	    (DrawStructure)(((Imager)RNAstructure.getCurrentWindow())
	    .getSketcher());

	// Initialize the Postscript writing wrapper. If it isn't initialized
	// correctly, show an error and return.
	Postscript_Wrapper wrap = null;
	boolean valid =
	    ( circular || radial ) && ( plain || probability || shape );

	if( valid ) { wrap = new Postscript_Wrapper(); }
	else {
	    RNAstructureInfoDialog.error(
                "<html><center>" +
		"Main structure type<br/>" +
		"(radial or circular)<br/>" +
		"must be specified with at least one annotation type<br/>" +
		"(plain, probability, or SHAPE)" +
		"." );
	    return;
	}

	// Initialize the unique identifiers for each written structure type.
	String identifier = null;
	String outFile = null;
	String probFile = sketcher.getProbabilityAnnotationFile();
	String shapeFile = sketcher.getSHAPEAnnotationFile();

	// If circular plain structure was requested, draw it.
	if( circular && plain ) {
	    identifier = "Unannotated Circular Structure";
	    outFile = save( identifier, type );

	    if( !outFile.equals( "" ) ) {
		showResult(
		    wrap.structureCircular( inFile, outFile ),
		    identifier
		);
	    }
	}

	// If circular probability annotated structure was requested, draw it.
	if( circular && probability ) {
            identifier = "Probability Annotated Circular Structure";
	    outFile = save( identifier, type );
            if( probFile.equals( "" ) ) { probFile = probability(); }

            if( !outFile.equals( "" ) && !probFile.equals( "" ) ) {
                showResult(
		    wrap.structureCircular_Probability(
			inFile, outFile, probFile
                    ),
                    identifier
		);
            }
	}

	// If circular SHAPE annotated structure was requested, draw it.
	if( circular && shape ) {
            identifier = "SHAPE Annotated Circular Structure";
            outFile = save( identifier, type );
            if( shapeFile.equals( "" ) ) { shapeFile = shape(); }

            if( !outFile.equals( "" ) && !shapeFile.equals( "" ) ) {
                showResult(
		     wrap.structureCircular_SHAPE(
			 inFile, outFile, shapeFile
		    ),
                    identifier
	        );
            }
	}

	// If radial plain structure was requested, draw it.
	if( radial && plain ) {
            identifier = "Unannotated Radial";
	    outFile = open( identifier, type );

            if( !outFile.equals( "" ) ) {
                showResult(
		    wrap.structureRadial( inFile, outFile ),
                    identifier
		);
            }
	}

	// If radial probability annotated structure was requested, draw it.
	if( radial && probability ) {
            identifier = "Probability Annotated Radial Structure";
            outFile = save( identifier, type );
            if( probFile.equals( "" ) ) { probFile = probability(); }

            if( !outFile.equals( "" ) && !probFile.equals( "" ) ) {
                showResult(
                    wrap.structureRadial_Probability(
			inFile, outFile, probFile
		    ),
                    identifier
	        );
            }
	}

	// If radial SHAPE annotated structure was requested, draw it.
	if( radial && shape ) {
            identifier = "SHAPE Annotated Radial Structure";
            outFile = save( identifier, type );
            if( shapeFile.equals( "" ) ) { shapeFile = shape(); }

            if( !outFile.equals( "" ) && !shapeFile.equals( "" ) ) {
                showResult(
                    wrap.structureRadial_SHAPE(
			inFile, outFile, shapeFile
		    ),
                    identifier
	        );
            }
	}
    }
}
