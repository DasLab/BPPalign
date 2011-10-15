/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.io.BufferedWriter;
import java.io.FileWriter;

import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * A window that deals with either reading or writing Dynalign alignment
 * constraints so the RNAstructure GUI can either input or output them.
 * @author Jessica Reuter
 */
public class DynalignAlignmentIOHandler {
    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -335671844754176548L;

    /**
     * Constructor
     * @param type  the type of constraints action to be taken
     *              1.  Read constraints from a file
     *              2.  Write constraints to a file
     *              3.  Remove all constraints
     *              4.  Show all constraints
     */
    public DynalignAlignmentIOHandler( int type ) {
        RNAstructure.setLabel( "For help, press F1." );

        if( type == 1 ) { runIO( 1 ); }
        else if( type == 2 ) { runIO( 2 ); }
        else if( type == 3 ) { resetConstraints(); }
        else if( type == 4 ) { showConstraints(); }
        else throw new IllegalArgumentException( "Undefined window type" );
    }

    /**
     * Reset all alignment constraints.
     */
    private void resetConstraints() {
	ActionBase clearing = new ActionBase() {
		public void act() {
		    ((DynalignFoldWindow)RNAstructure.getCurrentWindow())
			.clearForcedAlignments();
		}
	    };

        RNAstructureInfoDialog.confirm(
            "This will erase all alignment constraints.\nContinue?",
            new ChainedAction( clearing, new CloseAction() ) );
    }

    /**
     * Either read or write alignment constraints from or to a file.
     * @param type  1=read, 2=write
     */
    public void runIO( int type ) {
        String title = "Open";
        if( type == 2 ) { title = "Save"; }

        String file = new FileSelector( true, title, "Constraint Files,con" )
            .showSelector();

        if( !file.equals( "" ) ) {
            Dynalign_object dynalign =
		((DynalignWindow)RNAstructure.getCurrentWindow())
		.getDynalign();
            if( type == 1 ) { dynalign.ReadAlignmentConstraints( file ); }
            else { writeDynalignAlignmentConstraintsFile( file ); }
        }
    }

    /**
     * Show all alignment constraints
     */
    private void showConstraints() {
        // Set up variables and data structures
        StringBuilder constraints =
	    new StringBuilder( "<html>Alignment Constraints<br/><br/>" );

	Integer[][] alignments =
	    ((DynalignFoldWindow)RNAstructure.getCurrentWindow())
	    .getForcedAlignments();
	
	for( Integer[] entry: alignments ) {
	    constraints.append(
		"(1)" + Integer.toString( entry[0] ) + " -- " +
		"(2)" + Integer.toString( entry[1] ) + "<br/>" );
	}

        if( alignments.length == 0 ) { constraints.append( "None" ); }

        // Show final list dialog
        RNAstructureInfoDialog.message( constraints.toString() );
    }

    /**
     * Write a Dynalign alignment constraints file.
     * @param file  the file name
     */
    private void writeDynalignAlignmentConstraintsFile( String file ) {
	try {
	    BufferedWriter writer =
		new BufferedWriter( new FileWriter( file ) );

	    Integer[][] alignments =
		((DynalignFoldWindow)RNAstructure.getCurrentWindow())
		.getForcedAlignments();

	    for( Integer[] entry: alignments ) {
		String one = Integer.toString( entry[0] );
		String two = Integer.toString( entry[1] );
		writer.write( one + "\t" + two );
		writer.newLine();
		writer.flush();
	    }

	    writer.write( "-1\t-1" );
	    writer.newLine();
	    writer.close();
	} catch( Exception e ) {
	    RNAstructureInfoDialog.error(
		"Error writing Dynalign alignment file." );
	}
    }

    /**
     * An inner class which allows a ConstraintsIOHandler to be constructed
     * from the more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * The type of DynalignAlignmentIOHandler to create.
         */
        private int type;

        /**
         * Constructor
         * @param type  the type of window to create
         */
        public Factory( int type ) { this.type = type; }

        /**
         * Create a DynalignAlignmentIOHandler.
         * @return  a new DynalignAlignmentIOHandler
         */
        public DynalignAlignmentIOHandler createWindow() {
            return new DynalignAlignmentIOHandler( type );
        }
    }
}
