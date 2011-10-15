/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.io.BufferedWriter;
import java.io.FileWriter;

/**
 * A class that writes dot plot data to a file. This class can write dot plot
 * data files or Postscript files.
 * @author Jessica Reuter
 */
public class DotPlotFileWriter extends ActionBase {
    /**
     * The type of dot plot file being written
     */
    private String type;

    /**
     * Constructor
     * @param type  the type of file to write
     */
    private DotPlotFileWriter( String type ) { this.type = type; }

    /**
     * Write a dot plot to a file, if a file was properly selected.
     */
    public void act() {
        RNAstructure.setLabel( "For help, press F1." );
        FileSelector selector = new FileSelector( false, "save", type );
        String file = selector.showSelector();

        if( !file.equals( "" ) ) {
            if( file.endsWith( "dp" ) ) { writeDotPlotFile( file ); }
            else { writePostscriptFile( file ); }
        }
    }

    /**
     * Export a dot plot to a dot plot file. A dot plot file shows only those
     * dots that are currently visible.
     * @return  the DotPlotFileWriter that writes the dot plot file
     */
    public static DotPlotFileWriter exportToDotPlotFile() {
        return new DotPlotFileWriter( "Dot Plot Files,dp" );
    }

    /**
     * Export a dot plot to a Postscript file. A Postscript file generates an
     * image of the entire dot plot, whether dots are currently visible or not.
     * @return  the DotPlotFileWriter that writes the Postscript file
     */
    public static DotPlotFileWriter exportToPostscriptFile() {
        return new DotPlotFileWriter( "Postscript Files,ps" );
    }

    /**
     * Write the visible dots to a dot plot data file.
     * @param fileName  the name of the file to write to
     */
    private void writeDotPlotFile( String fileName ) {
        try {
	    /*
            BufferedWriter writer =
                new BufferedWriter( new FileWriter( fileName ) );
	    */
            Imager imager = (Imager)RNAstructure.getCurrentWindow();
            DrawDotPlot plotter = (DrawDotPlot)imager.getSketcher();
	    String title = imager.getFile();
            //int size = RNAstructure.getCurrentStrand().GetSequenceLength();

	    DotPlotHandler handler = new DotPlotHandler( title, fileName );

	    if( title.endsWith( ".sav" ) ) {
		handler.readFoldingData();
            } else if( title.endsWith( ".pfs" ) ) {
		handler.readPartitionData();
            } else if( title.endsWith( ".dsv" ) ) {
                int strand =
                    ((DrawDotPlotFromDynalignSave)plotter).getStrandNumber();

                if( strand == 1 ) {
		    handler.readDynalignSeq1Data();
                } else if( strand == 2 ) {
		    handler.readDynalignSeq2Data();
                }
            }

	    handler.setMinimum( plotter.getVisibleRange()[0] );
	    handler.setMaximum( plotter.getVisibleRange()[1] );
	    handler.writePlotFile();

	    /*
            int type = plotter.getPlotType();
            String label = ( type == 0 ) ? "kcal/mol" : "Probability";

            writer.write( Integer.toString( size ) );
            writer.newLine();
            writer.write( "i\tj\t" + label );
            writer.newLine();

            for( int i = 1; i < size; i++ ) {
                for( int j = i; j < size; j++ ) {
                    float value = plotter.getDotValue( j, i );

                    if( plotter.getDotIsVisible( value ) ) {
                        writer.write( i + "\t" + j + "\t" + value );
                        writer.newLine();
                        writer.flush();
                    }
                }
            }

            writer.close();
	    */

	    String message = ( !handler.isError() ) ?
		"Dot Plot File Written." : "";
            RNAstructureInfoDialog.message( "Error Writing Dot Plot File." );
        } catch( Exception e ) {
            RNAstructureInfoDialog.error( "Error Writing Dot Plot File." );
        }
    }

    /**
     * Write all dots to a Postscript image file.
     * @param fileName  the name of the file to write to
     */
    private void writePostscriptFile( String fileName ) {
        try {
            DrawDotPlot plotter =
                (DrawDotPlot)(((Imager)RNAstructure.getCurrentWindow())
                .getSketcher());
            //int entries = plotter.getEntries();

            String input =
		((Imager)RNAstructure.getCurrentWindow()).getTitle();
            //boolean energyPlot = input.endsWith( ".sav" );
	    float[] range = plotter.getVisibleRange();

	    DotPlotHandler plotHandler = new DotPlotHandler( input, fileName );
	    plotHandler.setEntries( plotter.getEntries() );
	    plotHandler.setMinimum( range[0] );
	    plotHandler.setMaximum( range[1] );

            Postscript_Wrapper post = new Postscript_Wrapper();
	    boolean ok = false;
            if( input.endsWith( ".sav" ) ) {
		ok = post.plotEnergy( plotHandler );
	    } else if( input.endsWith( ".pfs" ) ) {
		ok = post.plotProbability( plotHandler );
	    } else if( input.contains( ".dsv" ) ) {
		int strand =
		    ((DrawDotPlotFromDynalignSave)plotter).getStrandNumber();

		if( strand == 1 ) {
		    ok = post.plotDynalign1( plotHandler );
		} else if( strand == 2 ) {
		    ok = post.plotDynalign2( plotHandler );
		}
	    }

	    String message = ( ok ) ?
		"Postscript Plot Written." : "Error Writing Postscript Plot.";
            RNAstructureInfoDialog.message( message );
        } catch( Exception e ) {
            RNAstructureInfoDialog.error( "Error Writing Postscript File." );
        }
    }
}
