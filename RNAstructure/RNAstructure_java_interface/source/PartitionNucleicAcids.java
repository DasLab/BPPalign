/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.util.ArrayList;

/**
 * A class that handles calculation of the partition function for one or two
 * strands of nucleic acids.
 * @author Jessica Reuter
 */
public class PartitionNucleicAcids extends StrandStarter {
    /**
     * The type of partitioning being done, either single or double stranded
     */
    private String type;

    /**
     * Constructor
     * @param type  the type of partitioning being done
     */
    public PartitionNucleicAcids( String type ) { this.type = type; }

    /**
     * Calculate the partition function.
     */
    public void execute() {
        try {
	    StrandDataHolder holder =
		((InputWindow)RNAstructure.getCurrentWindow())
		.getDataHolder();
            ArrayList<Object> data = holder.getData();
            RNA strand = holder.getStrand();

            // Partition the strand(s)
            final String outputFile = data.get( 0 ).toString();
            int code = 0;

            if( type.equals( "Single Strand Partitioning" ) ) {
                code = strand.PartitionFunction( outputFile );
            } else {
                code = ((HybridRNA)strand)
                    .PartitionFunctionBimolecular( outputFile );
            }

            // Draw the dot plot
            if( !holder.isErrorStatus( code ) ) {
                DrawDotPlot plotter =
                    DrawDotPlot.drawProbabilityPlot( outputFile, 5 );

                new ChainedAction(
                    new CloseAction(),
                    plotter,
                    new WindowDisplayer( new Imager.DotPlotFactory(
                        plotter ) ) ).act();
            }
        } catch( Exception ex ) {
            RNAstructureInfoDialog.message( "Undefined Error During " +
                type + "." );
        }
    }
}
