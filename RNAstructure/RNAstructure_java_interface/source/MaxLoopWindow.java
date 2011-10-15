/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that creates a small window allowing for the maximum allowed size of
 * internal or bulge loops to be changed.
 * @author Jessica Reuter
 */
public class MaxLoopWindow {
    /**
     * Constructor
     */
    public MaxLoopWindow() {
        StrandDataHolder dataGroup =
	    ((InputWindow)RNAstructure.getCurrentWindow()).getDataHolder();

        int size = (Integer)dataGroup.getData().get( 0 );
        String defaultString = Integer.toString( size );

        String message = "Maximum number of unpaired nucleotides<br/>in " +
            "internal/bulge loops:";
        String input = RNAstructureInfoDialog.input( message, defaultString );

        int finalInput = 0;
        if( input.trim().equals( "" ) ) { finalInput = size; }
        else { finalInput = Integer.parseInt( input ); }

        dataGroup.setData( 0, finalInput );
    }

    /**
     * An inner class which allows a MaxLoopWindow to be constructed from the
     * more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * Create a MaxLoopWindow.
         * @return  a new MaxLoopWindow
         */
        public MaxLoopWindow createWindow() {
            return new MaxLoopWindow();
        }
    }
}
