/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.util.ArrayList;

/**
 * A class that holds a strand being analyzed, as well as associated data.
 * @author Jessica Reuter
 */
public class StrandDataHolder {
    /**
     * Data about the strand needed for calculation later on.
     * Index 0 is always maximum loop size.
     */
    private ArrayList<Object> data;

    /**
     * The strand of nucleic acids
     */
    private RNA strand;

    /**
     * Constructor
     */
    public StrandDataHolder() {
        data = new ArrayList<Object>();
        data.add( 30 );
    }

    /**
     * Add data to the data list for this strand.
     * @param points  the new data points
     */
    public void addData( Object... points ) {
        for( Object element: points ) { data.add( element ); }
    }

    /**
     * Get the data for this strand
     * @return  the strand data
     */
    public ArrayList<Object> getData() { return data; }

    /**
     * Get the data point at a particular index in this holder.
     * @return  the data point
     */
    public Object getDataAt( int index ) { return data.get( index ); }

    /**
     * Get the strand
     * @return  the strand
     */
    public RNA getStrand() { return strand; }

    /**
     * Report an error in this strand.
     * @return  true if an error occurred, false if not
     */
    public boolean isErrorStatus() {
        return isErrorStatus( strand.GetErrorCode() );
    }

    /**
     * Report an error in this strand.
     * @param code  the error code
     * @return  true if an error occurred, false if not
     */
    public boolean isErrorStatus( int code ) {
        if( code != 0 ) {
            RNAstructureInfoDialog.error( strand.GetErrorMessage( code ) );
            RNAstructure.getCurrentWindow().dispose();
            return true;
        }
        return false;
    }

    /**
     * Remove a specified piece of data from the holder.
     * @param index  the index in the list of the data to remove
     */
    public void removeData( int index ) { data.remove( index ); }

    /**
     * Set the data point at a particular index in the data list.
     * @param index  the index in the list
     * @param pt  the new data point
     */
    public void setData( int index, Object pt ) { data.set( index, pt ); }

    /**
     * Set the strand
     * @param strand  the strand to set
     */
    public void setStrand( RNA strand ) { this.strand = strand; }
}
