/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JButton;

/**
 * A class that holds all RNAstructure windows into which input is given.
 * @author Jessica Reuter
 */
public abstract class InputWindow extends MenuChangingWindow {
    /**
     * The data holder which holds all input into this window.
     */
    protected StrandDataHolder holder;

    /**
     * Boolean, true if progress should be monitored, false if not
     * Default is true.
     */
    protected boolean progress = true;

    /**
     * Boolean, true if progress should be monitored in indeterminate mode,
     * false if not
     * Default is false.
     */
    protected boolean progressIndeterminate = false;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = 5512655570078481655L;

    /**
     * Add a piece of data to the input for this window.
     * @param object  the new piece of input
     */
    protected void addData( Object object ) {
	holder.getData().add( object );
    }

    /**
     * Create a start button, and initialize the data holder.
     * @param action  the action when the start button is clicked
     * @return  the start button
     */
    protected JButton createStartButton( ActionBase action ) {
	holder = new StrandDataHolder();

        JButton button = new JButton( "START" );
        button.addActionListener( action );
        return button;
    }

    /**
     * Get the data holder for this window.
     * @return  the data holder
     */
    protected StrandDataHolder getDataHolder() { return holder; }

    /**
     * Get the number of input points for this window.
     * @return  the number of input points
     */
    protected int getInputDataNumber() {
        return holder.getData().size();
    }

    /**
     * Check if all necessary data has been entered; if the window is ready to
     * begin calculations
     */
    protected abstract boolean isReady();

    /**
     * Check if progress should be monitored on this window.
     * @return  true if progress should be monitored, false if not
     */
    protected boolean monitorProgress() { return progress; }

    /**
     * Check if progress should be monitored in indeterminate mode.
     * @return  true if progress should be monitored, false if not
     */
    protected boolean monitorProgressIndeterminate() {
	return progressIndeterminate;
    }

    /**
     * Read in SHAPE constraints for a strand's calculation.
     * @param index  the index at which SHAPE constraints start
     * @param length  the length of the data array
     * @return  true if addition of constraints was successful, false if not
     */
    protected boolean readSHAPEConstraints( int index, int length ) {
        for( int i = index; i < length; i++ ) {
            Object[] shape = (Object[])holder.getData().get( i );
            String file = shape[0].toString();
            double p1 = (Double)shape[1];
            double p2 = (Double)shape[2];

            boolean energy = ( shape.length == 3 );
            int code = holder.getStrand().ReadSHAPE( file, p1, p2, energy );

            if( holder.isErrorStatus( code ) ) { return false; }
        }
        return true;
    }

    /**
     * Remove a piece of data, by index, from this window.
     * @param index  the index of the data to be removed
     */
    protected void removeData( int index ) {
	holder.getData().remove( index );
    }

    /**
     * Add the data into the data list for the strand this window is creating.
     * @throws Exception  if setting of data ran into an error
     */
    protected abstract void setData() throws Exception;
}
