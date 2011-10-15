/* (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that handles manual exiting of the application
 * @author Jessica Reuter
 */
public class Exiter extends ActionBase {
    /**
     * Exit the application
     */
    protected void act() { System.exit( 0 ); }
}