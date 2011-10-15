/* 
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that resets a strand by removing all folding constraints.
 * @author Jessica Reuter
 */
public class StrandResetter extends Modifier {
    /**
     * Reset all folding constraints on the strand.
     */
    protected void modify() {
        RNAstructure.getCurrentStrand().RemoveConstraints();
    }
}