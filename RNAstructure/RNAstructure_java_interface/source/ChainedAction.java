/* 
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that links together ActionBases.
 * @author Jessica Reuter
 */
public class ChainedAction extends ActionBase {
    /**
     * The array of ActionBases to link
     */
    private ActionBase[] chain;

    /**
     * Constructor
     * @param chain  the chain of action listeners
     */
    public ChainedAction( ActionBase... actions ) { chain = actions; }

    /**
     * Do a chain of ActionBase actions
     */
    protected void act() {
        for( ActionBase element: chain ) { element.act(); }
    }

    /**
     * Get the chain of actions.
     * @return  the chain
     */
    protected ActionBase[] getChain() { return chain; }
}