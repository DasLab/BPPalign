/* (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that handles allowance of intermolecular pairs in a hybrid strand.
 * @author Jessica Reuter
 */
public class PermitMolecularity extends Modifier {
    /**
     * Change permissibility of intermolecular pairs.
     */
    protected void modify() {
        HybridRNA hybrid = (HybridRNA)RNAstructure.getCurrentStrand();
        boolean forbid = hybrid.GetForbidIntramolecular();
        hybrid.SetForbidIntramolecular( !forbid );
    }
}