/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JButton;

/**
 * A class that handles setting SHAPE constraints on a strand.
 * @author Jessica Reuter
 */
public class SHAPEConstraintSetter extends ActionBase {
    /**
     * Whether these are energy constraints (true) or not (false)
     */
    private boolean energy;

    /**
     * Constructor
     */
    private SHAPEConstraintSetter( boolean energy ) { this.energy = energy; }

    /**
     * Set SHAPE constraints
     */
    protected void act() {
        SHAPEWindow current =
            (SHAPEWindow)(((JButton)event.getSource()).getTopLevelAncestor());

        String file = current.panel.getFields()[0].getText();
        Double par1 =
            Double.parseDouble( current.panel2.getFields()[0].getText() );
        Double par2 =
            Double.parseDouble( current.panel2.getFields()[1].getText() );

        RNAstructure.getCurrentStrand().ReadSHAPE( file, par1, par2, energy );
    }

    /**
     * Create a listener that sets energy constraints.
     * @return  the energy constraints action
     */
    public static SHAPEConstraintSetter setEnergy() {
        return new SHAPEConstraintSetter( true );
    }

    /**
     * Create a listener that sets hard constraints.
     * @return  the hard constraints action
     */
    public static SHAPEConstraintSetter setHard() {
        return new SHAPEConstraintSetter( false );
    }
}