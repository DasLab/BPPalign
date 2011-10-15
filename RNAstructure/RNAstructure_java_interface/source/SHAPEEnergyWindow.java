/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A subclass of SHAPEWindow that handles the addition of pseudoenergy SHAPE
 * constraints to a strand of nucleic acids
 * @author Jessica Reuter
 */
public class SHAPEEnergyWindow extends SHAPEWindow {
    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -6269240604457152888L;

    /**
     * Constructor
     */
    public SHAPEEnergyWindow() {
        // Create the specific pseudoenergy constraints input panel
        String[] energy = { "Slope (kcal/mol):", "Intercept (kcal/mol):" };
        panel2.create( 2, 2, 5, energy );

        // Set the initial contents of the panel's text fields
        panel2.getFields()[0].setText( "2.6" );
        panel2.getFields()[1].setText( "-0.8" );

        // Set the specific action for the "OK" button
        ChainedAction energyAction =
            new ChainedAction(
                SHAPEConstraintSetter.setEnergy(), new CloseAction() );
        buttonPanel.setAction( energyAction, 0, "OK" );

        // Set components in their places.
        setComponents();
    }

    /**
     * An inner class which allows a SHAPEEnergyWindow to be constructed from
     * the more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * Create a SHAPEEnergyWindow.
         * @return  a new OligoScreenWindow
         */
        public SHAPEEnergyWindow createWindow() {
            return new SHAPEEnergyWindow();
        }
    }
}