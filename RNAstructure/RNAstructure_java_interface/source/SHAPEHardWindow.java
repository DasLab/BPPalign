/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A subclass of SHAPEWindow that handles additions of hard SHAPE constraints
 * to a strand of nucleic acids
 * @author Jessica Reuter
 */
public class SHAPEHardWindow extends SHAPEWindow {
    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -5015772989155905161L;

    /**
     * Constructor
     */
    public SHAPEHardWindow() {
        // Create the specific hard constraints input panel
        String[] hard = {
            "Threshold for Force Single Stranded:",
            "Threshold for Chemical Modification:"
        };
        panel2.create( 2, 2, 5, hard );

        // Set the initial contents of the panel's text fields
        panel2.getFields()[0].setText( "2" );
        panel2.getFields()[1].setText( "1" );

        // Set the specific action for the "OK" button
        ChainedAction energyAction =
            new ChainedAction(
                SHAPEConstraintSetter.setHard(), new CloseAction() );
        buttonPanel.setAction( energyAction, 0, "OK" );

        // Set components in their places.
        setComponents();
    }

    /**
     * An inner class which allows a SHAPEHardWindow to be constructed from the
     * more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * Create a SHAPEHardWindow.
         * @return  a new SHAPEHardWindow
         */
        public SHAPEHardWindow createWindow() {
            return new SHAPEHardWindow();
        }
    }
}