/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A window which extends PartitionWindow to partition two strands of
 * nucleic acids
 * @author Jessica Reuter
 */
public class PartitionDoubleWindow extends PartitionWindow {
    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -7011701877595773683L;

    /**
     * Constructor
     * @param acid  the nucleic acid type being partitioned
     */
    public PartitionDoubleWindow( String acid ) {
        super( acid + " Bimolecular Partition Function" );
        boolean isRNA = acid.equals( "RNA" );

        // Set the menu bar
        setSections( 1, 3, 4, 5, 7 );
        setMenus(
            RNAstructure.getBar().getMenuMaker().createTemperatureMenu() );
        RNAstructure.getBar().refresh( sections, menus );

        // Create the input panel
        createFilePanel( 3, "Sequence File 1", "Sequence File 2",
                         "Save File" );
        inputPanel.disableField( 1 );
        inputPanel.disableField( 2 );

        // Create and set actions for buttons
        // First: Select the first sequence file and set its name.
        // Second: Select the second sequence file and set its name.
        // Third: Create a hybrid strand from the two sequence files.
        // Fourth: Set the name of the default output file.
        // Fifth: Select an output file other than default and set its name.
        // "First" executes alone on button 0.
        // "Second", "Third", and "Fourth" execute in sequence on button 1.
        // "Fifth" executes alone on button 2.
        OpenSaveSetter first =
            OpenSaveSetter.open( 0, "Sequence Files", "seq" );
        OpenSaveSetter second =
            OpenSaveSetter.open( 1, "Sequence Files", "seq" );
        StrandCreator third = StrandCreator.createHybridSEQStrand( 1, isRNA );
        OutputFiller fourth = OutputFiller.fillMultiple( 2, "pfs" );
        OpenSaveSetter fifth =
            OpenSaveSetter.save( 2, "Partition Function Save Files", "pfs" );
        ChainedAction partitionSequence =
            new ChainedAction( second, third, fourth );

        inputPanel.setAction( 0, first );
        inputPanel.setAction( 1, partitionSequence );
        inputPanel.setAction( 2, fifth );

        // Set components in place and designate start button action
        String type = "Bimolecular Partitioning";
        setComponents( new PartitionNucleicAcids( type ) );
    }

    /**
     * Set the data for this window.
     */
    public void setData() {
        removeData( 0 );
        addData( inputPanel.getFields()[1].getText() );
        addData( inputPanel.getFields()[2].getText() );
    }

    /**
     * An inner class which allows a PartitionDoubleWindow to be constructed
     * from the more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * The nucleic acid the created window handles
         */
        private String acid;

        /**
         * Constructor
         * @param acid  the nucleic acid the window handles
         */
        public Factory( String acid ) { this.acid = acid; }

        /**
         * Create a PartitionDoubleWindow.
         * @return  a new PartitionDoubleWindow
         */
        public PartitionDoubleWindow createWindow() {
            return new PartitionDoubleWindow( acid );
        }
    }
}