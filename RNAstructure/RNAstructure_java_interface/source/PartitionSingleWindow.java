/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A window which extends PartitionWindow to partition a single strand of
 * nucleic acids
 * @author Jessica Reuter
 */
public class PartitionSingleWindow extends PartitionWindow {
    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -1537712297643684639L;

    /**
     * Constructor
     * @param acid  the nucleic acid type being partitioned
     */
    public PartitionSingleWindow( String acid ) {
        super( acid + " Single Strand Partition Function" );
        boolean isRNA = acid.equals( "RNA" );

        // Set the menu bar
        setSections( 1, 3, 4, 5, 7 );
        SmartMenuBar bar = RNAstructure.getBar();
        setMenus( bar.getMenuMaker().createTemperatureMenu(),
                  bar.getMenuMaker().createForceSingleMenu() );
        RNAstructure.getBar().refresh( sections, menus );

        // Create the input panel
        createFilePanel( 2, "Sequence File", "Save File" );
        inputPanel.disableField( 1 );

        // Create and set actions for buttons
        // First: Select the sequence file and set its name.
        // Second: Create a strand from the sequence file.
        // Third: Set the name of the default output file.
        // Fourth: Select an output file other than default and set its name.
        // "First", "Second", and "Third" execute in sequence on button 0.
        // "Fourth" executes alone on button 1.
        OpenSaveSetter first =
            OpenSaveSetter.open( 0, "Sequence Files", "seq" );
        StrandCreator second = StrandCreator.createSEQStrand( 0, isRNA );
        OutputFiller third = OutputFiller.fillSingle( 1, "pfs" );
        OpenSaveSetter fourth =
            OpenSaveSetter.save( 1, "Partition Function Save Files", "pfs" );
        ChainedAction partitionSequence =
            new ChainedAction( first, second, third );

        inputPanel.setAction( 0, partitionSequence );
        inputPanel.setAction( 1, fourth );

        // Set components in place and designate start button action
        String type = "Single Strand Partitioning";
        setComponents( new PartitionNucleicAcids( type ) );
    }

    /**
     * Set the data for this window.
     */
    public void setData() {
        removeData( 0 );
        addData( inputPanel.getFields()[1].getText() );
    }

    /**
     * An inner class which allows a PartitionSingleWindow to be constructed
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
         * Create a PartitionSingleWindow.
         * @return  a new PartitionSingleWindow
         */
        public PartitionSingleWindow createWindow() {
            return new PartitionSingleWindow( acid );
        }
    }
}