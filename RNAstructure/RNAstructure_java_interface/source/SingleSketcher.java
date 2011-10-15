/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.Dimension;

/**
 * A class that handles a sketcher which holds a specialized image panel and
 * legend, used with an Imager.
 * @author Jessica Reuter
 */
public abstract class SingleSketcher extends Sketcher {
    /**
     * The name of the file from which this drawing was created
     */
    protected String file;

    /**
     * The data holder for the sketched strand
     */
    protected StrandDataHolder holder;

    /**
     * The image panel this sketcher draws on
     */
    protected PrintablePanel imagePanel;

    /**
     * The legend for the sketched image, as an array.
     * In each row, indices are:
     *         0: text
     *         1: color
     */
    protected Object[][] legend;

    /**
     * The sequence of bases in the strand being sketched
     */
    protected char[] sequence;

    /**
     * The sequence length
     */
    protected int sequenceLength;

    /**
     * The strand type that needs to be made.
     */
    private int type;

    /**
     * Constructor
     * @param file  the file this sketcher is sketching
     */
    public SingleSketcher( String file ) {
        this.file = file;
        holder = new StrandDataHolder();
    }

    /**
     * Create a strand to be sketched, if necessary.
     */
    public void create() {
        if( file != null ) {
            if( file.equals( "" ) ) { file = selectFile(); }

            if( !file.equals( "" ) ) {
                holder.setStrand( new RNA( file, type, true ) );
            }
        }

        if( holder.getStrand() != null ) { inputSequence(); }
        else { preparationError = true; }
    }

    /**
     * Get the data holder this sketcher uses.
     * @return  the data holder
     */
    protected StrandDataHolder getDataHolder() { return holder; }

    /**
     * Get the file this sketcher draws
     * @return  the file
     */
    protected String getFile() { return file; }

    /**
     * Get the image panel this sketcher draws on.
     * @return  the image panel
     */
    protected PrintablePanel getImagePanel() { return imagePanel; }

    /**
     * Get the array which holds the legend.
     * @return  the legend array
     */
    protected Object[][] getLegend() { return legend; }

    /**
     * Get the preferred image panel size.
     * @return  the preferred image panel size
     */
    protected abstract Dimension getPreferredPanelSize();

    /**
     * Get the sequence for the strand being drawn.
     * @return  the sequence
     */
    protected char[] getSequenceArray() { return sequence; }

    /**
     * Get the size of the drawn sequence.
     * @return  the size of the sequence
     */
    protected int getSequenceSize() { return sequenceLength; }

    /**
     * Record the sequence from the RNA strand.
     */
    private void inputSequence() {
        RNA strand = holder.getStrand();
        sequenceLength = strand.GetSequenceLength();
        sequence = new char[sequenceLength];

        for( int i = 1; i <= sequenceLength; i++ ) {
            sequence[i-1] = strand.GetNucleotide( i );
            if( holder.isErrorStatus() ) {
                preparationError = true;
                return;
            }
        }
    }

    /**
     * Reset the file to its default, ""
     */
    protected void resetFile() { setFile( "" ); }

    /**
     * Set the file name.
     * @param file  the new file name
     */
    protected void setFile( String newName ) { file = newName; }

    /**
     * Set the strand type.
     * @param type  the type of strand to set
     */
    protected void setStrandType( int type ) { this.type = type; }
}