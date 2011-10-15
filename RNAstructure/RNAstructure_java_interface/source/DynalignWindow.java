/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A base class for Dynalign input windows.
 * @author Jessica Reuter
 */
public abstract class DynalignWindow extends InputWindow {
    /**
     * The current strand this window is working with, ie. for constraints
     */
    protected RNA currentStrand;

    /**
     * The Dynalign object that holds the strands for alignment
     */
    protected Dynalign_object dynalign;

    /**
     * The main input panel with file names
     */
    protected TextInputPanel inputPanel;

    /**
     * The input panel holding option data
     */
    protected TextInputPanel inputPanel2;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -2047610144631680095L;

    /**
     * Constructor
     * @param title  the title of the window
     */
    protected DynalignWindow( String title ) { setTitles( title ); }

    /**
     * Create the main input panel.
     * @param amount  the number of buttons in the panel
     * @param names  the list of button labels
     */
    protected void createFilePanel( int amount, String... names ) {
        inputPanel = new TextInputPanel();
        setBorder( inputPanel, 1, null );
        inputPanel.create( amount, 1, 30, names );
    }

    /**
     * Create the structure parameters panel.
     */
    protected void createStructurePanel() {
        String[] labels = {
            "Max % Energy Difference",
            "Max Number of Structures",
            "Structure Window Size",
            "Alignment Window Size"
        };

        inputPanel2 = new TextInputPanel();
        inputPanel2.create( 4, 2, 10, labels );
        setBorder( inputPanel2, 2, "Suboptimal Structure Parameters" );
    }

    /**
     * Get the Dynalign object.
     * @return  the Dynalign object.
     */
    protected Dynalign_object getDynalign() { return dynalign; }

    /**
     * Get the file panel.
     * @return  the file panel.
     */
    protected TextInputPanel getFilePanel() { return inputPanel; }

    /**
     * Get the parameter panel.
     * @return  the parameter panel
     */
    protected TextInputPanel getParameterPanel() { return inputPanel2; }

    /**
     * Get the strand that is currently being worked with.
     * @return  the strand
     */
    public RNA getWorkingStrand() { return currentStrand; }

    /**
     * Check if all required data fields have been filled in. This method is a
     * wrapper to call the TextInputPanel method checkValues; see
     * TextInputPanel.java for more information.
     * @return  true if all required data is filled, false if not
     */
    public boolean isReady() {
        return inputPanel.checkValues() && inputPanel2.checkValues();
    }

    /**
     * Get the necessary data from this Dynalign window.
     * This method doesn't do anything, since it isn't used in this class or
     * its subclasses. The Dynalign_object handles data for the strands, so the
     * StrandDataHolder isn't needed here.
     */
    protected void setData() { /* Do nothing. */ }

    /**
     * Set the dynalign object.
     * @param object  the Dynalign_object to set
     */
    protected void setDynalign( Dynalign_object object ) { dynalign = object; }

    /**
     * Set the working strand.
     * @param strand  the strand to work with
     */
    public void setWorkingStrand( RNA strand ) { currentStrand = strand; }
}