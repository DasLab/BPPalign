/* 
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A base class for all drawing actions.
 * @author Jessica Reuter
 */
public abstract class Sketcher extends ActionBase {
    /**
     * The file filters that allow for selection of proper file types
     */
    private String[] filters;

    /**
     * If preparation has encountered an error
     */
    protected boolean preparationError;

    /**
     * The scale at which the image is being viewed.
     * 1 is equal to full size, or 100% zoom.
     */
    protected double scale;

    /**
     * Constructor
     */
    protected Sketcher() { scale = 1; }

    /**
     * Begin sequence of events to create and display a drawing.
     */
    public void act() {
        RNAstructure.setLabel( "For help, press F1." );
        preparationError = false;
        create();
        if( !preparationError ) { prepare(); }
    }

    /**
     * Initialize data needed to draw an image.
     */
    public abstract void create();

    /**
     * Get the scale at which this sketcher is drawing.
     * @return  the scale
     */
    protected double getScale() { return scale; }

    /**
     * Prepare an image for drawing.
     */
    public abstract void prepare();

    /**
     * Select a nucleic acid strand to display.
     * @return  The new file name
     */
    public String selectFile() {
        FileSelector selector = new FileSelector( true, "Open", filters );
        return selector.showSelector();
    }

    /**
     * Set the filters for a specific sketcher type. The distinct filters in
     * the input string are separated by ";".
     * @param input  the input string
     */
    protected void setFilters( String input ) { filters = input.split( ";" ); }

    /**
     * Set the scale this sketcher is using.
     * @param scale  the new scale
     */
    protected void setScale( double scale ) { this.scale = scale; }
}