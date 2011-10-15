/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

public class ForceDistance extends ActionBase {
    /**
     * The panel from which new distance data comes
     */
    private TextInputPanel panel;

    /**
     * Constructor
     */
    public ForceDistance( TextInputPanel panel ) { this.panel = panel; }

    /**
     * Force a maximum pairing distance constraint.
     */
    protected void act() {
        String value = panel.getFields()[0].getText().trim();
        if( !value.equals( "" ) ) {
            try {
                int distance = Integer.parseInt( value );

                if( distance < 0 || distance > 600 ) {
                    RNAstructureInfoDialog.error( "Invalid maximum pairing " +
                        "distance given.");
                } else {
                    RNAstructure.getCurrentStrand()
                        .ForceMaximumPairingDistance( distance );
                }
            } catch( Exception ex ) {
                RNAstructureInfoDialog.error( "Error setting " +
                    "maximum pairing distance." );
            }
        }
    }
}