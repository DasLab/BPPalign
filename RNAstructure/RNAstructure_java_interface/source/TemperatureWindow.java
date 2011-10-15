/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that opens a window allowing the temperature to be changed.
 * @author Jessica Reuter
 */
public class TemperatureWindow {
    /**
     * Constructor
     */
    public TemperatureWindow() {
        double temperature = RNAstructure.getCurrentStrand().GetTemperature();
        String input = RNAstructureInfoDialog.input(
            "Temperature (degrees K):", Double.toString( temperature ) );

        if( !input.trim().equals( "" ) ) {
            try {
		double newTemp = Double.parseDouble( input );
                if( newTemp >= 0 ) {
		    RNAstructure.getCurrentStrand().SetTemperature( newTemp );
		} else { throw new IllegalArgumentException(); }
	    } catch( NumberFormatException ne ) {
		RNAstructureInfoDialog.error( "Temperature is not a number." );
	    } catch( Exception e ) {
		RNAstructureInfoDialog.error( "Temperature must be >= 0." );
	    }
        }
    }

    /**
     * An inner class which allows a TemperatureWindow to be constructed from
     * the more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * Create a TemperatureWindow.
         * @return  a new TemperatureWindow
         */
        public TemperatureWindow createWindow() {
            return new TemperatureWindow();
        }
    }
}