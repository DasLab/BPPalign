/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.Dimension;
import java.awt.Window;

import javax.swing.JComponent;
import javax.swing.JSpinner;

/**
 * A class that handles zooming of an image on an Imager imagePanel.
 * @author Jessica Reuter
 */
public class ZoomHandler extends ActionBase {
    /**
     * A boolean telling if the request for zooming came from clicking a
     * SpinnerWindow (true) or from another source (false)
     */
    private final boolean fromSpinner;

    /**
     * The sketcher linked to this zoom handler
     */
    private SingleSketcher sketch;

    /**
     * Default constructor
     */
    public ZoomHandler() {
	fromSpinner = true;
	sketch = null;
    }

    /**
     * Constructor
     * @param fromSpin  true if the request for zooming came from clicking a
     *                  SpinnerWindow, false if it came from another source
     */
    public ZoomHandler( SingleSketcher sketch ) {
	fromSpinner = false;
	this.sketch = sketch;
    } 

    /**
     * Zoom a drawn image.
     */
    protected void act() {
        if( fromSpinner ) {
	    JSpinner spinner =
                ((SpinnerWindow)((JComponent)event.getSource())
		.getTopLevelAncestor()).getSpinner();
	    Imager imager =
                (Imager)(((Window)spinner.getTopLevelAncestor()).getOwner());
	    sketch = imager.getSketcher();

            String value = spinner.getValue().toString();
            sketch.setScale( Double.parseDouble( value ) / 100 );
        }

	PrintablePanel panel = sketch.getImagePanel();
	double scale = sketch.getScale();

        Dimension preferred = sketch.getPreferredPanelSize();
        Dimension scaled = new Dimension(
            (int)(preferred.width * scale), (int)(preferred.height * scale) );

        panel.setPreferredSize( scaled );
        panel.setMinimumSize( panel.getPreferredSize() );
        panel.setMaximumSize( panel.getPreferredSize() );
        panel.setSize( panel.getPreferredSize() );

        if( sketch instanceof DrawDotPlot ) {
            ((DrawDotPlot)sketch).refreshCutoffs();
        }
    }
}
