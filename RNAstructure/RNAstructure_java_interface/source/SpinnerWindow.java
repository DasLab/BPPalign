/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JComponent;
import javax.swing.JSpinner;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;

/**
 * A class that holds a window which selects a specific value from a spinner.
 * @author Jessica Reuter
 */
public class SpinnerWindow extends DialogWindow {
    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = 3335586373451425783L;

    /**
     * The spinner this window holds
     */
    private JSpinner spinner;

    /**
     * Constructor
     * @param min  the minimum spinner value
     * @param max  the maximum spinner value
     * @param step  the step in the spinner values
     * @param start  the default spinner value
     * @param title  the title of the window
     * @param ok  the action on the "OK" button
     * @param cancel  the action on the "Cancel" button
     */
    public SpinnerWindow( int min, int max, int step, int start, String title,
                          ActionBase ok, ActionBase cancel ) {
        setTitle( title );

        SpinnerModel model = new SpinnerNumberModel( start, min, max, step );
        spinner = new JSpinner( model );
        ((JSpinner.DefaultEditor)spinner.getEditor()).getTextField()
            .setEditable( false );

        ButtonPanel buttonPanel = new ButtonPanel( 2 );
        buttonPanel.setAction( ok, 0, "OK" );
        buttonPanel.setAction( cancel, 1, "Cancel" );
        buttonPanel.addButtons();

        setInsets( 10, 0, 0, 0 );
        setPad( 20, 10 );
        placeComponent( 0, 0, spinner );

        setInsets( 0, 0, 0, 0 );
        setPad( 0, 0 );
        placeComponent( 0, 1, buttonPanel );
    }

    /**
     * Get the spinner this window holds.
     * @return  the spinner
     */
    public JSpinner getSpinner() { return spinner; }

    /**
     * An inner class which allows a SpinnerWindow to be constructed from the
     * more generic context of the RNAWindowFactory. This particular factory
     * creates a window used to set the number of colors in a dot plot.
     * @author Jessica Reuter
     */
    public static class ColorFactory implements RNAWindowFactory {
        /**
         * Create a SpinnerWindow for color selection.
         * @return  a new SpinnerWindow
         */
        public SpinnerWindow createWindow() {
            ChainedAction colors = new ChainedAction(
                PlotChanger.makeColorAction(), new CloseAction() );
            int entries =
		((DrawDotPlot)((Imager)RNAstructure.getCurrentWindow())
                .getSketcher()).getEntries();

            return new SpinnerWindow( 3, 15, 1, entries, "Choose Colors",
                colors, new CloseAction() );
        }
    }

    /**
     * An inner class which allows a SpinnerWindow to be constructed from the
     * more generic context of the RNAWindowFactory. This particular factory
     * creates a window which is used to set the structure viewed.
     * @author Jessica Reuter
     */
    public static class StructureChooserFactory implements RNAWindowFactory {
        /**
         * Create a SpinnerWindow for structure selection.
         * @return  a new SpinnerWindow
         */
        public SpinnerWindow createWindow() {
            final int structures =
                RNAstructure.getCurrentStrand().GetStructureNumber();

            ActionBase structureAction = new ActionBase() {
                public void act() {
                    String value = 
                        ((SpinnerWindow)((JComponent)event.getSource())
                        .getTopLevelAncestor()).getSpinner().getValue()
                        .toString();
                    int number = Integer.parseInt( value );

                    DrawStructure sketch = ((DrawStructure)((Imager)
			RNAstructure.getCurrentWindow()).getSketcher());
                    sketch.setStructureNumber( number );
                    sketch.getImagePanel().repaint();
                }
            };

            ChainedAction spinAction =
                new ChainedAction( structureAction, new CloseAction() );

            int current = ((DrawStructure)((Imager)
                RNAstructure.getCurrentWindow()).getSketcher())
                .getCurrentStructure();

            return new SpinnerWindow( 1, structures, 1, current,
                "Select Structure", spinAction, new CloseAction() );
        }
    }

    /**
     * An inner class which allows a SpinnerWindow to be constructed from the
     * more generic context of the RNAWindowFactory. This particular factory
     * creates a window which is used to set the zoom in an image.
     * @author Jessica Reuter
     */
    public static class ZoomFactory implements RNAWindowFactory {
        /**
         * Create a SpinnerWindow for zooming.
         * @return  a new SpinnerWindow
         */
        public SpinnerWindow createWindow() {
            Sketcher sketch = ((Imager)RNAstructure.getCurrentWindow())
                .getSketcher();
            int scale = (int)(sketch.getScale() * 100);

            ChainedAction zoom =
                new ChainedAction( new ZoomHandler(),
                new CloseAction() );

            return new SpinnerWindow( 5, 500, 1, scale, "Image Zoom",
                zoom, new CloseAction() );
        }
    }
}
