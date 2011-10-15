/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.Point;

import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

/**
 * A class that handles viewing a variety of images created by the RNAstructure
 * GUI from its calculations and data files.
 * @author Jessica Reuter
 */
public class Imager extends ResultsWindow {
    /**
     * The file this imager comes from
     */
    private String file;

    /**
     * A label which holds information on this particular image.
     */
    private JLabel infoLabel;

    /**
     * The panel that holds the legend.
     */
    private JPanel legendGrid;

    /**
     * The scroll pane that holds the main image.
     */
    private JScrollPane scroller;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = 1297732886935988738L;

    /**
     * The sketcher that gave rise to the image on this imager
     */
    private SingleSketcher sketcher;

    /**
     * Constructor
     * @param file  the file this image comes from
     * @param legend  the legend for this image
     * @param panel  the panel that holds an image
     */
    public Imager( String file, Object[][] legend, PrintablePanel panel ) {
        this.file = file;

        // Create the info label.
	JPanel infoPanel = new JPanel();
	infoLabel = new JLabel( "Imager" );
        infoLabel.setHorizontalAlignment( JLabel.CENTER );
	infoLabel.setVerticalAlignment( JLabel.CENTER );
	infoPanel.setPreferredSize( new Dimension( 600, 20 ) );
        infoPanel.add( infoLabel );
	getContentPane().add( infoPanel, BorderLayout.NORTH );

        // Add image panel to a scroll pane, then the pane to the main panel
        scroller = new JScrollPane( panel );
        scroller.setPreferredSize( new Dimension( 600, 600 ) );
        scroller.getViewport().setViewPosition( new Point( 600, 600 ) );
        resultsPanel.add( scroller, BorderLayout.CENTER );

        // Create the legend panel, then add it to the main panel
        if( legend != null ) {
            legendGrid = new JPanel();
            repaintLegend( legend );

            JScrollPane scroller2 = new JScrollPane( legendGrid );
            scroller2.setHorizontalScrollBarPolicy(
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER );
            scroller2.setPreferredSize( new Dimension( 600, 110 ) );
            resultsPanel.add( scroller2, BorderLayout.SOUTH );
        }

        // Add a zoom listener
        scroller.addKeyListener( new KeyZoomHandler() );
	
        // Add a focus listener for movement from the keyboard
        addFocusListener( new FocusAdapter() {
	    public void focusGained( FocusEvent e ) {
                scroller.requestFocusInWindow();
	    }
	});
    }

    /**
     * Get the file this imager came from.
     * @return  the file
     */
    public String getFile() { return file; }

    /**
     * Get the legend panel from this imager.
     * @return  the legend panel
     */
    public JPanel getLegendPanel() { return legendGrid; }

    /**
     * Get the main scroll pane from this imager
     */
    public JScrollPane getMainScroller() { return scroller; }

    /**
     * Get the sketcher from this imager
     * @return  the sketcher
     */
    public SingleSketcher getSketcher() { return sketcher; }

    /**
     * Repaint the legend when it changes.
     * @param legend  the legend text array
     */
    public void repaintLegend( Object[][] legend ) {
        int length = legend.length;

        legendGrid.removeAll();
        legendGrid.setLayout( new GridLayout( length, 1 ) );
        legendGrid.setPreferredSize( new Dimension( 600, 20 * length ) );
        for( int i = 0; i < length; i++ ) {
            JLabel label = new JLabel( legend[i][0].toString() );
            label.setHorizontalAlignment( JLabel.CENTER );
            label.setOpaque( true );
            label.setBackground( Color.white );
            label.setForeground( (Color)legend[i][1] );
            legendGrid.add( label );
        }
    }

    /**
     * Set defaults for imager viewing and menus
     * @param file  the file name the imager came from
     * @param menus  the menus a particular imager uses
     */
    public void setDefaults( SmartMenu... newMenus ) {
        setTitles( file );
        RNAstructure.getToolBar().setButtonsEnabled( false, true );
        setSections( 1, 3, 4, 5, 6, 7 );
        SmartMenuBar bar = RNAstructure.getBar();
        setMenus( newMenus );
        bar.refresh( sections, menus );
        RNAstructure.getBar().enableMenus();
    }

    /**
     * Set the information label.
     * @param text  the text to set on the label.
     */
    public void setInfoLabel( String text ) { infoLabel.setText( text ); }

    /**
     * Set the sketcher for this imager.
     * @param sketcher  the sketcher to set
     */
    public void setSketcher( SingleSketcher sketcher ) {
        this.sketcher = sketcher;
    }

    /**
     * An inner class which allows an Imager using a dot plot to be constructed
     * from the more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class DotPlotFactory implements RNAWindowFactory {
        /**
         * The DrawDotPlot sketcher being drawn
         */
        private DrawDotPlot plotter;

        /**
         * A constant value which designates an unknown plot type
         */
        public static final int UNKNOWN_PLOT = 0;

        /**
         * A constant value which designates an energy plot
         */
        public static final int ENERGY_PLOT = 1;

        /**
         * A constant value which designates a probability plot
         */
        public static final int PROBABILITY_PLOT = 2;

        /**
         * The type of plot being created; this value is compared to one of the
         * constant plot types
         */
        private int type;

        /**
         * Dot Plot Selection Constructor
         * @param type  the type of dot plot that is made
         */
        public DotPlotFactory( int type ) {
            plotter = null;
            this.type = type;
        }

        /**
         * Existing Dot Plot Constructor
         * @param action  the DrawDotPlot action to use
         */
        public DotPlotFactory( DrawDotPlot action ) {
            plotter = action;
            type = 0;
        }

        /**
         * Create an imager that displays a dot plot.
         * @return  the imager
         */
        public Object createWindow() {
            if( plotter == null || plotter.getFile().equals( "" ) ) {
                plotter =
                    ( type == ENERGY_PLOT ) ?
                        DrawDotPlot.drawEnergyPlot( "", 5 ) :
                    ( type == PROBABILITY_PLOT ) ?
                        DrawDotPlot.drawProbabilityPlot( "", 5 ) :
                    DrawDotPlotFromText.draw();
                plotter.act();
            }

            String file = plotter.getFile();
            if( !file.equals( "" ) ) {
                plotter.resetFile();
                Object[][] key = plotter.getLegend();
                PrintablePanel panel = plotter.getImagePanel();
                Imager imager = new Imager( file, key, panel );
                imager.setSketcher( plotter );

                MenuMaker maker = RNAstructure.getBar().getMenuMaker();
		SmartMenu plotMenu = ( type == ENERGY_PLOT ) ?
		    maker.createPlotOutputFoldingMenu() :
		    maker.createPlotOutputPartitionMenu();

		imager.setInfoLabel( "DOT PLOT" );
                imager.setDefaults( maker.createDrawMenu(), plotMenu );
                plotter.getImagePanel().addMouseListener(
                    new DotAnalyzer( plotter ) );

                return imager;
            } else return null;
        }
    }

    /**
     * An inner class which allows an Imager to be constructed from the more
     * generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class StructureFactory implements RNAWindowFactory {
        /**
         * The DrawStructure sketcher being drawn
         */
        private DrawStructure action;

        /**
         * Default constructor
         * Used when the DrawStructure instance for this factory must still be
         * selected manually
         */
        public StructureFactory() { action = null; }

        /**
         * DrawStructure constructor
         * Used when a particular specified structure is to be drawn.
         * @param draw  the sketcher to be drawn
         */
        public StructureFactory( DrawStructure draw ) { action = draw; }

        /**
         * Create an imager that shows a structure.
         * @return  the imager
         */
        public Object createWindow() {
            DrawStructure sketcher =
                ( action == null ) ? new DrawStructure( "" ) : action;
	    sketcher.create();

            String file = sketcher.getFile();
            if( !file.equals( "" ) ) {
		RNA current = new RNA( file, 1 );
		int structures = current.GetStructureNumber();
		boolean pseudoknot = false;
		for( int i = 1; i <= structures; i++ ) {
		    if( ( pseudoknot = current.ContainsPseudoknot( i ) ) ) {
			break;
		    }
		}

		if( pseudoknot ) {
		    sketcher = new DrawStructureCircular( "" );
		    sketcher.setFile( file );
		    sketcher.create();
		    ((DrawStructureCircular)sketcher).prepare();
		} else {
		    ((DrawStructure)sketcher).prepare();
		}

                if( !sketcher.containsPairs() ) {
                    RNAstructureInfoDialog.error(
                        "This structure contains no pairs." );
                    return null;
                }

                sketcher.resetFile();
                Object[][] key = sketcher.getLegend();
                PrintablePanel panel = sketcher.getImagePanel();

                Imager imager = new Imager( file, key, panel );
                imager.setSketcher( sketcher );
                imager.getMainScroller()
                    .addKeyListener( new StructureMover() );

                MenuMaker maker = RNAstructure.getBar().getMenuMaker();
                imager.setDefaults(
                    maker.createStructureAnnotationMenu(),
                    maker.createStructureMenu() );

                return imager;
            } else return null;
        }
    }

    /**
     * An inner class which constructs two separate but related Imagers.
     * @author Jessica Reuter
     */
    public static class TwoImagersFactory extends ActionBase {
	/**
	 * The first file from which an image is drawn
	 */
	private String file1;

	/**
	 * The second file from which an image is drawn
	 */
	private String file2;

	/**
	 * Default constructor
	 * Used for Dynalign dot plots
	 */
	public TwoImagersFactory() {
	    file1 = null;
	    file2 = null;
	}

	/**
	 * Constructor
	 * @param file1  the first file
	 * @param file2  the second file
	 */
	public TwoImagersFactory( String file1, String file2 ) {
	    this.file1 = file1;
	    this.file2 = file2;
	}

        /**
         * Create two imagers
         */
        public void act() {
	    RNAWindowFactory one = null;
	    RNAWindowFactory two = null;

	    if( file1 != null && file2 != null ) {
		one = new StructureFactory( new DrawStructure( file1 ) );
		two = new StructureFactory( new DrawStructure( file2 ) );
	    } else {
		FileSelector selector = new FileSelector( false, "open",
		    "Dynalign Save Files,dsv" );
		String file = selector.showSelector();

		if( !file.equals( "" ) ) {
		    DrawDotPlotFromDynalignSave first =
			DrawDotPlotFromDynalignSave.drawFirstPlot( file );
		    first.act();

		    DrawDotPlotFromDynalignSave second =
			DrawDotPlotFromDynalignSave.drawSecondPlot( file );
		    second.act();

		    one = new DotPlotFactory( first );
		    two = new DotPlotFactory( second );
		}
	    }

	    if( one != null && two != null ) {
		new WindowDisplayer( one ).act();
		new WindowDisplayer( two ).act();
	    }
        }
    }
}
