/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.BorderLayout;
import java.awt.CardLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridLayout;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

/**
 * A class that holds a window which contains the results of an Oligowalk
 * calculation, as well as various ways to access them.
 * @author Jessica Reuter
 */
public class OligoWalkResultsWindow extends ResultsWindow {
    /**
     * The array of colors used for the bars inside the graph on this window
     */
    private final Color[] colors = {
        Color.BLUE,                 // Overall
        new Color( 114, 67, 212 ),  // Target Break
        new Color( 17, 167, 77 ),   // Duplex
        Color.MAGENTA,              // Self-Oligo
        Color.BLACK,                // Temperature
        new Color( 139, 25, 14 )    // Oligo-Oligo
    };

    /**
     * The main panel which holds the graph
     */
    private JPanel graph;

    /**
     * The array of the current and maximum index selected
     */
    private int[] indices;
    
    /**
     * The array of labels at the top of the window
     */
    private JLabel[] labels;

    /**
     * The most stable oligo
     */
    private int mostStable;

    /**
     * The length of the oligos
     */
    private final int oligoLength;

    /**
     * The Oligowalk_object this window is describing
     */
    private Oligowalk_object strand;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = 5731622965900211137L;

    /**
     * Constructor
     * @param strand  the Oligowalk_object this window describes
     * @param oligoLength  the oligo length
     * @param concentration  the oligo concentration
     * @param title  the title of the window
     */
    public OligoWalkResultsWindow( Oligowalk_object strand, int oligoLength,
                                   String concentration, String title ) {
	// Set menus and sections
	setSections( 1, 3, 4, 5, 7 );
	SmartMenuBar bar = RNAstructure.getBar();
	setMenus( bar.getMenuMaker().createOligoMenu() );
	bar.refresh( sections, menus );

	// Set titles, strand, oligo length, and most stable oligo
	setTitles( title );
	setResizable( false );

	this.strand = strand;
	this.oligoLength = oligoLength;
	int size = strand.GetSequenceLength() - oligoLength + 1;

	double maxOligoValue = Double.MAX_VALUE;
	for( int i = 1; i <= size; i++ ) {
	    double next =
		OligoWalkResultsWindow.getOligoValue( strand, 3, i );
	    if( next < maxOligoValue ) {
		maxOligoValue = next;
		mostStable = i;
	    }
	}

	// Create labels and indices arrays
        labels = new JLabel[10];
        
        indices = new int[2];
        indices[0] = 1;
        indices[1] = size;

        // Create the panel which holds oligo concentration and current number
        JPanel numberPanel = new JPanel( new GridLayout( 2, 1 ) );
        labels[0] = createLabel( "Oligo: 1" );
        numberPanel.add( labels[0] );

        numberPanel.add( createLabel(
	    "Oligo Concentration: " + concentration ) );

        // Create labels that hold specific oligo deltaG data and their panel
        String deltaGString = "\u0394G\u00B037";
        String[] gData = { "Overall", "Break Target", "Duplex",
            "Oligo-Self", "Tm", "Oligo-Oligo" };
        
        JPanel deltaGPanel = new JPanel( new GridLayout( 3, 2 ) );
        for( int i = 1; i <= 6; i++ ) {
            String text = gData[i-1] + " " + deltaGString + " = 0.0";
            labels[i] = createLabel( text, colors[i-1] );
            deltaGPanel.add( labels[i] );
        }
        
        // Create the encompassing panel for the above label panels
        JPanel mainLabelPanel = new JPanel();
        mainLabelPanel.setBackground( Color.white );
        mainLabelPanel.add( numberPanel );
        mainLabelPanel.add( deltaGPanel );
        
        // Create panel that moves oligo bar graph
        JPanel buttonPanel = new JPanel( new GridLayout( 1, 5 ) );
        Object[] buttons =
	    { "<<", -10, "<", -1, "Go...", 0, ">", 1, ">>", 10 };
        
        for( int i = 1; i <= 5; i++ ) {
            JButton button = new JButton( buttons[(2*i)-2].toString() );
	    button.setOpaque( true );
            button.setBackground( Color.LIGHT_GRAY );
            
            int increment = (Integer)buttons[(2*i)-1];
            button.addActionListener(
                GraphListener.createIncrementMover( increment ) );
            
            JPanel mini = new JPanel();
	    mini.setBackground( Color.WHITE );
            mini.setBorder( BorderFactory.createEmptyBorder( 5, 5, 5, 5 ) );
            mini.add( button );
            buttonPanel.add( mini );
        }

	// Create two empty labels to hold the graph range
        labels[8] = createLabel( "" );
        labels[9] = createLabel( "" );
        for( int i = 8; i <= 9; i++ ) {
            labels[i].setFont( new Font( Font.MONOSPACED, 0, 16 ) );
            labels[i].setHorizontalAlignment( JLabel.RIGHT );
        }

        
        // Create panel that holds oligo bar graph
        graph = new JPanel( new CardLayout() );

        StackedOligoGraphPanel one =
	    new StackedOligoGraphPanel( strand, oligoLength );
        one.createStackedGraph();
        graph.add( one, "and" );
        
        for( int i = 1; i <= 6; i++ ) {
	    OligoGraphPanel next = new OligoGraphPanel( strand, oligoLength );
	    if( i != 5 ) {
	        next.createGraph( i, colors[i-1] );
	        graph.add( next, Integer.toString( i ) );
	    }
        }
                
        // Create labels which hold the sequence and oligos
        JLabel seqLabel = createLabel( "" );
        seqLabel.setFont( new Font( Font.MONOSPACED, 0, 16 ) );
        seqLabel.setHorizontalAlignment( JLabel.LEFT );
        StringBuilder seq = new StringBuilder( "<html>&nbsp;" );

	int fullSize = strand.GetSequenceLength();
	for( int i = 1; i < fullSize; i++ ) {
	    char base = strand.GetNucleotide( i + 1 );
	    if( strand.GetErrorCode() != 0 ) {
		RNAstructureInfoDialog.error(
	            strand.GetErrorMessage( strand.GetErrorCode() ) );
		return;
	    } else {
		if( strand.GetPair( i ) != 0 ) {
		    seq.append( "<FONT COLOR=\"FF0000\">" + base + "</FONT>" );
		} else { seq.append( base ); }
	    }
	}

        seqLabel.setText( seq.toString() );
        
        JLabel oligoLabel = createLabel( "" );
        oligoLabel.setFont( new Font( Font.MONOSPACED, 0, 16 ) );
        oligoLabel.setHorizontalAlignment( JLabel.LEFT );
        labels[7] = oligoLabel;
        setOligo( 1 );

        // Create the panel for the sequence and oligo labels
        JPanel nucLabelPanel = new JPanel( new GridLayout( 2, 1 ) );
        nucLabelPanel.add( oligoLabel );
        nucLabelPanel.add( seqLabel );
        
        // Create spacers for the edges
        JLabel[] spacers = new JLabel[2];
        for( int i = 0; i < 2; i++ ) {
            JLabel label = createLabel( "" );
            label.setPreferredSize( new Dimension( 10, 500 ) );
            spacers[i] = label;
        }
        
        // Put together the main panel and add it to scroll pane
        JPanel mainPanel = new JPanel( new BorderLayout() );
        mainPanel.add( nucLabelPanel, BorderLayout.NORTH );
        mainPanel.add( spacers[0], BorderLayout.WEST );
        mainPanel.add( spacers[1], BorderLayout.EAST );
        mainPanel.add( graph );
        
        JScrollPane scroller = new JScrollPane( mainPanel );
        Dimension graphSize = new Dimension( 875, 565 );
        scroller.setPreferredSize( graphSize );
        scroller.setMaximumSize( graphSize );
        scroller.setAutoscrolls( true );
        
        // Create the bounds panel
        JPanel boundsPanel = new JPanel( new BorderLayout() );
        boundsPanel.setBackground( Color.white );
        boundsPanel.add( labels[8], BorderLayout.NORTH );
        
        JLabel spacer = createLabel( "" );
        spacer.setPreferredSize( new Dimension( 15, 470 ) );
        boundsPanel.add( spacer );
        boundsPanel.add( labels[9], BorderLayout.SOUTH );
        
        // Create the completed graph panel with bounds label added
        JPanel completePanel = new JPanel();
        completePanel.setBackground( Color.white );
        completePanel.add( boundsPanel );
        completePanel.add( scroller );
        
        // Place components in their proper places
        add( mainLabelPanel, BorderLayout.NORTH );
        add( buttonPanel );
        add( completePanel, BorderLayout.SOUTH );

	// Set the labels to the first oligo
	setLabels( 1 );
    }

    /**
     * Create a centered label on a white background with the default text
     * color, black.
     * @param text  the text on the label
     * @return  the label
     */
    private JLabel createLabel( String text ) {
    	return createLabel( text, Color.BLACK );
    }
    
    /**
     * Create a centered label on a white background with the specified text
     * color.
     * @param text  the text on the label
     * @param color  the color of the label
     * @return  the label
     */
    private JLabel createLabel( String text, Color color ) {
        JLabel label = new JLabel( text );
        label.setBackground( Color.WHITE );
        label.setForeground( color );
        label.setHorizontalAlignment( JLabel.CENTER );
        label.setOpaque( true );
        label.setFont( new Font( label.getFont().getName(), 0, 24 ) );
        return label;
    }
    
    /**
     * Get the main graph panel
     * @return  the label
     */
    public JPanel getGraphPanel() { return graph; }

    /**
     * Get the index of the most stable oligo
     * @return  the index
     */
    public int getMostStableOligo() { return mostStable; }

    /**
     * Get the oligo length.
     * @return  the oligo length
     */
    public int getOligoLength() { return oligoLength; }
    
    /**
     * Get the value of an oligo, based on its type.
     * @param result  the result strand the value comes from
     * @param type  the type of data to get
     * @param index  the index of the oligo
     * @return  the value
     */
    public static double getOligoValue( Oligowalk_object result, int type,
					int current ) {
	double newValue =
	    ( type == 1 ) ? result.GetOverallDG( current ) :
	    ( type == 2 ) ? result.GetBreakTargetDG( current ) :
	    ( type == 3 ) ? result.GetDuplexDG( current ) :
	    ( type == 4 ) ? result.GetOligoSelfDG( current ) :
	    ( type == 5 ) ? result.GetTm( current ) :
	    ( type == 6 ) ? result.GetOligoOligoDG( current ) :
	    0.0;

	if( result.GetErrorCode() != 0 ) {
	    RNAstructureInfoDialog.error(
		result.GetErrorMessage( result.GetErrorCode() ) );
	    return 1000000;
	}

	if( type >= 1 && type <= 6 ) { return newValue; }
	else { throw new IllegalArgumentException( "Invalid result type." ); }
    }


    /**
     * Get the visible panel of the graph, ie. the graph currently showing.
     */
    public OligoGraphPanel getVisiblePanel() {
    	Component[] components = graph.getComponents();
    	for( Component element: components ) {
	    if( element.isVisible() ) { return (OligoGraphPanel)element; }
    	}
    	return (OligoGraphPanel)components[0];
    }
    
    /**
     * Set the labels above the graph, depending on which oligo is selected.
     * @param  current  the selected oligo
     */
    public void setLabels( int current ) {
        String idxLabel = labels[0].getText();
        String oldIndex =
            idxLabel.substring( idxLabel.lastIndexOf( ":" ) + 2 ).trim();
        idxLabel =
            idxLabel.replace( oldIndex, Integer.toString( current ) );
        labels[0].setText( idxLabel );
        
        for( int i = 1; i <= 6; i++ ) {
            String text = labels[i].getText();
            double newValue =
		OligoWalkResultsWindow.getOligoValue( strand, i, current );
            
	    if( newValue != 1000000 ) {
		String oldValue =
		    text.substring( text.lastIndexOf( "=" ) + 2 );
		text = text.replace( oldValue, Double.toString( newValue ) );
		labels[i].setText( text );
	    } else return;
        }

	double[] bounds = getVisiblePanel().getGraphBounds();
	labels[9].setText( "<html>" + Double.toString( bounds[0] ) );
	labels[8].setText( "<html>" + Double.toString( bounds[1] ) );
    }
    
    /**
     * Set the current oligo on the graph.
     * @param index  the index of the oligo
     */
    public void setOligo( int index ) {
	int adjustedLength = strand.GetSequenceLength() - oligoLength + 1;
        if( adjustedLength >= index ) {
            StringBuilder fwd = new StringBuilder();
            for( int i = 0; i < oligoLength - 1; i++ ) {
		char base = strand.GetNucleotide( index + i + 1 );
		if( strand.GetErrorCode() != 0 ) {
		    RNAstructureInfoDialog.error(
			strand.GetErrorMessage( strand.GetErrorCode() ) );
		    return;
		} else { fwd.append( base ); }
	    }
            
            StringBuilder rev = new StringBuilder();
            for( int i = 1; i < index; i++ ) { rev.append( " " ); }

            rev.append( "3" );
            for( int i = 0; i < oligoLength - 1; i++ ) {
                char next = fwd.charAt( i );
                
                switch( next ) {
                    case 'A':
                    	if( strand.getIsrna() ) { rev.append( 'U' ); }
                    	else { rev.append( 'T' ); }
                    	break;
                    case 'C':
                    	rev.append( 'G' );
                    	break;
                    case 'G':
                    	rev.append( 'C' );
                    	break;
                    case 'T':
                    case 'U':
                    	rev.append( 'A' );
                    	break;
                    case 'X':
                    	rev.append( 'X' );
                    	break;
                    default:
                    	RNAstructureInfoDialog.error( "Not a base: " + next );
			return;
                };
            }
            rev.append( "5" );
            labels[7].setText( rev.toString() );
        } else throw new IllegalArgumentException( "Oligo doesn't exist." );
    }   
}
