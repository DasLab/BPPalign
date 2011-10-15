/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 *  
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dialog;
import java.util.Iterator;
import java.util.Map;

import javax.swing.JTextPane;
import javax.swing.text.BadLocationException;
import javax.swing.text.Style;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyleContext;
import javax.swing.text.StyledDocument;

/**
 * A class for windows inside the RNAstructure GUI which display color codes
 * for various annotation keys.
 * @author Jessica Reuter
 */
public class AnnotationKeyWindow extends DialogWindow {
    /**
     * The JTextPane that this window uses to display its annotation key.
     */
    private JTextPane pane;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -8688294811937117629L;

    /**
     * Constructor
     * @param title  the title of the color key
     */
    public AnnotationKeyWindow( String title ) {
        setModalityType( Dialog.ModalityType.APPLICATION_MODAL );
        setTitle( title );

        initialize();
        createKey( title );
        finalizeKey();
    }

    /**
     * Add a key to the window, highlighted with a particular color
     * @param text  the text of the key
     * @param color the color to highlight the text with
     */
    private void addKey( String text, Color color ) {
        try {
            // Create appropriate style for addition
            StyledDocument doc = pane.getStyledDocument();
            Style def = StyleContext.getDefaultStyleContext().
                getStyle( StyleContext.DEFAULT_STYLE );

            Style style = doc.addStyle( color.toString(), def );
            StyleConstants.setForeground( style, color );
            StyleConstants.setFontSize( style, 24 );
            StyleConstants.setBold( style, true );

            // Add the key to the pane
            doc.insertString( doc.getLength(), text + "\n", style );
        } catch( BadLocationException e ) {
            RNAstructureInfoDialog.error( "Error creating annotation key." );
        }
    }

    /**
     * Create the key by adding to the frame.
     * @param title  the title of the key
     */
    private void createKey( String title ) {
        // Add entries to the key, depending on which key to construct
        Map<String, Color> map = 
            title.equals( "Probability Annotation Key" ) ?
                AnnotationMapHolder.getChanceMap() :
                AnnotationMapHolder.getSHAPEMap();
        Iterator<Map.Entry<String, Color>> iterator =
            map.entrySet().iterator();

        while( iterator.hasNext() ) {
            Map.Entry<String, Color> current = iterator.next();
            addKey( current.getKey(), current.getValue() );
        }
    }

    /**
     * Initialize the text pane and dialog.
     */
    private void initialize() {
        // Set layout and create text pane
        setLayout( new BorderLayout() );
        pane = new JTextPane();
        pane.setBackground( Color.WHITE );
        pane.setEditable( false );
        add( pane );
    }

    /**
     * Create the final packing and viewing of the window.
     */
    private void finalizeKey() {
        setSize( 400, 300 );
        setResizable( false );
    }

    /**
     * An inner class which allows an AnnotationKeyWindow to be created from
     * the more generic context of the RNAWindowFactory.
     * @author  Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * The title of the window
         */
        private String title;

        /**
         * Constructor
         * @param title  the title of the window
         */
        public Factory( String title ) { this.title = title; }

        /**
         * Create a new AnnotationKeyWindow.
         * @return  the new window
         */
        public Object createWindow() {
            return new AnnotationKeyWindow( title );
        }
    }
}