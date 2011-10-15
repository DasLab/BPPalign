/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.text.DefaultEditorKit;

public class MovementAction extends ActionBase {
    /**
     * The type of action being done
     */
    private final int type;

    /**
     * Constructor, private to ensure no edit of it
     * @param type  the type of action
     */
    private MovementAction( int type ) { this.type = type; }

    /**
     * Activate a movement action.
     */
    protected void act() {
        RNAstructure.setLabel( "For help, press F1." );
        
        switch( type ) {
            case 0:
                new DefaultEditorKit.CutAction().actionPerformed( event );
                break;
            case 1:
                new DefaultEditorKit.CopyAction().actionPerformed( event );
                break;
            case 2:
                new DefaultEditorKit.PasteAction().actionPerformed( event );
        };
    }

    /**
     * Create a text copy action.
     */
    protected static MovementAction createTextCopy() {
        return new MovementAction( 1 );
    }

    /**
     * Create a text cut action.
     */
    protected static MovementAction createTextCut() {
        return new MovementAction( 0 );
    }

    /**
     * Create a text paste action.
     */
    protected static MovementAction createTextPaste() {
        return new MovementAction( 2 );
    }
}