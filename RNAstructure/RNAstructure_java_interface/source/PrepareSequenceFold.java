/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JButton;

/**
 * A class that spawns a FoldSingleWindow from a new sequence window.
 * @author Jessica Reuter
 */
public class PrepareSequenceFold extends ActionBase {
    /**
     * The type of nucleic acid to fold
     */
    private String acid;

    /**
     * Constructor
     * @param acid  the type of nucleic acid
     */
    public PrepareSequenceFold( String acid ) { this.acid = acid; }

    /**
     * Spawn the folding window.
     */
    protected void act() {
        NewSequenceWindow window =
	    (NewSequenceWindow)RNAstructure.getCurrentWindow();
        boolean goToFold = true;

        String file = window.getTitle();
        if( file.equals( "New Sequence" ) ) {
            SequenceSaver saver = new SequenceSaver( true );
            saver.act();
            file = saver.getSaveFile();
            if( saver.isError() ) { goToFold = false; }
        }

        if( goToFold ) {
            FoldSingleWindow fold =
                new FoldSingleWindow.Factory( acid ).createWindow();
            fold.setFileNames( file, file.replace( ".seq", ".ct" ) );
            fold.pack();
            fold.setLocationRelativeTo( RNAstructure.getFrame() );
            fold.setVisible( true );
            RNAstructure.setCurrentWindow( fold );

            JButton source = fold.inputPanel.getButtons()[0];

            ChainedAction foldAction =
                (ChainedAction)source.getActionListeners()[0];
            ActionBase[] chain = foldAction.getChain();
            ActionBase.event.setSource( source ); 
            for( int i = 1; i < chain.length; i++ ) { chain[i].act(); }
        }
    }
}
