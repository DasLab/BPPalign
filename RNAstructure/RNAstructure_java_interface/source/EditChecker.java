/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that monitors a sequence to see if it's been edited.
 * @author Jessica Reuter
 */
public class EditChecker extends WindowActionBase {
    /**
     * Check a window for edits
     */
    public void closeWindowDoDuring() {
	final NewSequenceWindow window = (NewSequenceWindow)event.getSource();

        // If the window has been edited, ask the user about saving it
        if( window.isEdited() ) {
            boolean isNew = window.getTitle().equals( "New Sequence" );

	    ActionBase closer = new ActionBase() {
		public void act() {
		    window.dispose();
		}
	    };

            ChainedAction save = new ChainedAction(
                new CloseAction(), new SequenceSaver( isNew ), closer );

            ChainedAction discard = new ChainedAction(
                new CloseAction(), closer );

            // Query the user about saving or discarding the sequence
            RNAstructureInfoDialog.choose( "The sequence has been " +
                "modified.\nSave Changes?", save, discard );
        } else { window.dispose(); }
    }
}
