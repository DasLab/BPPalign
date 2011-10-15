/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JButton;
import javax.swing.JDialog;

/**
 * A class that handles the closing of dialogs by clicking on a button.
 * @author Jessica Reuter
 */
public class CloseAction extends ActionBase {
    /**
     * Close the dialog on which a button panel is placed.
     */
    protected void act() {
        ((JDialog)((JButton)event.getSource()).getTopLevelAncestor())
            .dispose();
    }
}