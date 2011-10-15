/* 
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import javax.swing.JButton;
import javax.swing.JComponent;

/**
 * A base class for all actions that are spawned from start buttons.
 * @author Jessica Reuter
 */
public abstract class StrandStarter extends ActionBase {
    /**
     * Execute an action that starts on a start button
     */
    protected void act() {
        new Thread() {
            public void run() {
		try {
		    InputWindow window =
			(InputWindow)(((JComponent)event.getSource())
			.getTopLevelAncestor());

		    if( window.isReady() ) {
			window.setData();

			window.setEnabled( false );
			RNAstructure.getBar().disableVariableMenus();

			RNAProgressDialog runner = null;

			if( window.monitorProgress() ) {
			    runner = new RNAProgressDialog();

			    if( !window.monitorProgressIndeterminate() ) {
				runner.runDialog();
			    }
			}

			execute();
			
			if( window.monitorProgress() &&
			    !window.monitorProgressIndeterminate() ) {
			    runner.setDone();
			}

			window.setEnabled( true );
		    }
		} catch( Exception ex ) {
		    RNAstructureInfoDialog.error(
                        "Error during module initialization." );
		}
            }
        }.start();
    }

    /**
     * Execute a particular action. This action is only implemented in the
     * subclasses of this class.
     */
    protected abstract void execute();
}
