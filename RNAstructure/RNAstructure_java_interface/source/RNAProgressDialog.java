/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.Dimension;
import java.awt.Window;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;

import javax.swing.BorderFactory;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.SwingWorker;

/**
 * A class that creates progress dialogs for actions begun from start buttons.
 * @author Jessica Reuter
 */
public class RNAProgressDialog extends DialogWindow {
    /**
     * The progress bar
     */
    private JProgressBar bar;

    /**
     * A boolean that's true if monitoring is done, false if not
     */
    private boolean isDone;

    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -3220519658651768712L;

    /**
     * Constructor
     */
    public RNAProgressDialog() {
        setTitle( "Calculation in Progress..." );
	isDone = false;

        bar = new JProgressBar( 0, 100 );
        bar.setPreferredSize( new Dimension( 400, 50 ) );
        bar.setValue( 0 );

        JPanel panel = new JPanel();
        panel.add( bar );
        panel.setBorder( BorderFactory.createEmptyBorder( 10, 10, 10, 10 ) );

        placeComponent( 0, 0, panel );
        pack();
        setVisible( true );
    }

    /**
     * Run the progress dialog.
     * @param window  the window that gave input information that will be
     *                monitored with this dialog
     */
    public void runDialog() {
        SwingWorker<Void, Void> worker = new SwingWorker<Void, Void>() {
            public Void doInBackground() {
		try{
                    RNA strand = RNAstructure.getCurrentStrand();
                    ProgressMonitor monitor = new ProgressMonitor();
                    strand.SetProgress( new TProgressDialog( monitor ) );

                    setProgress( 0 );
                    int progress = 0;

                    while( !isDone ) {
                        try { Thread.sleep( 1000 ); }
                        catch( InterruptedException e ) {}

                        progress = monitor.getMonitorProgress();
                        setProgress( Math.min( progress, 100 ) );
                    }

                    strand.StopProgress();
		} catch( Exception e ) { e.printStackTrace(); }
                return null;
            }

            public void done() {
                ((Window)(bar.getTopLevelAncestor())).dispose();
            }
        };

        worker.addPropertyChangeListener( new PropertyChangeListener() {
            public void propertyChange( PropertyChangeEvent e ) {
                if ( "progress" == e.getPropertyName() ) {
                    int progress = (Integer)e.getNewValue();
                    bar.setValue( progress );
                } 
            }
        });

        worker.execute();
    }

    /**
     * Set the flag that this progress monitoring is done.
     */
    public void setDone() { isDone = true; }
}
