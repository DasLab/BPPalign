/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.print.PageFormat;
import java.awt.print.Printable;
import java.awt.print.PrinterException;
import java.awt.print.PrinterJob;

import javax.swing.JPanel;

/**
 * A class which holds a panel whose contents are capable of being printed.
 * @author Jessica Reuter
 */
public class PrintablePanel extends JPanel implements Printable {
    /**
     * The job that prints out this panel
     */
    private PrinterJob job = PrinterJob.getPrinterJob();

    /**
     * A serialized ID, required when extending GUI classes
     */
    private static final long serialVersionUID = -489847909728935586L;

    /**
     * Get the current print job.
     * @return  the print job
     */
    public PrinterJob getPrintJob() { return job; }

    /**
     * Render the panel for the print job (inherited from Printable)
     * @param graphics  the graphics object
     * @param format  the page format
     * @param index  the page index
     */
    public int print( Graphics graphics, PageFormat format, int index ) {
        // The panel only takes up one page, so if the index says we're on page
        // 2, something went wrong.
        if ( index > 0 ) { return NO_SUCH_PAGE; }

        // Print the panel
        Graphics2D g2d = (Graphics2D)graphics;
        g2d.translate( format.getImageableX(), format.getImageableY() );
        printAll( graphics );

        return PAGE_EXISTS;
    }

    /**
     * Print if the user decides to, then reset the printer job properties
     */
    public void printPanel() {
        job.setPrintable( this );
        if( job.printDialog() ) {
            try { 
                job.print();
                job = PrinterJob.getPrinterJob();
            } catch( PrinterException e ) {
                e.printStackTrace();
            }
        }
    }
}