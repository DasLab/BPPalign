/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.io.BufferedWriter;
import java.io.FileWriter;

/**
 * A class that writes different kinds of structure files.
 * @author Jessica Reuter
 */
public abstract class StructureFileWriter extends ActionBase {
    /**
     * The type of structure file being written
     */
    protected String type;

    /**
     * Write structure data to some form of output file.
     */
    protected void act() {
        RNAstructure.setLabel( "For help, press F1." );
        writeFile();
    }

    /**
     * Write a structure file.
     */
    protected abstract void writeFile();
}
