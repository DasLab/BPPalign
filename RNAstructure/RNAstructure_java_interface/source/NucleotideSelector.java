/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that selects for nucleotides based on key typing.
 * @author Jessica Reuter
 */
public class NucleotideSelector extends KeyActionBase {
    /**
     * Select a nucleotide.
     */
    protected void typeKey() {
        char c = event.getKeyChar();

        boolean valid = c == 'a' || c == 'A' ||
                        c == 'c' || c == 'C' ||
                        c == 'g' || c == 'G' ||
                        c == 't' || c == 'T' ||
                        c == 'u' || c == 'U' ||
                        c == 'x' || c == 'X' ||
                        c == ' ';

        if( !valid ) { event.consume(); }
    }
}