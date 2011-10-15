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
 * A class that creates a nucleic acid strand.
 * @author Jessica Reuter
 */
public class StrandCreator extends DataSelector {
    /**
     * The array of data needed to create a strand
     */
    protected Object[] data;

    /**
     * Constructor
     * @param data  the data needed to create a strand
     */
    private StrandCreator( Object... data ) { this.data = data; }

    /**
     * Create an action that allows for creation of a strand from a CT file.
     * @param index  the index of the button on the panel
     * @param isRNA  true if the strand is RNA, false if not
     * @return  the SEQ strand creator
     */
    public static StrandCreator createCTStrand( int index, boolean isRNA ) {
        return new StrandCreator( index, 1, isRNA, false );
    }

    /**
    * Create an action that allows for creation of a hybrid strand from two
    * SEQ files (This ends up making one big strand with a linker).
    * @param idx  the index of the button on the panel
    * @param rna  true if the strand is RNA, false if not
    * @return  the SEQ strand creator
    */
   public static StrandCreator createHybridSEQStrand( int idx, boolean rna ) {
       return new StrandCreator( idx, 2, rna, true );
   }

   /**
    * Create an action that allows for creation of a strand from a partition
    * function save file (PFS). This can be either a single or hybrid strand
    * once created, but it's treated as a single strand because only one file
    * is involved.
    * @param index  the index of the button on the panel
    * @param isRNA  true if the strand is RNA, false if not
    * @return  the SEQ strand creator
    */
   public static StrandCreator createPFSStrand( int index, boolean isRNA ) {
       return new StrandCreator( index, 3, isRNA, false );
   }

   /**
    * Create an action that allows for creation of a strand from a SAV (save)
    * file. This can be either a single or hybrid strand once created but it's
    * treated as a single strand because only one file is involved.
    * @param index  the index of the button on the panel
    * @param isRNA  true if the strand is RNA, false if not
    * @return  the SAV strand creator
    */
   public static StrandCreator createSAVStrand( int index, boolean isRNA ) {
       return new StrandCreator( index, 4, isRNA, false );
   }

   /**
    * Create an action that allows for creation of a strand from a SEQ file.
    * @param index  the index of the button on the panel
    * @param isRNA  true if the strand is RNA, false if not
    * @return  the SEQ strand creator
    */
   public static StrandCreator createSEQStrand( int index, boolean isRNA ) {
       return new StrandCreator( index, 2, isRNA, false );
   }

    /**
     * Create a strand using the appropriate input.
     * If this method creates a hybrid strand it assumes that both strands are
     * of the same type.
     */
    protected void selectInput() {
        int index = (Integer)data[0];
        int type = (Integer)data[1];
        boolean isRNA = (Boolean)data[2];
        boolean multi = (Boolean)data[3];
        TextInputPanel panel =
            (TextInputPanel)(((JButton)event.getSource()).getParent());

        String file = panel.getFields()[index].getText();

        // If a valid file was selected, create the strand
        if( !file.equals( "" ) ) {
            RNA strand;

            // Create a strand, either a single strand or hybrid of two.
            if( !multi ) { strand = new RNA( file, type, isRNA ); }
            else {
                String prev = panel.getFields()[index-1].getText();
                strand = new HybridRNA( prev, type, file, type, isRNA );
            }

            // If the strand was created successfully, enable its menus
            int code = strand.GetErrorCode();
            if( code == 0 ) {
                ((InputWindow)(((JComponent)event.getSource())
		    .getTopLevelAncestor())).getDataHolder()
		    .setStrand( strand );
                RNAstructure.getBar().enableMenus();
            } else {
                RNAstructureInfoDialog.error(
		    strand.GetErrorMessage( code ) );
            }
        }
    }
}
