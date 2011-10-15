/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * A window that deals with either reading or writing folding constraints so
 * the RNAstructure GUI can either input or output them.
 * @author Jessica Reuter
 */
public class ConstraintsIOHandler {
    /**
     * A serialized ID, required for extending GUI classes
     */
    private static final long serialVersionUID = -335674004754176548L;

    /**
     * Constructor
     * @param type  the type of constraints action to be taken
     *              1.  Read constraints from a file
     *              2.  Write constraints to a file
     *              3.  Remove all constraints
     *              4.  Show all constraints
     */
    public ConstraintsIOHandler( int type ) {
        RNAstructure.setLabel( "For help, press F1." );

        if( type == 1 ) { runIO( 1 ); }
        else if( type == 2 ) { runIO( 2 ); }
        else if( type == 3 ) { resetConstraints(); }
        else if( type == 4 ) { showConstraints(); }
        else throw new IllegalArgumentException( "Undefined window type" );
    }

    /**
     * Reset all folding constraints.
     */
    private void resetConstraints() {
        RNAstructureInfoDialog.confirm(
            "This will erase all constraints.\nContinue?",
            new ChainedAction(
                new StrandResetter(), new CloseAction() ) );
    }

    /**
     * Either read or write folding constraints from or to a file.
     * @param type  1=read, 2=write
     */
    public void runIO( int type ) {
        String title = "Open";
        if( type == 2 ) { title = "Save"; }

        String file = new FileSelector( true, title, "Constraint Files,con" )
            .showSelector();

        if( !file.equals( "" ) ) {
            RNA strand = RNAstructure.getCurrentStrand();
            if( type == 1 ) { strand.ReadConstraints( file ); }
            else { strand.WriteConstraints( file ); }
        }
    }

    /**
     * Show all folding constraints
     */
    private void showConstraints() {
        // Set up variables and data structures
        StringBuilder constraints = new StringBuilder( "<html>" );
        RNA strand = RNAstructure.getCurrentStrand();

        LinkedHashMap<String, Integer> typeMap =
            new LinkedHashMap<String, Integer>();
        TreeMap<Integer, Integer> map = new TreeMap<Integer, Integer>();
        TreeSet<Integer> set = new TreeSet<Integer>();

        typeMap.put( "Forced Base Pairs", strand.GetNumberOfForcedPairs() );
        typeMap.put( "Forced Prohibited Pairs",
            strand.GetNumberOfForcedProhibitedPairs() );
        typeMap.put( "Forced Modifications", 
            strand.GetNumberOfForcedModifications() );
        typeMap.put( "Forced FMN Cleavages",
            strand.GetNumberOfForcedFMNCleavages() );
        typeMap.put( "Forced Single Stranded",
            strand.GetNumberOfForcedSingleStranded() );
        typeMap.put( "Forced Double Stranded", 
            strand.GetNumberOfForcedDoubleStranded() );

        // Go through all constraints and list them
        Iterator<Map.Entry<String, Integer>> iterator =
            typeMap.entrySet().iterator();

        int counter = 1;
        while( iterator.hasNext() ) {
            Map.Entry<String, Integer> current = iterator.next();
            constraints.append( current.getKey() ).append( "<br/>" );
            Integer value = current.getValue();
            if( value == 0 ) {
                constraints.append( "None" );
                if( iterator.hasNext() ) {
                    constraints.append( "<br/><br/>" );
                }
            } else {
                for( int i = 0; i < value; i++ ) {
                    switch( counter ) {
                    case 1:
                        map.put( strand.GetForcedPair( i, true ), 
                            strand.GetForcedPair( i, false ) );
                        break;
                    case 2:
                        map.put( strand.GetForcedProhibitedPair( i, true ), 
                            strand.GetForcedProhibitedPair( i, false ) );
                        break;
                    case 3:
                        set.add( strand.GetForcedModification( i ) );
                        break;
                    case 4:
                        set.add( strand.GetForcedFMNCleavage( i ) );
                        break;
                    case 5:
                        set.add( strand.GetForcedSingleStranded( i ) );
                        break;
                    case 6:
                        set.add( strand.GetForcedDoubleStranded( i ) );
                        break;
                    }
                }

                int counter2 = 1;
                if( map.size() > 0 ) {
                    Iterator<Map.Entry<Integer, Integer>> iterator2 =
                        map.entrySet().iterator();

                    while( iterator2.hasNext() ) {
                        Map.Entry<Integer, Integer> current2 =
                            iterator2.next();

                        constraints.append( current2.getKey() );
                        constraints.append( "-" );
                        constraints.append( current2.getValue() );

                        if( iterator2.hasNext() ) {
                            constraints.append( ", " );
                            if( counter2 % 5 == 0 ) {
                                constraints.append( "<br/>" );
                            }
                        } else if( iterator.hasNext() ) {
                            constraints.append( "<br/><br/>" );
                        }
                        counter2++;
                    }
                } else if( set.size() > 0 ) {
                    for( Integer element: set ) {
                        constraints.append( element );

                        if( counter2 != set.size() ) {
                            constraints.append( ", " );
                            if( counter2 % 5 == 0 ) {
                                constraints.append( "<br/>" );
                            }
                        } else if( iterator.hasNext() ) {
                            constraints.append( "<br/><br/>" );
                        }
                        counter2++;
                    }
                }

                map.clear();
                set.clear();
            }
            counter++;
        }

        // Show final list dialog
        RNAstructureInfoDialog.message( constraints.toString() );
    }

    /**
     * An inner class which allows a ConstraintsIOHandler to be constructed
     * from the more generic context of the RNAWindowFactory.
     * @author Jessica Reuter
     */
    public static class Factory implements RNAWindowFactory {
        /**
         * The type of ConstraintsIOHandler to create.
         */
        private int type;

        /**
         * Constructor
         * @param type  the type of window to create
         */
        public Factory( int type ) { this.type = type; }

        /**
         * Create a ConstraintsIOHandler.
         * @return  a new ConstraintsIOHandler
         */
        public ConstraintsIOHandler createWindow() {
            return new ConstraintsIOHandler( type );
        }
    }
}