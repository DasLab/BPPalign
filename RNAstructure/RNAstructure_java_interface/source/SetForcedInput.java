/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that handles the action of inputting a forced constraint.
 * @author Jessica Reuter
 */
public class SetForcedInput extends Modifier {
    /**
     * The constraint being forced
     */
    private String constraint;

    /**
     * Constructor
     * @param constraint  the constraint to be forced
     */
    public SetForcedInput( String constraint ) {
        this.constraint = constraint;
    }

    /**
     * Force a constraint.
     */
    public void modify() {
        RNA strand = RNAstructure.getCurrentStrand();
        int error = 0;

        // Constraints that don't deal with pairing have simple constraints.
        boolean complex = ( constraint.equals( "Force Pair" ) ||
            constraint.equals( "Prohibit Base Pairs" ) );
        boolean dynalign =
            constraint.equals( "Force Alignment" );

        if( !complex && !dynalign ) {
            Integer c = ForceInputWindow.getSimpleConstraint();

            if( constraint.equals( "Chemically Modified" ) ) {
                error = strand.ForceModification( c );
            } else if( constraint.equals( "Force Double" ) ) {
                error = strand.ForceDoubleStranded( c );
            } else if( constraint.equals( "Force Single" ) ) {
                error = strand.ForceSingleStranded( c );
            } else if( constraint.equals( "U in GU Pair" ) ) {
                error = strand.ForceFMNCleavage( c );
            }
        }

        // Constraints that deal with pairing have complex constraints.
        else if( complex ) {
            Integer[] constraints =
                ForceInputWindow.getComplexConstraints();
            int helixLength = constraints[2];

            for( int i = 0; i < helixLength; i++ ) {
                if( constraint.equals( "Force Pair" ) ) {
                    error = strand.ForcePair(
                        constraints[0] + i, constraints[1] - i );
                } else if( constraint.equals( "Prohibit Base Pairs" ) ) {
                    error = strand.ForceProhibitPair(
                        constraints[0] + i, constraints[1] - i );
                }

                if( error != 0 ) { break; }
            }
        }

        // Dynalign alignment constraints use the Dynalign_object.
        else {
            Integer[] aligns = ForceInputWindow.getDynalignConstraints();
            ((DynalignFoldWindow)RNAstructure.getCurrentWindow())
                .addForcedAlignment( aligns[0], aligns[1] );
        }

        if( error != 0 ) {
            RNAstructureInfoDialog.error( strand.GetErrorMessage( error ) );
        } else { ForceInputWindow.reset(); }
    }
}
