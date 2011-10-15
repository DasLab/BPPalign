/*
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * 
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program.
 */

package source;

/**
 * A class that handles the creation of menus on the main RNAstructure menu bar
 * @author Jessica Reuter
 */
public class MenuMaker {
    /**
     * Create the DNA menu.
     * @return  the DNA menu
     */
    public SmartMenu createDNAMenu() { return makeNucleicAcidMenu( "DNA" ); }

    /**
     * Create the draw menu.
     * @return  the draw menu
     */
    public SmartMenu createDrawMenu() {
        SmartMenu menu = new SmartMenu( "Draw" );

        menu.createItem( "Zoom",
            new WindowDisplayer(
                new SpinnerWindow.ZoomFactory() ),
                    "Select the zoom factor for a dot plot." );

        menu.createItem( "Choose Colors",
            new WindowDisplayer(
                new SpinnerWindow.ColorFactory() ),
                    "Select the number of colors in a dot plot." );

        menu.createItem( "Plot Range",
            new WindowDisplayer(
                new PlotRangeWindow.Factory() ),
            "Select the range of values displayed on the plot." );

        return menu;
    }

    /**
     * Create a Dynalign alignment constraints menu.
     * @return  the alignment constraints menu
     */
    public SmartMenu createDynalignAlignmentMenu() {
        SmartMenu menu = makeConstraintMenu( "Constraints for Alignment" );
        menu.setSectionsVisible( 7 );
        return menu;
    }

    /**
     * Create a Dynalign constraints menu.
     * @param sequence  the sequence these constraints will be applied to
     * @return  the constraints menu
     */
    public SmartMenu createDynalignConstraintsMenu( int sequence ) {
        String menuString =
            "Constraints for Sequence " + Integer.toString( sequence );
        SmartMenu menu = makeConstraintMenu( menuString );
        menu.setSectionsVisible( 1, 5, 6 );
        return menu;
    }

    /**
     * Create a menu which can do cut, copy, and paste of text
     * @return  the edit menu
     */
    public SmartMenu createEditMenu() {
        SmartMenu menu = new SmartMenu( "Edit" );

        // Initialize variables
        ActionBase[] listeners = {
            MovementAction.createTextCut(),
            MovementAction.createTextCopy(),
            MovementAction.createTextPaste()
        };

        String[] labels = { "Cut", "Copy", "Paste" };
        String[] strokes = { "control X", "control C", "control V" };

        String[] texts = new String[3];
        for( int i = 0; i < listeners.length; i++ ) {
            String moveType = ( i != 2 ) ? "to" : "from";
            texts[i] = labels[i] + " a block of text " + moveType +
                " the clipboard.";
        }

        // Create menu items
        for( int i = 0; i < listeners.length; i++ ) {
            menu.createItem( labels[i], listeners[i], texts[i], strokes[i] );
        }

        return menu;
    }

    /**
     * Create the file menu, with all possible options.
     * @return  the file menu
     */
    public SmartMenu createFileMenu() {
        SmartMenu menu = new SmartMenu( "File" );

        menu.createItem( "New Sequence",
            new WindowDisplayer(
                new NewSequenceWindow.Factory( "New Sequence" ) ),
                "Create a new sequence.", "control N" );

        menu.createItem( "Open Sequence",
            new SequenceOpener(),
            "Open an existing sequence.", "control O" );

        menu.addSeparator();
        menu.addSectionStart( 3 );

        menu.createItem( "Save Sequence",
            new SequenceSaver( false ),
            "Save a sequence with its existing name." );

        menu.createItem( "Save Sequence As...",
            new SequenceSaver( true ),
            "Save a sequence with a new name." );

        menu.addSeparator();
        menu.addSectionStart( 6 );

        menu.createItem( "OligoScreen", 
            new WindowDisplayer( new OligoScreenWindow.Factory() ),
            "Calculate thermodynamic parameters for a set of oligonucleotides." );

        menu.addSeparator();
        menu.addSectionStart( 8 );

        menu.createItem( "Draw",
            new WindowDisplayer( new Imager.StructureFactory() ),
            "Draw a secondary structure." );

        menu.createItem( "Dot Plot",
            new WindowDisplayer( new Imager.DotPlotFactory(
                Imager.DotPlotFactory.ENERGY_PLOT ) ),
            "Display the energy dot plot for a sequence that was previously folded." );

        menu.createItem( "Dot Plot Partition Function",
            new WindowDisplayer( new Imager.DotPlotFactory(
                Imager.DotPlotFactory.PROBABILITY_PLOT ) ),
            "Display base pairing probabilities for a previously calculated " +
                "sequence." );

        menu.createItem( "Dot Plot Dynalign",
	    new Imager.TwoImagersFactory(),
            "Generate a Dynalign dot plot for two sequences." );

        menu.createItem( "Dot Plot From Text File",
            new WindowDisplayer( new Imager.DotPlotFactory(
                Imager.DotPlotFactory.UNKNOWN_PLOT ) ),
            "Draw a dot plot from a text file." );

        menu.addSeparator();
        menu.addSectionStart( 14 );

        menu.createItem( "Refold From Save File",
            new WindowDisplayer( new RefoldWindow.Factory() ),
            "Refold a sequence from its save file." );

	menu.createItem( "Refold From Dynalign Save File",
	    new WindowDisplayer( new DynalignRefoldWindow.Factory() ),
            "Refold from a Dynalign calculation." );

        menu.addSeparator();
        menu.addSectionStart( 17 );

        menu.createItem( "Print Setup",
            new PrintCurrentResult(),
            "Print the selected content." );

        menu.addSeparator();
        menu.addSectionStart( 19 );

        menu.createItem( "Exit",
            new Exiter(),
            "Exit the application." );

        return menu;
    }

    /**
     * Create the force menu for double stranded folding.
     * @return  the force menu for double stranded folding.
     */
    public SmartMenu createForceDoubleMenu() {
        SmartMenu menu = makeConstraintMenu( "Force" );
        menu.setSectionsVisible( 2 );
        return menu;
    }

    /**
     * Create the force menu for strand refolding.
     * @return  the force menu for strand refolding.
     */
    public SmartMenu createForceRefoldMenu() {
        SmartMenu menu = makeConstraintMenu( "Force" );
        menu.setSectionsVisible( 5 );
        return menu;
    }

    /**
     * Create the menu for forcing constraints on a single strand. The menu can
     * be used both for single strand folding and single strand partitioning.
     * @return  the force menu for a single strand
     */
    public SmartMenu createForceSingleMenu() {
        SmartMenu menu = makeConstraintMenu( "Force" );
        menu.setSectionsVisible( 1, 3, 4, 5, 6 );
        return menu;
    }

    /**
     * Create the force menu for suboptimal structure folding.
     * @return  the force menu for suboptimal structure folding.
     */
    public SmartMenu createForceSuboptimalMenu() {
        SmartMenu menu = makeConstraintMenu( "Force" );
        menu.setSectionsVisible( 1, 3, 4, 5, 6 );
        return menu;
    }

    /**
     * Create the RNAstructure help menu
     * @return  the help menu
     */
    public SmartMenu createHelpMenu() {
        SmartMenu menu = new SmartMenu( "Help" );

        menu.createItem( "Help Topics",
            new WindowDisplayer( new HelpWindow.Factory( true ) ),
            "Get online help.", "F1" );

        menu.addSeparator();

        menu.createItem( "About RNAstructure...",
            new WindowDisplayer( new HelpWindow.Factory( false ) ),
            "Display program information, version number, and copyright." );

        return menu;
    }

    /**
     * Create maximum loop menu.
     * @return  the maximum loop menu.
     */
    public SmartMenu createMaxLoopMenu() {
        SmartMenu menu = new SmartMenu( "Maximum Loop" );

        menu.createItem( "Set Maximum Loop Size",
            new WindowDisplayer( new MaxLoopWindow.Factory() ),
            "Set the maximum number of unpaired nucleotides allowed." );

        return menu;
    }

    /**
     * Create the menu that switches between views of OligoWalk results.
     * @return  the oligo menu
     */
    public SmartMenu createOligoMenu() {
	SmartMenu menu = new SmartMenu( "Oligo Graph" );

	menu.createItem( "Overall and Duplex",
			 OligoGraphSwitcher.overallAndDuplex(),
			 "Show overall and duplex oligo results." );

	menu.addSeparator();

	menu.createItem( "Overall",
			 OligoGraphSwitcher.overall(),
			 "Show overall oligo results." );

	menu.createItem( "Duplex",
			 OligoGraphSwitcher.duplex(),
			 "Show duplex oligo results." );

	menu.createItem( "Broken Target Structure",
			 OligoGraphSwitcher.breakTargetStructure(),
			 "Show broken target structure oligo results." );

	menu.createItem( "Oligomer Unimolecular",
			 OligoGraphSwitcher.unimolecular(),
			 "Show unimolecular oligo results." );

	menu.createItem( "Oligomer Bimolecular",
			 OligoGraphSwitcher.bimolecular(),
			 "Show bimolecular oligo results." );

	return menu;
    }

    /**
     * Create the dot plot output menu for folding save files.
     * @return  the dot plot output menu
     */
    public SmartMenu createPlotOutputFoldingMenu() {
        SmartMenu menu = new SmartMenu( "Output Plot" );

        menu.createItem( "Write to Dot Plot File",
            DotPlotFileWriter.exportToDotPlotFile(),
            "Write dot plot to a file. This format saves only visible dots." );

        menu.createItem( "Write to Postscript Image File",
            DotPlotFileWriter.exportToPostscriptFile(),
            "Write dot plot to a Postscript image, which shows all dots, " +
                "not only those that are visible." );

        return menu;
    }

    /**
     * Create the dot plot output menu for partition function save file.
     * @return  the dot plot output menu
     */
    public SmartMenu createPlotOutputPartitionMenu() {
	SmartMenu menu = createPlotOutputFoldingMenu();

	menu.createItem( "Write Probable Structures",
	    ProbableStructureFileWriter.exportToProbableFile(),
            "Write a ct file containing structures which are composed of " +
	         "probable pairs of different levels." );

	return menu;
    }

    /**
     * Create out loud reading menu.
     * @return  the reading menu.
     */
    public SmartMenu createReadMenu() {
        SmartMenu menu = new SmartMenu( "Read" );

        menu.createCheckItem( "Read While Typing",
            ReadAloud.readWhileTyping(),
            "Read a sequence out loud as it is typed into the keyboard." );

        return menu;
    }

    /**
     * Create the RNA menu.
     * @return  the RNA menu
     */
    public SmartMenu createRNAMenu() {
        SmartMenu menu = makeNucleicAcidMenu( "RNA" );

        menu.createItem( "OligoWalk",
            new WindowDisplayer( new OligoWalkWindow.Factory( "RNA" ) ),
            "RNA" + " OligoWalk calculation." );

        menu.createItem( "Break RNA Pseudoknots",
            new WindowDisplayer( new PseudoknotWindow.Factory() ),
            "Break pseudoknots in a structure, leaving the lowest free energy pseudoknot-free structure." );

        return menu;
    }

    /**
     * Create one of the various menus for sequence constraints.
     * @param num  the number of the sequence
     * @return  a sequence constraints menu
     */
    public SmartMenu createSequenceConstraintsMenu( int num ) {
        String title = "Constraints for Sequence " + Integer.toString( num );
        SmartMenu menu = makeConstraintMenu( title );
        menu.setSectionsVisible( 1, 5, 6 );
        return menu;
    }

    /**
     * Create the structure annotation menu.
     * @return  the structure annotation menu
     */
    public SmartMenu createStructureAnnotationMenu() {
        SmartMenu menu = new SmartMenu( "Annotations" );

        menu.createItem( "Add Probability Annotation",
            AnnotationWriter.probability(),
            "Add probability color annotation to the structure." );

        String key1 = "Probability Annotation Key";
        menu.createItem( "Show " + key1,
            new WindowDisplayer( new AnnotationKeyWindow.Factory( key1 ) ),
            "Show the probability color annotation key," );

        menu.addSeparator();

        menu.createItem( "Add SHAPE Annotation",
            AnnotationWriter.shape(),
            "Add SHAPE color annotation to the structure." );

        String key2 = "SHAPE Annotation Key";
        menu.createItem( "Show " + key2,
            new WindowDisplayer( new AnnotationKeyWindow.Factory( key2 ) ),
            "Show the SHAPE color annotation key," );

        menu.addSeparator();

        menu.createItem( "Remove Annotation",
            AnnotationWriter.remove(),
            "Remove annotation from this structure." );

        return menu;
    }

    /**
     * Create the structure output menu.
     * @return  the structure output menu
     */
    public SmartMenu createStructureMenu() {
        SmartMenu menu = new SmartMenu( "Draw" );

        menu.createItem( "Go to Structure...",
            new WindowDisplayer( new SpinnerWindow.StructureChooserFactory() ),
            "Switch to a selected structure." );

        menu.createItem( "Zoom",
            new WindowDisplayer( new SpinnerWindow.ZoomFactory() ),
            "Zoom this structure." );

        menu.addSeparator();

	menu.createCheckItem( "Render Clockwise/Counterclockwise",
	    FlipStructure.flip(),
	     "Render a structure clockwise (if checked) or counterclockwise " +
                  "(if not checked)" );

	menu.addSeparator();

        menu.createItem( "Write Dot Bracket File",
            BracketStructureFileWriter.exportToDotBracketFile(),
            "Write a dot bracket file." );

        menu.createItem( "Write Helix (Text) File",
            HelixStructureFileWriter.exportToHelixFile(),
            "Write a text file of the helices in this structure." );

        menu.createItem( "Write Structure Group Postscript File(s)",
            new WindowDisplayer( new PostscriptDialogWindow.Factory() ),
            "Write the group of structures to Postscript files, with " +
                "optional annotation if desired." );

        return menu;
    }

    /**
     * Create the temperature menu.
     * @return  the temperature menu
     */
    public SmartMenu createTemperatureMenu() {
        SmartMenu menu = new SmartMenu( "Temperature" );

        menu.createItem( "Set Temperature",
            new WindowDisplayer( new TemperatureWindow.Factory() ),
            "Set the temperature at which calculations occur." );

        return menu;
    }

    /**
     * Make a menu for one of the various constraint menus
     * @param title  the title of the menu
     */
    private SmartMenu makeConstraintMenu( String title ) {
        SmartMenu menu = new SmartMenu( title );

        menu.createItem( "Base Pair",
            new WindowDisplayer(
                new ForceInputWindow.Factory( "Force Pair" ) ),
            "Force two bases to be paired in a folded structure." );

        menu.createItem( "Prohibit Base Pairs",
            new WindowDisplayer(
                new ForceInputWindow.Factory( "Prohibit Base Pairs" ) ),
            "Mandate that certain bases be unable to pair in a structure." );

        menu.createItem( "Chemical Modification",
            new WindowDisplayer(
                new ForceInputWindow.Factory( "Chemically Modified" ) ),
            "Modify specified bases chemically prior to folding." );

        menu.createItem( "FMN Cleavage",
            new WindowDisplayer(
                new ForceInputWindow.Factory( "U in GU Pair") ),
            "Cleave a strand of interest." );

        menu.createItem( "Single Stranded",
            new WindowDisplayer(
                new ForceInputWindow.Factory( "Force Single" ) ),
            "Force particular regions to be single stranded." );

        menu.createItem( "Double Stranded",
            new WindowDisplayer(
                new ForceInputWindow.Factory( "Force Double" ) ),
            "Force particular regions to be double stranded." );

        menu.addSeparator();
        menu.addSectionStart( 7 );

        menu.createCheckItem( "Forbid Unimolecular Pairs",
            new PermitMolecularity(),
            "Forbid pairs from forming between bases of the same molecule." );

        menu.addSeparator();
        menu.addSectionStart( 9 );

        menu.createItem( "Maximum Pairing Distance",
            new WindowDisplayer( new MaxPairWindow.Factory() ),
            "Limit the maximum distance allowed between paired bases." );

        menu.addSeparator();
        menu.addSectionStart( 11 );

        menu.createItem( "Read SHAPE Reactivity -- Hard Constraints",
            new WindowDisplayer( new SHAPEHardWindow.Factory() ),
            "Input constraints on SHAPE data to be used in folding." );

        menu.createItem( "Read SHAPE Reactivity -- Pseudo-Energy Constraints",
            new WindowDisplayer( new SHAPEEnergyWindow.Factory() ),
            "Input constraints on SHAPE data to be used in folding." );

        menu.addSeparator();
        menu.addSectionStart( 14 );

        menu.createItem( "Show Current Constraints",
            new WindowDisplayer( new ConstraintsIOHandler.Factory( 4 ) ),
            "View constraints that are currently applied to this sequence." );

        menu.createItem( "Reset Current Constraints",
            new WindowDisplayer( new ConstraintsIOHandler.Factory( 3 ) ),
            "Reset the constraints currently applied to this sequence." );

        menu.addSeparator();
        menu.addSectionStart( 17 );

        menu.createItem( "Save Constraints",
            new WindowDisplayer( new ConstraintsIOHandler.Factory( 2 ) ),
            "Save folding constraints to a save file." );

        menu.createItem( "Restore Constraints",
            new WindowDisplayer( new ConstraintsIOHandler.Factory( 1 ) ),
            "Get a set of folding constraints from a save file." );

        menu.addSeparator();
        menu.addSectionStart( 20 );

        menu.createItem( "Force Alignment",
            new WindowDisplayer(
                new ForceInputWindow.Factory( "Force Alignment" ) ),
            "Force the alignment of particular bases." );

        menu.addSeparator();

        menu.createItem( "Show Current Alignment Constraints",
            new WindowDisplayer( new DynalignAlignmentIOHandler.Factory( 4 ) ),
            "Show constraints set for this particular alignment." );

        menu.createItem( "Reset Alignment Constraints",
            new WindowDisplayer( new DynalignAlignmentIOHandler.Factory( 3 ) ),
            "Reset the constraints that are applied to this alignment." );

        menu.addSeparator();

        menu.createItem( "Save Alignment",
            new WindowDisplayer( new DynalignAlignmentIOHandler.Factory( 2 ) ),
            "Save the alignment to a save file." );

        menu.createItem( "Restore Alignment From File",
            new WindowDisplayer( new DynalignAlignmentIOHandler.Factory( 1 ) ),
            "Get an alignment from a save file." );

        return menu;
    }

    /**
     * Make a menu for one of the two nucleic acid menus
     * @param acid  the nucleic acid type
     */
    private SmartMenu makeNucleicAcidMenu( String acid ) {
        SmartMenu menu = new SmartMenu( acid );

        menu.createItem( "Fold " + acid + " Single Strand",
            new WindowDisplayer( new FoldSingleWindow.Factory( acid ) ),
            "Predict a secondary structure for a single strand using " + acid +
                " thermodynamics." );

        menu.createItem( "Fold " + acid + " Bimolecular",
            new WindowDisplayer( new FoldDoubleWindow.Factory( acid ) ),
            "Predict a secondary structure for two strands using " + acid +
                " thermodynamics." );

        menu.addSeparator();

        menu.createItem( "Partition Function " + acid,
            new WindowDisplayer( new PartitionSingleWindow.Factory( acid ) ),
            "Predict base pairing probabilities for all " + acid + " pairs." );

        menu.createItem( "Partition Function " + acid + " Bimolecular",
            new WindowDisplayer( new PartitionDoubleWindow.Factory( acid ) ),
            "Predict base pairing probabilities for all " + acid + "pairs" +
                " for two strands." );

        menu.addSeparator();

        menu.createItem( "Generate All Suboptimal " + acid + " Structures",
            new WindowDisplayer( new FoldSuboptimalWindow.Factory( acid ) ),
            "Generate all structures within a given increment of " +
                "the lowest free energy structure." );

        menu.createItem( "Stochastic " + acid + " Sampling",
            new WindowDisplayer( new StochasticWindow.Factory( acid ) ),
            "Predict " + acid + " structures using stochastic sampling." );

	menu.createItem(
            "MaxExpect: Predict " + acid + " MEA Structure",
	    new WindowDisplayer( new MaxExpectWindow.Factory( acid ) ),
            "Predict the maximum expected accuracy" + acid +
                " secondary structure." );

	menu.createItem( "ProbKnot: Predict " + acid + " Structures Including " +     	"Pseudoknots",
             new WindowDisplayer( new ProbKnotWindow.Factory() ),
             "Predict pseudoknots in a " + acid + " structure." );

        menu.createItem( "Efn2 " + acid,
            new WindowDisplayer( new Efn2Window.Factory( acid ) ),
            "Calculate the free energy of a(n) " + acid + " structure." );

        menu.addSeparator();

        menu.createItem( acid + " Dynalign",
            new WindowDisplayer( new DynalignFoldWindow.Factory( acid ) ),
            "Find a common secondary structure for two " + acid + " sequences." );

        return menu;
    }
}
