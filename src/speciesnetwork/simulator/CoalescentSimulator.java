package speciesnetwork.simulator;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import beast.app.seqgen.*;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Runnable;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.*;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import speciesnetwork.EmbeddedTree;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;
import speciesnetwork.SanityChecks;

/**
 * @author Chi Zhang
 */

@Description("Simulate gene trees given a species network (multispecies coalescent).")
public class CoalescentSimulator extends Runnable {
    public final Input<State> startStateInput =
            new Input<>("state", "elements of the state space");

    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "Species network for embedding the gene trees.", Validate.REQUIRED);
    public final Input<BirthHybridSimulator> networkSimulatorInput =
            new Input<>("networkSimulator", "Species network simulator.", Validate.XOR, speciesNetworkInput);
    public final Input<RealParameter> popSizesInput =
            new Input<>("popSizes", "Constant per-branch population sizes.", Validate.REQUIRED);
    public final Input<TaxonSet> taxonSuperSetInput =
            new Input<>("taxonset", "Super-set of taxon sets mapping lineages to species.", Validate.REQUIRED);

    public final Input<List<EmbeddedTree>> geneTreesInput =
            new Input<>("geneTree", "Gene tree embedded in the species network.", new ArrayList<>());
    public final Input<RealParameter> ploidiesInput =
            new Input<>("ploidy", "Ploidy (copy number) for each gene (default is 2).");
    public final Input<List<SequenceSimulator>> seqSimulatorsInput =
            new Input<>("sequenceSimulator", "Sequence simulator.", new ArrayList<>());

    public final Input<String> outputFileNameInput =
            new Input<>("outputFileName", "If provided, write to this file rather than to standard out.");
    public final Input<Integer> iterationsInput =
            new Input<>("iterations", "Number of iterations to simulate (default is 1).");
    public final Input<Boolean> networkOperatorInput =
            new Input<>("networkOperator", "Whether to write network topology operators (default false).", false);
    public final Input<String> initMethodInput =
            new Input<>("initMethod", "Initializing method (point, random, user).", "random");

    private Network speciesNetwork;
    private RealParameter popSizes;
    private List<EmbeddedTree> geneTrees;
    private RealParameter ploidies;

    private int nrOfGeneTrees;
    private Multimap<NetworkNode, Node> networkNodeGeneLineagesMap = HashMultimap.create();
    private int nodeIndex;  // gene tree internal node index

    private List<SequenceSimulator> seqSimulators;
    private List<Alignment> alignments = new ArrayList<>();

    @Override
    public void initAndValidate() {
        geneTrees = geneTreesInput.get();
        // sanity check
        if (geneTrees == null)
            throw new RuntimeException("Check gene tree input!");
        nrOfGeneTrees = geneTrees.size();

        if ((ploidies = ploidiesInput.get()) == null)
            ploidies = new RealParameter("2.0");  // default
        ploidies.setDimension(nrOfGeneTrees);

        seqSimulators = seqSimulatorsInput.get();
    }

    @Override
    public void run() throws IOException {
        // initialize state nodes, essential
        State state = startStateInput.get();
        if (state == null)
            throw new RuntimeException("Input 'state' must be specified!");
        state.initialise();

        final int nrOfIterations;
        if (iterationsInput.get() == null)
            nrOfIterations = 1;
        else
            nrOfIterations = iterationsInput.get();
        for (int iteration = 0; iteration < nrOfIterations; iteration++) {
            simulate();

            String outputFileName = outputFileNameInput.get();
            if (nrOfIterations == 1) {
                writeXMLOutput(outputFileName);  // generate an XML file for a single iteration
            } else {
                writeGeneTrees(outputFileName);  // otherwise, only output the gene trees
            }
        }
    }

    public void simulate() {
        if (speciesNetworkInput.get() == null)
            speciesNetwork = networkSimulatorInput.get().simulate();  // simulate a species network
        else
            speciesNetwork = speciesNetworkInput.get();
        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());

        // set popSizes dimension
        popSizes = popSizesInput.get();
        final int speciesBranchCount = speciesNetwork.getBranchCount();
        popSizes.setDimension(speciesBranchCount);

        final int traversalNodeCount = speciesNetwork.getTraversalNodeCount();
        // simulate each gene tree and alignment
        for (int ig = 0; ig < nrOfGeneTrees; ig++) {
            EmbeddedTree geneTree = geneTrees.get(ig);

            // initialize embedding matrix to -1 (no traversal)
            geneTree.resetEmbedding(traversalNodeCount, -1);

            networkNodeGeneLineagesMap.clear();
            // generate map of tip names to tip nodes
            final Map<String, NetworkNode> speciesNodeMap = new HashMap<>();
            for (NetworkNode leafNode : speciesNetwork.getLeafNodes()) {
                final String speciesName = leafNode.getLabel();
                speciesNodeMap.put(speciesName, leafNode);
            }
            final Map<String, Node> geneNodeMap = new HashMap<>();
            for (Node leafNode : geneTree.getExternalNodes()) {
                final String geneName = leafNode.getID();
                geneNodeMap.put(geneName, leafNode);
            }
            // multimap of species network tip node to gene tree tip nodes
            final TaxonSet taxonSuperSet = taxonSuperSetInput.get();
            for (Taxon speciesTip : taxonSuperSet.taxonsetInput.get()) {
                final NetworkNode speciesNode = speciesNodeMap.get(speciesTip.getID());
                final TaxonSet speciesTaxonSet = (TaxonSet) speciesTip;
                for (Taxon geneTip : speciesTaxonSet.taxonsetInput.get()) {
                    final Node geneNode = geneNodeMap.get(geneTip.getID());
                    networkNodeGeneLineagesMap.put(speciesNode, geneNode);
                }
            }

            // reset visited indicator
            speciesNetwork.resetAllVisited();
            // simulate the gene tree
            nodeIndex = 0;
            simulateGeneTree(speciesNetwork.getRoot(), geneTree, ploidies.getValue(ig));

            // simulate alignment on the gene tree
            if (seqSimulators.size() > ig) {
                alignments.add(seqSimulators.get(ig).simulate());
            }
        }
    }

    // recursively simulate lineages coalescent in each population
    private void simulateGeneTree(NetworkNode snNode, EmbeddedTree geneTree, double ploidy) {
        if (snNode.isVisited())
            return;
        for (NetworkNode c: snNode.getChildren()) {
            simulateGeneTree(c, geneTree, ploidy);
        }

        snNode.setVisited(true);  // set visited indicator

        final Collection<Node> lineagesAtBottom = networkNodeGeneLineagesMap.get(snNode);

        if (snNode.isReticulation()) {
            // assign lineages at the bottom to the left and right populations
            final Collection<Node> lineagesAtLBottom = new HashSet<>();
            final Collection<Node> lineagesAtRBottom = new HashSet<>();
            for (Node lineage : lineagesAtBottom) {
                if (Randomizer.nextDouble() < snNode.getGammaProb())
                    lineagesAtLBottom.add(lineage);
                else
                    lineagesAtRBottom.add(lineage);
            }

            final double bottomHeight = snNode.getHeight();
            final Integer lBranchNumber = snNode.gammaBranchNumber;
            NetworkNode lParent = snNode.getParentByBranch(lBranchNumber);
            final double lPopSize = popSizes.getValue(lBranchNumber);
            final double lTopHeight = lParent.getHeight();
            List<Node> lineagesAtLTop =
                    simulateCoalescentEvents(lineagesAtLBottom, bottomHeight, lTopHeight, ploidy*lPopSize, geneTree);
            final Integer rBranchNumber = snNode.gammaBranchNumber + 1;
            NetworkNode rParent = snNode.getParentByBranch(rBranchNumber);
            final double rPopSize = popSizes.getValue(rBranchNumber);
            final double rTopHeight = rParent.getHeight();
            List<Node> lineagesAtRTop =
                    simulateCoalescentEvents(lineagesAtRBottom, bottomHeight, rTopHeight, ploidy*rPopSize, geneTree);

            networkNodeGeneLineagesMap.putAll(lParent, lineagesAtLTop);
            networkNodeGeneLineagesMap.putAll(rParent, lineagesAtRTop);
            // update embedding
            final int traversalLParentNr = lParent.getTraversalNumber();
            for (final Node geneNode : lineagesAtLTop)
            	geneTree.setEmbedding(geneNode.getNr(), traversalLParentNr, lBranchNumber);
            final int traversalRParentNr = rParent.getTraversalNumber();
            for (final Node geneNode : lineagesAtRTop)
            	geneTree.setEmbedding(geneNode.getNr(), traversalRParentNr, rBranchNumber);
        }
        else {
            final double bottomHeight = snNode.getHeight();
            final Integer sBranchNumber = snNode.gammaBranchNumber;
            NetworkNode sParent = snNode.getParentByBranch(sBranchNumber);
            final double popSize = popSizes.getValue(sBranchNumber);
            final double topHeight;
            if (sParent.isOrigin())  // network root
                topHeight = Double.POSITIVE_INFINITY;
            else
                topHeight = sParent.getHeight();

            List<Node> lineagesAtTop =
                    simulateCoalescentEvents(lineagesAtBottom, bottomHeight, topHeight, ploidy*popSize, geneTree);
            if (sParent.isOrigin()) {
                geneTree.setRoot(lineagesAtTop.get(0));
            } else {
                networkNodeGeneLineagesMap.putAll(sParent, lineagesAtTop);
                // update embedding
                final int traversalParentNr = sParent.getTraversalNumber();
                for (final Node geneNode : lineagesAtTop)
                	geneTree.setEmbedding(geneNode.getNr(), traversalParentNr, sBranchNumber);
            }
        }
    }

    private List<Node> simulateCoalescentEvents(Collection<Node> lineages, double bottomHeight,
                                                double topHeight, double pNu, EmbeddedTree geneTree) {
        // start from the lineages at the bottom
        List<Node> currentLineages = new ArrayList<>(lineages);
        double currentHeight = bottomHeight;

        List<Node> internalNodes = geneTree.getInternalNodes();

        // then go up backward in time
        while (currentLineages.size() > 1 && currentHeight < topHeight) {
            // generate a coalescent waiting time
            final int nLineage = currentLineages.size();
            final double coalescentRate = nLineage * (nLineage - 1) / (2 * pNu);
            final double waitingTime = Randomizer.nextExponential(coalescentRate);
            currentHeight += waitingTime;

            if (currentHeight < topHeight) {
                // randomly pick two lineages to coalescence
                int rnd = Randomizer.nextInt(nLineage);
                final Node left = currentLineages.get(rnd);
                currentLineages.remove(left);
                rnd = Randomizer.nextInt(nLineage - 1);
                final Node right = currentLineages.get(rnd);
                currentLineages.remove(right);
                // deal with the parent of the two picked nodes
                final Node node = internalNodes.get(nodeIndex++);
                left.setParent(node); right.setParent(node);
                node.setLeft(left);   node.setRight(right);
                node.setHeight(currentHeight);
                currentLineages.add(node);
            }
        }

        return currentLineages;
    }

    private void writeXMLOutput(String outputFileName) throws IOException {
        PrintStream out;  // where to print
        if (outputFileName == null) {
            out = System.out;
        } else {
            String msg = "Writing";
            if (new File(outputFileName).exists())
                msg = "Warning: Overwriting";
            System.err.println(msg + " file " + outputFileName);
            out = new PrintStream(outputFileName);
        }

        DecimalFormat df = new DecimalFormat("#.########");
        // print header
        out.println("<?xml version='1.0' encoding='UTF-8'?>");
        out.println("<beast namespace=\"beast.core:beast.core.util:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.evolution.operators:" +
                    "beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.branchratemodel:beast.evolution.likelihood\" version=\"2.4\">");
        // print sequence data
        for (int i = 0; i < nrOfGeneTrees; i++) {
            out.println("    <data id=\"gene" + (i+1) + "\" name=\"alignment\">");
            if (seqSimulators.size() > i) {  // have simulated alignments
                Alignment alignment = alignments.get(i);
                List<Sequence> sequences = alignment.sequenceInput.get();
                for (Sequence seq : sequences)
                    out.println("        <sequence taxon=\"" + seq.getTaxon() + "\" value=\"" + seq.getData() + "\"/>");
            } else {
            	EmbeddedTree geneTree = geneTrees.get(i);
                for (Node leaf : geneTree.getExternalNodes())
                    out.println("        <sequence taxon=\"" + leaf.getID() + "\" totalcount=\"4\" value=\"-\"/>");
            }
            out.println("    </data>");
        }
        out.println("");  // mappings
        out.println("    <map name=\"Uniform\">beast.math.distributions.Uniform</map>\n" +
                    "    <map name=\"Exponential\">beast.math.distributions.Exponential</map>\n" +
                    "    <map name=\"LogNormal\">beast.math.distributions.LogNormalDistributionModel</map>\n" +
                    "    <map name=\"Normal\">beast.math.distributions.Normal</map>\n" +
                    "    <map name=\"Beta\">beast.math.distributions.Beta</map>\n" +
                    "    <map name=\"Gamma\">beast.math.distributions.Gamma</map>\n" +
                    "    <map name=\"LaplaceDistribution\">beast.math.distributions.LaplaceDistribution</map>\n" +
                    "    <map name=\"InverseGamma\">beast.math.distributions.InverseGamma</map>\n" +
                    "    <map name=\"OneOnX\">beast.math.distributions.OneOnX</map>\n" +
                    "    <map name=\"prior\">beast.math.distributions.Prior</map>\n");
        // print initial species network
        out.println("    <init spec=\"beast.util.TreeParser\" id=\"newick:species\" IsLabelledNewick=\"true\" adjustTipHeights=\"false\"\n" +
                    "          newick=\"" + speciesNetwork.getOrigin().toString(df, true) + "\"/>");
        out.println("    <run id=\"mcmc\" spec=\"MCMC\" chainLength=\"50000000\" storeEvery=\"5000\">");  // MCMC block
        out.println("        <state id=\"state\">");  // states
        // print state nodes
        out.println("            <stateNode id=\"network:species\" spec=\"speciesnetwork.NetworkParser\" tree=\"@newick:species\">");
        out.println("                <taxonset id=\"taxonsuperset\" spec=\"TaxonSet\">");
        final TaxonSet taxonsuperset = taxonSuperSetInput.get();
        for (Taxon speciesTip : taxonsuperset.taxonsetInput.get()) {
            out.println("                    <taxon id=\"" + speciesTip.getID() + "\" spec=\"TaxonSet\">");
            final TaxonSet speciesTaxonSet = (TaxonSet) speciesTip;
            for (Taxon geneTip : speciesTaxonSet.taxonsetInput.get())
                out.println("                        <taxon id=\"" + geneTip.getID() + "\" spec=\"Taxon\"/>");
            out.println("                    </taxon>");
        }
        out.println("                </taxonset>");
        out.println("            </stateNode>");
        out.println("            <parameter id=\"originTime:species\" lower=\"0.0\" name=\"stateNode\">" + df.format(speciesNetwork.getOrigin().getHeight()) + "</parameter>");
        out.println("            <parameter id=\"netDivRate:species\" lower=\"0.0\" name=\"stateNode\">1.0</parameter>");
        out.println("            <parameter id=\"turnOverRate:species\" lower=\"0.0\" upper=\"1.0\" name=\"stateNode\">0.5</parameter>");
        for (int i = 0; i < nrOfGeneTrees; i++) {
            out.println("            <stateNode id=\"tree:gene" + (i+1) + "\" spec=\"speciesnetwork.EmbeddedTree\">");
            out.println("                <taxonset id=\"taxonset:gene" + (i + 1) + "\" spec=\"TaxonSet\" alignment=\"@gene" + (i + 1) + "\"/>");
            out.println("            </stateNode>");
        }
        out.println("        </state>\n");  // end of states
        // print initial/true gene trees
        out.println("        <!--");
        for (int i = 0; i < nrOfGeneTrees; i++) {
            EmbeddedTree geneTree = geneTrees.get(i);
            out.println("        <init spec=\"beast.util.TreeParser\" id=\"newick:gene" + (i+1) + "\" initial=\"@tree:gene" + (i+1) + "\" taxa=\"@gene" + (i+1) + "\" IsLabelledNewick=\"true\"\n" +
                        "              newick=\"" + geneTree.getRoot().toNewick() + "\"/>");
        }
        out.println("        -->");
        // starbeast initializer
        final String initMethod = initMethodInput.get();
        out.println("        <init id=\"SNI\" spec=\"speciesnetwork.SpeciesNetworkInitializer\" estimate=\"false\" method=\"" + initMethod + "\" speciesNetwork=\"@network:species\" origin=\"@originTime:species\">");
        for (int i = 0; i < nrOfGeneTrees; i++)
            out.println("            <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
        out.println("            <rebuildEmbedding id=\"initEmbed\" spec=\"speciesnetwork.operators.RebuildEmbedding\" speciesNetwork=\"@network:species\" taxonset=\"@taxonsuperset\" weight=\"0.0\">");
        for (int i = 0; i < nrOfGeneTrees; i++)
            out.println("                <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
        out.println("            </rebuildEmbedding>");
        out.println("            <coalescentSimulator id=\"coalSim\" spec=\"speciesnetwork.simulator.CoalescentSimulator\" speciesNetwork=\"@network:species\" taxonset=\"@taxonsuperset\">");
        out.println("                <parameter id=\"popSizes\" estimate=\"false\" name=\"popSizes\">" + popSizes.getValue() + "</parameter>");
        for (int i = 0; i < nrOfGeneTrees; i++)
            out.println("                <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
        out.println("            </coalescentSimulator>");
        out.println("        </init>\n");
        // print posterior, prior, and likelihood stuff
        out.println("        <distribution id=\"posterior\" spec=\"util.CompoundDistribution\">");
        // coalescent
        out.println("            <distribution id=\"coalescent\" spec=\"speciesnetwork.MultispeciesCoalescent\" speciesNetwork=\"@network:species\">");
        for (int i = 0; i < nrOfGeneTrees; i++) {
            out.println("                <geneTreeWithin id=\"geneTree:gene" + (i+1) + "\" spec=\"speciesnetwork.GeneTreeInSpeciesNetwork\" speciesNetwork=\"@network:species\" geneTree=\"@tree:gene" + (i+1) + "\"/>");
        }
        out.println("                <!-- populationModel id=\"popModel\" popSizes=\"@popSizes\" spec=\"speciesnetwork.ConstantPopulation\"/ -->");
        out.println("                <populationModel id=\"popModel\" alpha=\"10.0\" mean=\"" + popSizes.getValue() + "\" spec=\"speciesnetwork.ConstantPopIntegrated\"/>");
        out.println("            </distribution>");
        out.println("            <distribution id=\"prior\" spec=\"util.CompoundDistribution\">");  // prior
        // network prior
        out.println("                <distribution id=\"networkPrior\" spec=\"speciesnetwork.BirthHybridizationModel\" network=\"@network:species\" netDiversification=\"@netDivRate:species\" turnOver=\"@turnOverRate:species\" betaShape=\"1.0\"/>");
        out.println("                <prior id=\"networkOrigin\" name=\"distribution\" x=\"@originTime:species\">");
        out.println("                    <Exponential id=\"exponential.0\" name=\"distr\" mean=\"" + df.format(speciesNetwork.getOrigin().getHeight()) + "\"/>");
        out.println("                </prior>");
        out.println("                <prior id=\"netDivPrior\" name=\"distribution\" x=\"@netDivRate:species\">");
        out.println("                    <Exponential id=\"exponential.01\" name=\"distr\" mean=\"10.0\"/>");
        out.println("                </prior>");
        out.println("                <prior id=\"turnOverPrior\" name=\"distribution\" x=\"@turnOverRate:species\">");
        out.println("                    <Beta id=\"betadistr.01\" name=\"distr\" alpha=\"1.0\" beta=\"1.0\"/>");
        out.println("                </prior>");
        out.println("            </distribution>");
        // likelihood
        out.println("            <distribution id=\"likelihood\" spec=\"util.CompoundDistribution\">");  // likelihood
        for (int i = 0; i < nrOfGeneTrees; i++) {
            out.println("                <distribution id=\"likelihood:gene" + (i + 1) + "\" data=\"@gene" + (i+1) + "\" tree=\"@tree:gene" + (i+1) + "\" spec=\"TreeLikelihood\">");
            out.println("                    <siteModel id=\"siteModel:gene" + (i+1) + "\" mutationRate=\"1.0\" spec=\"SiteModel\">");
            out.println("                        <substModel id=\"jc:gene" + (i+1) + "\" spec=\"JukesCantor\"/>");
            out.println("                    </siteModel>");
            out.println("                    <branchRateModel id=\"strictClock:gene" + (i+1) + "\" clock.rate=\"1.0\" spec=\"StrictClockModel\"/>");
            out.println("                </distribution>");
        }
        out.println("            </distribution>");
        out.println("        </distribution>\n");
        // print operators
        // gene tree operators
        for (int i = 0; i < nrOfGeneTrees; i++) {
            out.println("        <operator id=\"scaleAndEmbed:gene" + (i+1) + "\" spec=\"speciesnetwork.operators.RebuildEmbedding\" speciesNetwork=\"@network:species\" taxonset=\"@taxonsuperset\" weight=\"3.0\">");
            out.println("            <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
            out.println("            <operator id=\"scale:gene" + (i+1) + "\" spec=\"ScaleOperator\" scaleFactor=\"0.5\" tree=\"@tree:gene" + (i+1) + "\" weight=\"0.0\"/>");
            out.println("        </operator>");
            out.println("        <operator id=\"scaleRootAndEmbed:gene" + (i+1) + "\" spec=\"speciesnetwork.operators.RebuildEmbedding\" speciesNetwork=\"@network:species\" taxonset=\"@taxonsuperset\" weight=\"3.0\">");
            out.println("            <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
            out.println("            <operator id=\"scaleRoot:gene" + (i+1) + "\" spec=\"ScaleOperator\" rootOnly=\"true\" scaleFactor=\"0.5\" tree=\"@tree:gene" + (i+1) + "\" weight=\"0.0\"/>");
            out.println("        </operator>");
            out.println("        <operator id=\"uniformAndEmbed:gene" + (i+1) + "\" spec=\"speciesnetwork.operators.RebuildEmbedding\" speciesNetwork=\"@network:species\" taxonset=\"@taxonsuperset\" weight=\"30.0\">");
            out.println("            <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
            out.println("            <operator id=\"uniform:gene" + (i+1) + "\" spec=\"Uniform\" tree=\"@tree:gene" + (i+1) + "\" weight=\"0.0\"/>");
            out.println("        </operator>");
            out.println("        <operator id=\"subSlideAndEmbed:gene" + (i+1) + "\" spec=\"speciesnetwork.operators.RebuildEmbedding\" speciesNetwork=\"@network:species\" taxonset=\"@taxonsuperset\" weight=\"15.0\">");
            out.println("            <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
            out.println("            <operator id=\"subSlide:gene" + (i+1) + "\" spec=\"SubtreeSlide\" tree=\"@tree:gene" + (i+1) + "\" weight=\"0.0\"/>");
            out.println("        </operator>");
            out.println("        <operator id=\"narrowAndEmbed:gene" + (i+1) + "\" spec=\"speciesnetwork.operators.RebuildEmbedding\" speciesNetwork=\"@network:species\" taxonset=\"@taxonsuperset\" weight=\"15.0\">");
            out.println("            <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
            out.println("            <operator id=\"narrow:gene" + (i+1) + "\" spec=\"Exchange\" tree=\"@tree:gene" + (i+1) + "\" weight=\"0.0\"/>");
            out.println("        </operator>");
            out.println("        <operator id=\"wideAndEmbed:gene" + (i+1) + "\" spec=\"speciesnetwork.operators.RebuildEmbedding\" speciesNetwork=\"@network:species\" taxonset=\"@taxonsuperset\" weight=\"5.0\">");
            out.println("            <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
            out.println("            <operator id=\"wide:gene" + (i+1) + "\" spec=\"Exchange\" isNarrow=\"false\" tree=\"@tree:gene" + (i+1) + "\" weight=\"0.0\"/>");
            out.println("        </operator>");
            out.println("        <operator id=\"WilsonBaldingAndEmbed:gene" + (i+1) + "\" spec=\"speciesnetwork.operators.RebuildEmbedding\" speciesNetwork=\"@network:species\" taxonset=\"@taxonsuperset\" weight=\"5.0\">");
            out.println("            <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
            out.println("            <operator id=\"WilsonBalding:gene" + (i+1) + "\" spec=\"WilsonBalding\" tree=\"@tree:gene" + (i+1) + "\" weight=\"0.0\"/>");
            out.println("        </operator>");
        }
        // species network operators
        out.println("");
        out.println("        <operator id=\"divrRateScale:species\" spec=\"ScaleOperator\" parameter=\"@netDivRate:species\" scaleFactor=\"0.75\" weight=\"3.0\"/>");
        out.println("        <operator id=\"turnOverScale:species\" spec=\"ScaleOperator\" parameter=\"@turnOverRate:species\" scaleFactor=\"0.75\" weight=\"3.0\"/>");
        out.println("        <operator id=\"gammaProbUniform:species\" spec=\"speciesnetwork.operators.GammaProbUniform\" speciesNetwork=\"@network:species\" weight=\"10.0\"/>");
        out.println("        <operator id=\"gammaProbRndWalk:species\" spec=\"speciesnetwork.operators.GammaProbRndWalk\" speciesNetwork=\"@network:species\" weight=\"10.0\"/>\n");
        out.println("        <operator id=\"networkMultiplier:species\" spec=\"speciesnetwork.operators.NetworkMultiplier\" speciesNetwork=\"@network:species\" weight=\"15.0\">");
        for (int i = 0; i < nrOfGeneTrees; i++)
            out.println("            <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
        out.println("        </operator>");
        out.println("        <operator id=\"originMultiplier:species\" spec=\"speciesnetwork.operators.OriginMultiplier\" speciesNetwork=\"@network:species\" origin=\"@originTime:species\" weight=\"5.0\"/>");

        out.println("        <operator id=\"nodeUniformAndEmbed:species\" spec=\"speciesnetwork.operators.RebuildEmbedding\" speciesNetwork=\"@network:species\" taxonset=\"@taxonsuperset\" weight=\"" + 3*(nrOfGeneTrees+20) + "\">");
        for (int i = 0; i < nrOfGeneTrees; i++)
            out.println("            <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
        out.println("            <operator id=\"nodeUniform:species\" spec=\"speciesnetwork.operators.NodeUniform\" speciesNetwork=\"@network:species\" weight=\"0.0\"/>");
        out.println("        </operator>");
        out.println("        <operator id=\"nodeSliderAndEmbed:species\" spec=\"speciesnetwork.operators.RebuildEmbedding\" speciesNetwork=\"@network:species\" taxonset=\"@taxonsuperset\" weight=\"" + 4*(nrOfGeneTrees+20) + "\">");
        for (int i = 0; i < nrOfGeneTrees; i++)
            out.println("            <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
        out.println("            <operator id=\"nodeSlider:species\" spec=\"speciesnetwork.operators.NodeSlider\" speciesNetwork=\"@network:species\" origin=\"@originTime:species\" isNormal=\"true\" sigma=\"0.004\" weight=\"0.0\"/>");
        out.println("        </operator>");
        /*
        out.println("        <operator id=\"coordNodeUniformAndEmbed:species\" spec=\"speciesnetwork.operators.RebuildEmbedding\" speciesNetwork=\"@network:species\" taxonset=\"@taxonsuperset\" weight=\"60.0\">");
        for (int i = 0; i < nrOfGeneTrees; i++)
            out.println("            <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
        out.println("            <operator id=\"coordNodeUniform:species\" spec=\"speciesnetwork.operators.CoordinatedNodeUniform\" speciesNetwork=\"@network:species\" weight=\"0.0\">");
        for (int i = 0; i < nrOfGeneTrees; i++)
            out.println("                <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
        out.println("            </operator>");
        out.println("        </operator>");
        out.println("        <operator id=\"coordNodeSliderAndEmbed:species\" spec=\"speciesnetwork.operators.RebuildEmbedding\" speciesNetwork=\"@network:species\" taxonset=\"@taxonsuperset\" weight=\"120.0\">");
        for (int i = 0; i < nrOfGeneTrees; i++)
            out.println("            <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
        out.println("            <operator id=\"coordNodeSlider:species\" spec=\"speciesnetwork.operators.CoordinatedNodeSlider\" speciesNetwork=\"@network:species\" origin=\"@originTime:species\" isNormal=\"true\" sigma=\"0.004\" weight=\"0.0\">");
        for (int i = 0; i < nrOfGeneTrees; i++)
            out.println("                <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
        out.println("            </operator>");
        out.println("        </operator>"); */
        // whether or not to write network topology operators
        if (!networkOperatorInput.get())  out.println("        <!--");
        out.println("        <operator id=\"branchRelocateWideAndEmbed:species\" spec=\"speciesnetwork.operators.RebuildEmbedding\" speciesNetwork=\"@network:species\" taxonset=\"@taxonsuperset\" weight=\"200.0\">");
        for (int i = 0; i < nrOfGeneTrees; i++)
            out.println("            <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
        out.println("            <operator id=\"branchRelocatorW:species\" spec=\"speciesnetwork.operators.BranchRelocator\" speciesNetwork=\"@network:species\" isWide=\"true\" weight=\"0.0\"/>");
        out.println("        </operator>");
        out.println("        <operator id=\"branchRelocateNarrowAndEmbed:species\" spec=\"speciesnetwork.operators.RebuildEmbedding\" speciesNetwork=\"@network:species\" taxonset=\"@taxonsuperset\" weight=\"60.0\">");
        for (int i = 0; i < nrOfGeneTrees; i++)
            out.println("            <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
        out.println("            <operator id=\"branchRelocatorN:species\" spec=\"speciesnetwork.operators.BranchRelocator\" speciesNetwork=\"@network:species\" isWide=\"false\" weight=\"0.0\"/>");
        out.println("        </operator>");
        out.println("        <operator id=\"addReticulationAndEmbed:species\" spec=\"speciesnetwork.operators.RebuildEmbedding\" speciesNetwork=\"@network:species\" taxonset=\"@taxonsuperset\" weight=\"120.0\">");
        for (int i = 0; i < nrOfGeneTrees; i++)
            out.println("            <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
        out.println("            <operator id=\"addReticulation:species\" spec=\"speciesnetwork.operators.AddReticulation\" speciesNetwork=\"@network:species\" weight=\"0.0\"/>");
        out.println("        </operator>");
        out.println("        <operator id=\"deleteReticulationAndEmbed:species\" spec=\"speciesnetwork.operators.RebuildEmbedding\" speciesNetwork=\"@network:species\" taxonset=\"@taxonsuperset\" weight=\"120.0\">");
        for (int i = 0; i < nrOfGeneTrees; i++)
            out.println("            <geneTree idref=\"tree:gene" + (i+1) + "\"/>");
        out.println("            <operator id=\"deleteReticulation:species\" spec=\"speciesnetwork.operators.DeleteReticulation\" speciesNetwork=\"@network:species\" weight=\"0.0\"/>");
        out.println("        </operator>");
        if (!networkOperatorInput.get())  out.println("        -->");
        // print loggers
        out.println("");
        out.println("        <logger id=\"screenlog\" logEvery=\"20000\" model=\"@posterior\">");
        out.println("            <log idref=\"posterior\"/>");
        out.println("            <log id=\"ESS.0\" spec=\"util.ESS\" arg=\"@posterior\"/>");
        out.println("            <log idref=\"likelihood\"/>");
        out.println("            <log idref=\"prior\"/>");
        out.println("            <log idref=\"coalescent\"/>");
        out.println("        </logger>");
        out.println("        <logger id=\"tracelog\" fileName=\"" + outputFileName + ".trace.log\" logEvery=\"2000\" model=\"@posterior\" sort=\"smart\">");
        out.println("            <log idref=\"posterior\"/>");
        out.println("            <log idref=\"likelihood\"/>");
        out.println("            <log idref=\"prior\"/>");
        out.println("            <log idref=\"coalescent\"/>");
        out.println("            <log idref=\"netDivRate:species\"/>");
        out.println("            <log idref=\"turnOverRate:species\"/>");
        out.println("            <log idref=\"originTime:species\"/>");
        out.println("            <log id=\"height:species\" speciesNetwork=\"@network:species\" spec=\"speciesnetwork.NetworkStatLogger\"/>");
        for (int i = 0; i < nrOfGeneTrees; i++) {
            out.println("            <log id=\"height:gene" + (i+1) + "\" tree=\"@tree:gene" + (i+1) + "\" spec=\"beast.evolution.tree.TreeStatLogger\"/>");
        }
        out.println("        </logger>");
        out.println("        <logger id=\"specieslog\" fileName=\"" + outputFileName + ".species.trees\" logEvery=\"2000\" mode=\"tree\">");
        out.println("            <log id=\"networkLogger:species\" spec=\"speciesnetwork.NetworkWithMetaDataLogger\" speciesNetwork=\"@network:species\"/>");
        out.println("        </logger>");
        for (int i = 0; i < nrOfGeneTrees; i++) {
            out.println("        <logger id=\"treelog:gene" + (i + 1) + "\" fileName=\"" + outputFileName + ".gene" + (i+1) + ".trees\" logEvery=\"2000\" mode=\"tree\">");
            out.println("            <log id=\"treeLogger:gene" + (i+1) + "\" tree=\"@tree:gene" + (i+1) + "\" spec=\"beast.evolution.tree.TreeWithMetaDataLogger\"/>");
            out.println("        </logger>");
        }
        out.println("    </run>");  // end of MCMC
        out.println("</beast>");
    }

    private void writeGeneTrees(String outputFileName) throws IOException {
        if (outputFileName == null) {
            for (int i = 0; i < nrOfGeneTrees; i++) {
            	EmbeddedTree geneTree = geneTrees.get(i);
                System.out.println(geneTree.getRoot().toNewick() + ";");
            }
        } else {
            FileWriter fw = new FileWriter(outputFileName, true);
            for (int i = 0; i < nrOfGeneTrees; i++) {
            	EmbeddedTree geneTree = geneTrees.get(i);
                fw.write(geneTree.getRoot().toNewick() + ";\n");
            }
            fw.close();
        }
    }
}
