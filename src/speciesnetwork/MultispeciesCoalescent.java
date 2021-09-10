package speciesnetwork;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multiset;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;

/**
 * @author Huw Ogilvie
 * @author Chi Zhang
 */

@Description("Calculates probability of gene trees conditioned on a species network (multispecies network coalescent, MSNC).")
public class MultispeciesCoalescent extends Distribution {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);
    public final Input<List<EmbeddedTree>> geneTreesInput =
            new Input<>("geneTree", "Gene tree embedded in the species network.", new ArrayList<>());
    public final Input<RealParameter> ploidiesInput =
            new Input<>("ploidy", "Ploidy (copy number) array for all genes (default is 2).");
    public final Input<PopulationSizeModel> populationModelInput =
            new Input<>("populationModel", "The species network population model.", Validate.REQUIRED);

    private int nGeneTrees;
    private RealParameter ploidies;

    // list of map(species branch -> [coalescent times]) for each gene
    private final List<ListMultimap<Integer, Double>> coalescentTimes = new ArrayList<>();
    // list of set(species branch number) counting lineages at the tipward end of each species branch
    private final List<Multiset<Integer>> bottomLineageCounts = new ArrayList<>();
    // sum of the log inheritance probabilities (log(Lambda), part of log[f(G|Psi)])
    private double logGammaSum;

    private final List<int[]> allLineageCounts = new ArrayList<>();
    private final List<int[]> allEventCounts = new ArrayList<>();
    private final List<List<Double[]>> allCoalescentTimes = new ArrayList<>();

    @Override
    public void initAndValidate() {
        final List<EmbeddedTree> geneTrees = geneTreesInput.get();
        // sanity check
        if (geneTrees == null)
            throw new RuntimeException("Check gene tree input!");
        nGeneTrees = geneTrees.size();

        if ((ploidies = ploidiesInput.get()) == null)
            ploidies = new RealParameter("2.0");  // default
        ploidies.setDimension(nGeneTrees);

        final Network speciesNetwork = speciesNetworkInput.get();
        final int speciesBranchCount = speciesNetwork.getBranchCount();
        final PopulationSizeModel populationModel = populationModelInput.get();
        populationModel.initPopSizes(speciesBranchCount);
    }

    @Override
    public double calculateLogP() {
        logP = coalescentProb();
        return logP;
    }

    /**
     * @return the coalescent probability of gene trees embedded in the species network
     */
    public double coalescentProb() {
        final Network speciesNetwork = speciesNetworkInput.get();
        final List<EmbeddedTree> geneTrees = geneTreesInput.get();
        // SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin()); // species network should not be insane

        final int speciesBranchCount = speciesNetwork.getBranchCount();
        final NetworkNode speciesNetworkRoot = speciesNetwork.getRoot();
        final Integer speciesRootBranchNumber = speciesNetworkRoot.gammaBranchNumber;

        coalescentTimes.clear();
        bottomLineageCounts.clear();
        logGammaSum = 0.0;
        // collect the coalescent times of each gene tree in each species branch, also calculate log(Lambda)
        for (EmbeddedTree geneTree : geneTrees) {
            final Embedding embedding = geneTree.embedding;
            final Node geneTreeRoot = geneTree.getRoot();

            ListMultimap<Integer, Double> coalescentTimesJ = ArrayListMultimap.create();
            Multiset<Integer> bottomLineageCountsJ = HashMultiset.create();  // empty
            /* The recursion starts from the root of gene tree and root of species network, and moves forward in time.
               Typically, the root age of gene tree is larger than the root age of species network, but it is not always
               the case due to reticulations in the network or incomplete sampling of individuals in the gene tree. */
            try {
                recurseCoalescentEvents(geneTreeRoot, speciesNetworkRoot, speciesRootBranchNumber, Double.POSITIVE_INFINITY,
                                        embedding, coalescentTimesJ, bottomLineageCountsJ);
            } catch (Exception e) {
                e.printStackTrace();
            }

            coalescentTimes.add(coalescentTimesJ);
            bottomLineageCounts.add(bottomLineageCountsJ);
        }

        // reset all che counters
        allLineageCounts.clear();
        allEventCounts.clear();
        allCoalescentTimes.clear();
        for (int i = 0; i < speciesBranchCount; i++) {
            allLineageCounts.add(new int[nGeneTrees]);
            allEventCounts.add(new int[nGeneTrees]);
            allCoalescentTimes.add(new ArrayList<>());
        }

        // transpose gene-branch list of lists to branch-gene list of lists
        final double[] genePloidy = new double[nGeneTrees];
        for (int j = 0; j < nGeneTrees; j++) {  // gene tree "j"
            for (int i = 0; i < speciesBranchCount; i++) {  // species network branch "i"
                // number of lineages at the tipward end of species branch "i"
                final int lineageCount = bottomLineageCounts.get(j).count(i);

                // number of coalescent events in species branch "i"
                final List<Double> timesView = coalescentTimes.get(j).get(i);
                final int eventCount = timesView.size();

                // add branch start and end times to the coalescent times
                final Double[] coalTimes = new Double[eventCount + 2];
                timesView.toArray(coalTimes);
                final int networkNodeNr = speciesNetwork.getNodeNumber(i);
                final NetworkNode snNode = speciesNetwork.getNode(networkNodeNr);
                final NetworkNode parentNode = snNode.getParentByBranch(i);
                if (parentNode.isOrigin())
                    coalTimes[eventCount] = Double.POSITIVE_INFINITY;
                else
                    coalTimes[eventCount] = parentNode.getHeight();
                coalTimes[eventCount + 1] = snNode.getHeight();
                // sort times of coalescent events in ascending order
                Arrays.sort(coalTimes);

                // collect things together
                allEventCounts.get(i)[j] = eventCount;
                allLineageCounts.get(i)[j] = lineageCount;
                allCoalescentTimes.get(i).add(coalTimes);
            }
            genePloidy[j] = ploidies.getValue(j);
        }

        // now calculate coalescent prob. by looping over the species branches
        final PopulationSizeModel populationModel = populationModelInput.get();

        double logProb = logGammaSum;
        for (int i = 0; i < speciesBranchCount; i++) {
            final List<Double[]> branchCoalescentTimes = allCoalescentTimes.get(i);
            final int[] branchLineageCounts = allLineageCounts.get(i);
            final int[] branchEventCounts = allEventCounts.get(i);

            logProb += populationModel.branchLogP(i, genePloidy, branchCoalescentTimes, branchLineageCounts, branchEventCounts);
        }

        return logProb;
    }

    private void recurseCoalescentEvents(Node geneTreeNode, NetworkNode speciesNetworkNode, Integer speciesBranchNumber, double lastHeight,
                                         Embedding embedding, ListMultimap<Integer, Double> coalTimes, Multiset<Integer> bottomBrNrs) {
        final double geneNodeHeight = geneTreeNode.getHeight();
        final double speciesNodeHeight = speciesNetworkNode.getHeight();
        final int geneTreeNodeNumber = geneTreeNode.getNr();

        if (geneTreeNode.isLeaf() && speciesNetworkNode.isLeaf()) {
            // reach the tip with height >= 0, gene tree tip height == species tip height
            // speciesOccupancy[geneTreeNodeNumber][speciesBranchNumber] = lastHeight - speciesNodeHeight;
            bottomBrNrs.add(speciesBranchNumber);
        }
        else if (geneNodeHeight <= speciesNodeHeight) {
            // current gene tree node occurs in a descendant branch of current species node
            // speciesOccupancy[geneTreeNodeNumber][speciesBranchNumber] = lastHeight - speciesNodeHeight;
            bottomBrNrs.add(speciesBranchNumber);
            if (speciesNetworkNode.isReticulation()) {
                final double gamma = speciesNetworkNode.inheritProb;
                if (speciesNetworkNode.gammaBranchNumber.equals(speciesBranchNumber)) {
                    logGammaSum += Math.log(gamma);
                } else {
                    logGammaSum += Math.log(1.0 - gamma);
                }
            }
            // move on to the descendant species node (traversal direction forward in time)
            final int traversalNodeNumber = speciesNetworkNode.getTraversalNumber();
            final Integer nextSpeciesBranchNumber = embedding.getDirection(geneTreeNodeNumber, traversalNodeNumber);
            assert (nextSpeciesBranchNumber >= 0);
            final NetworkNode nextSpeciesNode = speciesNetworkNode.getChildByBranch(nextSpeciesBranchNumber);
            assert nextSpeciesNode != null;
            recurseCoalescentEvents(geneTreeNode, nextSpeciesNode, nextSpeciesBranchNumber, speciesNodeHeight,
                                    embedding, coalTimes, bottomBrNrs);
        } else {
            // current gene tree node occurs above current species node
            // speciesOccupancy[geneTreeNodeNumber][speciesBranchNumber] = lastHeight - geneNodeHeight;
            coalTimes.put(speciesBranchNumber, geneNodeHeight);
            // move on to the descendant gene tree nodes (traversal direction forward in time)
            for (Node geneChildNode : geneTreeNode.getChildren()) {
                recurseCoalescentEvents(geneChildNode, speciesNetworkNode, speciesBranchNumber, geneNodeHeight,
                                        embedding, coalTimes, bottomBrNrs);
            }
        }
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
    }
}
