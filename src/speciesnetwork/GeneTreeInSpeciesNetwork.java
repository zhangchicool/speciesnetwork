package speciesnetwork;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multiset;
import com.google.common.collect.HashMultiset;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;

/**
 * @author Huw Ogilvie
 * @author Chi Zhang
 */

public class GeneTreeInSpeciesNetwork extends CalculationNode {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "Species network for embedding the gene tree.", Validate.REQUIRED);
    public final Input<EmbeddedTree> geneTreeInput =
            new Input<>("geneTree", "Gene tree embedded in the species network.", Validate.REQUIRED);
    public final Input<Double> ploidyInput =
            new Input<>("ploidy", "Ploidy (copy number) for this gene (default is 2).", 2.0);
    protected double ploidy;

    private boolean needsUpdate;
    private EmbeddedTree geneTree;
    private Network speciesNetwork;

    // the coalescent times of this gene tree in each species branch
    protected ListMultimap<Integer, Double> coalescentTimes = ArrayListMultimap.create();
    protected ListMultimap<Integer, Double> storedCoalescentTimes = ArrayListMultimap.create();
    // the number of lineages at the tipward end of each species branch
    protected Multiset<Integer> coalescentLineageCounts = HashMultiset.create();
    protected Multiset<Integer> storedCoalescentLineageCounts = HashMultiset.create();
    // sum of the log inheritance probabilities (log(Lambda), part of log[f(G|Psi)])
    protected double logGammaSum;
    protected double storedLogGammaSum;

    @Override
    public boolean requiresRecalculation() {
        needsUpdate = geneTreeInput.isDirty() || speciesNetworkInput.isDirty();
        return needsUpdate;
    }

    @Override
    public void store() {
        storedCoalescentTimes.clear();
        storedCoalescentLineageCounts.clear();

        storedCoalescentTimes.putAll(coalescentTimes);
        storedCoalescentLineageCounts.addAll(coalescentLineageCounts);
        storedLogGammaSum = logGammaSum;

        super.store();
    }

    @Override
    public void restore() {
        ListMultimap<Integer, Double> tmpCoalescentTimes = coalescentTimes;
        Multiset<Integer> tmpCoalescentLineageCounts = coalescentLineageCounts;
        double tmpLogGammaSum = logGammaSum;

        coalescentTimes = storedCoalescentTimes;
        coalescentLineageCounts = storedCoalescentLineageCounts;
        logGammaSum = storedLogGammaSum;

        storedCoalescentTimes = tmpCoalescentTimes;
        storedCoalescentLineageCounts = tmpCoalescentLineageCounts;
        storedLogGammaSum = tmpLogGammaSum;

        super.restore();
    }

    public void initAndValidate() {
        ploidy = ploidyInput.get();
        geneTree = geneTreeInput.get();
        speciesNetwork = speciesNetworkInput.get();
        needsUpdate = true;
    }

    protected void computeCoalescentTimes() {
        if (needsUpdate) {
            update();
        }
    }

    private void update() {
        // reset coalescent arrays as these values need to be recomputed after any changes to the species or gene tree
        coalescentLineageCounts.clear();
        coalescentTimes.clear();
        logGammaSum = 0.0;

        geneTree = geneTreeInput.get();
        final Node geneTreeRoot = geneTree.getRoot();
        speciesNetwork = speciesNetworkInput.get();
        final NetworkNode speciesNetworkRoot = speciesNetwork.getRoot();
        final Integer speciesRootBranchNumber = speciesNetworkRoot.gammaBranchNumber;
        /* The recursion starts from the root of gene tree and root of species network, and moves forward in time.
           Typically, the root age of gene tree is larger than the root age of species network, but it is not always
           the case due to reticulations in the network or incomplete sampling of individuals in the gene tree. */
        try {
            recurseCoalescentEvents(geneTreeRoot, speciesNetworkRoot, speciesRootBranchNumber, Double.POSITIVE_INFINITY);
        } catch (Exception e) {
            e.printStackTrace();
        }

        needsUpdate = false;
    }

    private void recurseCoalescentEvents(final Node geneTreeNode, final NetworkNode speciesNetworkNode,
                                         final Integer speciesBranchNumber, final double lastHeight) {
        final double geneNodeHeight = geneTreeNode.getHeight();
        final double speciesNodeHeight = speciesNetworkNode.getHeight();
        final int geneTreeNodeNumber = geneTreeNode.getNr();

        if (geneTreeNode.isLeaf() && speciesNetworkNode.isLeaf()) {
            // reach the tip with height >= 0, gene tree tip height == species tip height
            // speciesOccupancy[geneTreeNodeNumber][speciesBranchNumber] = lastHeight - speciesNodeHeight;
            coalescentLineageCounts.add(speciesBranchNumber);
        }
        else if (geneNodeHeight <= speciesNodeHeight) {
            // current gene tree node occurs in a descendant branch of current species node
            // speciesOccupancy[geneTreeNodeNumber][speciesBranchNumber] = lastHeight - speciesNodeHeight;
            coalescentLineageCounts.add(speciesBranchNumber);
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
            final Integer nextSpeciesBranchNumber = geneTree.embedding.getDirection(geneTreeNodeNumber, traversalNodeNumber);
            assert (nextSpeciesBranchNumber >= 0);
            final NetworkNode nextSpeciesNode = speciesNetworkNode.getChildByBranch(nextSpeciesBranchNumber);
            assert nextSpeciesNode != null;
            recurseCoalescentEvents(geneTreeNode, nextSpeciesNode, nextSpeciesBranchNumber, speciesNodeHeight);
        }
        else {
            // current gene tree node occurs above current species node
            // speciesOccupancy[geneTreeNodeNumber][speciesBranchNumber] = lastHeight - geneNodeHeight;
            coalescentTimes.put(speciesBranchNumber, geneNodeHeight);
            // move on to the descendant gene tree nodes (traversal direction forward in time)
            for (Node geneChildNode : geneTreeNode.getChildren()) {
                recurseCoalescentEvents(geneChildNode, speciesNetworkNode, speciesBranchNumber, geneNodeHeight);
            }
        }
    }
}
