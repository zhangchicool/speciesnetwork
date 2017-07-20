package speciesnetwork;

import java.util.List;

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

    // the coalescent times of this gene tree for all species branches
    protected ListMultimap<Integer, Double> coalescentTimes = ArrayListMultimap.create();
    protected ListMultimap<Integer, Double> storedCoalescentTimes = ArrayListMultimap.create();
    // the number of lineages at the tipward end of each species branch
    protected Multiset<Integer> coalescentLineageCounts = HashMultiset.create();
    protected Multiset<Integer> storedCoalescentLineageCounts = HashMultiset.create();

    protected double[][] speciesOccupancy;
    protected double[][] storedSpeciesOccupancy;
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

        storedSpeciesOccupancy = new double[speciesOccupancy.length][speciesOccupancy[0].length];
        System.arraycopy(speciesOccupancy, 0, storedSpeciesOccupancy, 0, speciesOccupancy.length);
        
        storedLogGammaSum = logGammaSum;

        super.store();
    }

    @Override
    public void restore() {
        ListMultimap<Integer, Double> tmpCoalescentTimes = coalescentTimes;
        Multiset<Integer> tmpCoalescentLineageCounts = coalescentLineageCounts;
        double[][] tmpSpeciesOccupancy = speciesOccupancy;
        double tmpLogGammaSum = logGammaSum;

        coalescentTimes = storedCoalescentTimes;
        coalescentLineageCounts = storedCoalescentLineageCounts;
        speciesOccupancy = storedSpeciesOccupancy;
        logGammaSum = storedLogGammaSum;

        storedCoalescentTimes = tmpCoalescentTimes;
        storedCoalescentLineageCounts = tmpCoalescentLineageCounts;
        storedSpeciesOccupancy = tmpSpeciesOccupancy;
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
        speciesNetwork = speciesNetworkInput.get();
        geneTree = geneTreeInput.get();
        logGammaSum = 0.0;

        final int geneTreeNodeCount = geneTree.getNodeCount();
        final int speciesBranchCount = speciesNetwork.getBranchCount();
        speciesOccupancy = new double[geneTreeNodeCount][speciesBranchCount];

        // reset coalescent arrays as these values need to be recomputed after any changes to the species or gene tree
        coalescentLineageCounts.clear();
        coalescentTimes.clear();

        final Node geneTreeRoot = geneTree.getRoot();
        final NetworkNode speciesNetworkRoot = speciesNetwork.getRoot();
        final Integer speciesRootBranchNumber = speciesNetworkRoot.gammaBranchNumber;
        try {
            recurseCoalescentEvents(geneTreeRoot, speciesNetworkRoot, speciesRootBranchNumber, Double.POSITIVE_INFINITY);
        } catch (Exception e) {
            e.printStackTrace();
        }

        needsUpdate = false;
    }

    // forward in time recursion, unlike StarBEAST 2
    private void recurseCoalescentEvents(final Node geneTreeNode, final NetworkNode speciesNetworkNode,
                                         final Integer speciesBranchNumber, final double lastHeight) throws Exception {
        final double geneNodeHeight = geneTreeNode.getHeight();
        final double speciesNodeHeight = speciesNetworkNode.getHeight();
        final int geneTreeNodeNumber = geneTreeNode.getNr();

        // check if coalescent node occurs in a descendant species network branch
        if (geneNodeHeight < speciesNodeHeight) {
            speciesOccupancy[geneTreeNodeNumber][speciesBranchNumber] += lastHeight - speciesNodeHeight;
            coalescentLineageCounts.add(speciesBranchNumber);
            if (speciesNetworkNode.isReticulation()) {
                final double gamma = speciesNetworkNode.inheritProb;
                if (speciesNetworkNode.gammaBranchNumber.equals(speciesBranchNumber)) {
                    logGammaSum += Math.log(gamma);
                } else {
                    logGammaSum += Math.log(1.0 - gamma);
                }
            }
            // traversal direction forward in time
            final int traversalNodeNumber = speciesNetworkNode.getTraversalNumber();
            final Integer nextSpeciesBranchNumber = geneTree.getEmbedding(geneTreeNodeNumber, traversalNodeNumber);
            assert (nextSpeciesBranchNumber >= 0);
            final NetworkNode nextSpeciesNode = speciesNetworkNode.getChildByBranch(nextSpeciesBranchNumber);
            assert nextSpeciesNode != null;
            recurseCoalescentEvents(geneTreeNode, nextSpeciesNode, nextSpeciesBranchNumber, speciesNodeHeight);
        } else if (geneTreeNode.isLeaf()) { // assumes tip node heights are always zero
            speciesOccupancy[geneTreeNodeNumber][speciesBranchNumber] += lastHeight;
            coalescentLineageCounts.add(speciesBranchNumber);
        } else {
            speciesOccupancy[geneTreeNodeNumber][speciesBranchNumber] += lastHeight - geneNodeHeight;
            coalescentTimes.put(speciesBranchNumber, geneNodeHeight);
            final Node leftChild = geneTreeNode.getLeft();
            final Node rightChild = geneTreeNode.getRight();
            recurseCoalescentEvents(leftChild, speciesNetworkNode, speciesBranchNumber, geneNodeHeight);
            recurseCoalescentEvents(rightChild, speciesNetworkNode, speciesBranchNumber, geneNodeHeight);
        }
    }

    public double[][] getSpeciesOccupancy() {
        if (needsUpdate) update();
        return speciesOccupancy;
    }

    /**
     * @return the first tip node which is descendant of
     * @param gTreeNode
     * this can be in Tree.java as gTreeNode.getGeneTreeTipDescendant()
     */
    public Node getGeneNodeDescendantTip(Node gTreeNode) {
        final List<Node> gTreeTips = geneTree.getExternalNodes();  // tips
        for (Node tip : gTreeTips) {
            Node node = tip;
            while(node != null && !node.equals(gTreeNode)) {
                node = node.getParent();
            }
            if (node != null)
                return tip;  // find you!
        }
        return null;  // looped all the tips but nothing found
    }
}
