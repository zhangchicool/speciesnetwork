package speciesnetwork.operators;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import beast.evolution.tree.Node;
import speciesnetwork.EmbeddedTree;
import speciesnetwork.Embedding;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;
import speciesnetwork.SanityChecks;

/**
 * @author Chi Zhang
 */

@Description("Randomly select an internal network node and move its height using a sliding window." +
        "Also update the affected gene tree node heights at each locus to maintain the compatibility.")
public class CoordinatedNodeSlider extends Operator {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);
    public final Input<List<EmbeddedTree>> geneTreesInput = new Input<>("geneTree",
            "The gene tree within the species network.", new ArrayList<>());
    public final Input<RealParameter> originInput =
            new Input<>("origin", "The time when the process started.", Validate.REQUIRED);
    public final Input<Boolean> isNormalInput =
            new Input<>("isNormal", "Using normal proposal (default).", true);
    public final Input<Double> sigmaInput =
            new Input<>("sigma", "Standard deviation of the normal proposal (default is 0.01).", 0.01);
    public final Input<Double> windowSizeInput =
            new Input<>("windowSize", "Window size of the uniform proposal (default is 0.02).", 0.02);

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        final Network speciesNetwork = speciesNetworkInput.get();

        // pick an internal node randomly, including origin
        final NetworkNode[] internalNodes = speciesNetwork.getInternalNodesWithOrigin();
        final int randomIndex = Randomizer.nextInt(internalNodes.length);
        final NetworkNode pickedNode = internalNodes[randomIndex];

        // determine the lower and upper bounds
        double upper = Double.MAX_VALUE;
        for (NetworkNode p: pickedNode.getParents()) {
            upper = Math.min(upper, p.getHeight());
        }
        double lower = 0.0;
        for (NetworkNode c: pickedNode.getChildren()) {
            lower = Math.max(lower, c.getHeight());
        }
        if (lower >= upper)
            throw new RuntimeException("Developer ERROR: lower bound >= upper bound!");

        // propose a new height, reflect it back if it's outside the boundary
        double oldHeight = pickedNode.getHeight();
        double newHeight;
        if (isNormalInput.get()) {
            final double sigma = sigmaInput.get();
            newHeight = oldHeight + Randomizer.nextGaussian() * sigma;
        } else {
            final double windowSize = windowSizeInput.get();
            newHeight = oldHeight + (Randomizer.nextDouble() - 0.5) * windowSize;
        }
        while (newHeight < lower || newHeight > upper) {
            if (newHeight < lower)
                newHeight = 2.0 * lower - newHeight;
            if (newHeight > upper)
                newHeight = 2.0 * upper - newHeight;
        }

        // update the new node height
        if (pickedNode.isOrigin()) {
            final RealParameter originTime = originInput.get();
            if (outsideBounds(newHeight, originTime))
                return Double.NEGATIVE_INFINITY;
            originTime.setValue(newHeight);
        }
        speciesNetwork.startEditing(this);
        pickedNode.setHeight(newHeight);
        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());

        // update gene tree node heights, return proposal ratio of this update
        // this move does not break gene tree embeddings thus no need to rebuild
        return updateRubberBand(pickedNode, oldHeight, newHeight, lower, upper);
    }

    private boolean outsideBounds(final double value, final RealParameter param) {
        final Double l = param.getLower();
        final Double h = param.getUpper();

        return (value < l || value > h);
    }

    /**
     * The RubberBand algorithm of Rannala & Yang, 2003 (Appendix Step 4)
     */
    protected double updateRubberBand(NetworkNode networkNode, final double oldHeight, final double newHeight,
                                      final double lower, final double upper) {
        final List<EmbeddedTree> geneTrees = geneTreesInput.get();

        int m = 0;  // # gene node heights changed relative to 'upper'
        int n = 0;  // # gene node heights changed relative to 'lower'
        if (networkNode.isSpeciation()) {
            final Integer speciesBrNr = networkNode.gammaBranchNumber;
            final NetworkNode parentNode = networkNode.getParentByBranch(speciesBrNr);
            for (EmbeddedTree geneTree : geneTrees) {
                Embedding embedding = geneTree.embedding;
                geneTree.startEditing(this);  // *all* gene trees will be edited
                // update the gene tree node heights
                for (Node gNode : geneTree.getInternalNodes()) {
                    final double gNodeHeight = gNode.getHeight();
                    if (oldHeight <= gNodeHeight && gNodeHeight < upper) {
                        if (networkNode.isRoot() ||
                            isWithinChildBranch(parentNode, Collections.singletonList(speciesBrNr), gNode, embedding, upper)) {
                            // update the node height relative to 'upper'
                            final double gNewNodeHeight = upper - (upper - gNodeHeight) * (upper - newHeight) / (upper - oldHeight);
                            gNode.setHeight(gNewNodeHeight);
                            m++;
                        }
                    } else if (lower < gNodeHeight && gNodeHeight < oldHeight) {
                        if (isWithinChildBranch(networkNode, networkNode.childBranchNumbers, gNode, embedding, oldHeight)) {
                            // update the node height relative to 'lower'
                            final double gNewNodeHeight = lower + (gNodeHeight - lower) * (newHeight - lower) / (oldHeight - lower);
                            gNode.setHeight(gNewNodeHeight);
                            n++;
                        }
                    }
                }
            }
        }
        else if (networkNode.isReticulation()) {
            final Integer snLeftBrNr = networkNode.gammaBranchNumber;
            final Integer snRightBrNr = networkNode.gammaBranchNumber + 1;
            final NetworkNode parentLNode = networkNode.getParentByBranch(snLeftBrNr);
            final NetworkNode parentRNode = networkNode.getParentByBranch(snRightBrNr);
            for (EmbeddedTree geneTree : geneTrees) {
                Embedding embedding = geneTree.embedding;
                geneTree.startEditing(this);  // *all* gene trees will be edited
                // update the gene tree node heights
                for (Node gNode : geneTree.getInternalNodes()) {
                    final double gNodeHeight = gNode.getHeight();
                    if (oldHeight <= gNodeHeight && gNodeHeight < upper) {
                        if (isWithinChildBranch(parentLNode, Collections.singletonList(snLeftBrNr),  gNode, embedding, parentLNode.getHeight()) ||
                            isWithinChildBranch(parentRNode, Collections.singletonList(snRightBrNr), gNode, embedding, parentRNode.getHeight())) {
                            // update the node height relative to 'upper'
                            final double gNewNodeHeight = upper - (upper - gNodeHeight) * (upper - newHeight) / (upper - oldHeight);
                            gNode.setHeight(gNewNodeHeight);
                            m++;
                        }
                    } else if (lower < gNodeHeight && gNodeHeight < oldHeight) {
                        if (isWithinChildBranch(networkNode, networkNode.childBranchNumbers, gNode, embedding, oldHeight)) {
                            // update the node height relative to 'lower'
                            final double gNewNodeHeight = lower + (gNodeHeight - lower) * (newHeight - lower) / (oldHeight - lower);
                            gNode.setHeight(gNewNodeHeight);
                            n++;
                        }
                    }
                }
            }
        }

        return m * Math.log((upper - newHeight)/(upper - oldHeight)) +
               n * Math.log((newHeight - lower)/(oldHeight - lower));
    }

    /* check if a gene tree node is within certain child branches of the network node */
    private boolean isWithinChildBranch(NetworkNode snNode, List<Integer> childBrNrs, Node gNode,
                                        Embedding embedding, final double upper) {
        final int traversalNodeNr = snNode.getTraversalNumber();
        int withinBrNr;
        Node treNode = gNode;
        do {
            withinBrNr = embedding.getDirection(treNode.getNr(), traversalNodeNr);
            treNode = treNode.getParent();
        }
        while (withinBrNr < 0 && treNode != null && treNode.getHeight() < upper);

        return childBrNrs.contains(withinBrNr);
    }
}
