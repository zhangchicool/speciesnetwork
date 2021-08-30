package speciesnetwork.operators;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.util.Randomizer;
import beast.evolution.tree.Node;
import speciesnetwork.EmbeddedTree;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;
import speciesnetwork.SanityChecks;

/**
 * @author Chi Zhang
 */

@Description("Randomly select an internal network node and move its height uniformly." +
        "Also update the affected gene tree node heights at each locus to maintain the compatibility.")
public class CoordinatedNodeUniform extends Operator {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);
    public final Input<List<EmbeddedTree>> geneTreesInput = new Input<>("geneTree",
            "The gene tree within the species network.", new ArrayList<>());

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        final Network speciesNetwork = speciesNetworkInput.get();

        // pick an internal node randomly
        final NetworkNode[] internalNodes = speciesNetwork.getInternalNodes();
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

        // propose a new height uniformly
        double oldHeight = pickedNode.getHeight();
        double newHeight = Randomizer.nextDouble() * (upper - lower) + lower;
        speciesNetwork.startEditing(this);
        pickedNode.setHeight(newHeight);
        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());

        // update gene tree node heights (necessary only when speciation network node)
        // return proposal ratio of this update
        return updateRubberBand(pickedNode, oldHeight, newHeight, lower, upper);
    }

    /**
     * The RubberBand algorithm of Rannala & Yang, 2003 (Appendix Step 4)
     */
    protected double updateRubberBand(NetworkNode networkNode, final double oldHeight, final double newHeight,
                                      final double lower, final double upper) {
        if (!networkNode.isSpeciation())  // necessary only when speciation network node
            return 0.0;

        final Integer speciesBrNr = networkNode.gammaBranchNumber;
        final NetworkNode parentNode = networkNode.getParentByBranch(speciesBrNr);

        final List<EmbeddedTree> geneTrees = geneTreesInput.get();
        final int nLoci = geneTrees.size();  // if nLoci == 0 (no input), gene trees are not changed

        int m = 0;  // # gene node heights changed relative to 'upper'
        int n = 0;  // # gene node heights changed relative to 'lower'
        for (EmbeddedTree geneTree : geneTrees) {  // loop over all loci
            geneTree.startEditing(this);  // Tell BEAST that *all* gene trees will be edited

            // update the gene tree node heights
            for (Node gNode : geneTree.getInternalNodes()) {
                final double gNodeHeight = gNode.getHeight();

                if (oldHeight <= gNodeHeight && gNodeHeight < upper) {
                    if (networkNode.isRoot() ||
                            isWithinChildBranch(parentNode, Collections.singletonList(speciesBrNr), gNode, geneTree, upper)) {
                        // update the node height relative to 'upper'
                        final double gNewNodeHeight = upper - (upper - gNodeHeight) * (upper - newHeight) / (upper - oldHeight);
                        gNode.setHeight(gNewNodeHeight);
                        m++;
                    }
                } else if (lower < gNodeHeight && gNodeHeight < oldHeight) {
                    if (isWithinChildBranch(networkNode, networkNode.childBranchNumbers, gNode, geneTree, upper)) {
                        // update the node height relative to 'lower'
                        final double gNewNodeHeight = lower + (gNodeHeight - lower) * (newHeight - lower) / (oldHeight - lower);
                        gNode.setHeight(gNewNodeHeight);
                        n++;
                    }
                }
            }
        }

        return m * Math.log((upper - newHeight)/(upper - oldHeight)) +
               n * Math.log((newHeight - lower)/(oldHeight - lower));
    }

    /* check if a gene tree node is within certain child branches of the network node */
    private boolean isWithinChildBranch(NetworkNode snNode, List<Integer> childBrNrs,
                                        Node gNode, EmbeddedTree geneTree, final double upper) {
        final int traversalNodeNr = snNode.getTraversalNumber();
        int withinBrNr;
        Node treNode = gNode;
        do {
            withinBrNr = geneTree.embedding.getDirection(treNode.getNr(), traversalNodeNr);
            treNode = treNode.getParent();
        }
        while (withinBrNr < 0 && treNode != null && treNode.getHeight() < upper);

        return childBrNrs.contains(withinBrNr);
    }
}
