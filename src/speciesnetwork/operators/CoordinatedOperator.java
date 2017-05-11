package speciesnetwork.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.evolution.tree.Node;
import speciesnetwork.EmbeddedTree;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Chi Zhang
 */

public abstract class CoordinatedOperator extends Operator {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Input.Validate.REQUIRED);
    public final Input<List<EmbeddedTree>> geneTreesInput = new Input<>("geneTree",
            "The gene tree within the species network.", new ArrayList<>());

    @Override
    public void initAndValidate() {
    }

    /**
     * The RubberBand algorithm of Rannala & Yang, 2003 (Appendix Step 4)
     */
    protected double updateRubberBand(NetworkNode networkNode, final double oldHeight, final double newHeight,
                                      final double lower, final double upper) {
        if (!networkNode.isSpeciation())  // necessary only when speciation node
            return 0.0;

        final Integer speciesBrNr = networkNode.gammaBranchNumber;
        final NetworkNode parentNode = networkNode.getParentByBranch(speciesBrNr);

        final int nLoci = geneTreesInput.get().size();  // if nLoci == 0 (no input), gene trees are not changed

        final List<EmbeddedTree> geneTrees = geneTreesInput.get();
        int m = 0;  // # gene node heights changed relative to 'upper'
        int n = 0;  // # gene node heights changed relative to 'lower'
        for (int i = 0; i < nLoci; i++) {  // loop over all loci
            EmbeddedTree geneTree = geneTrees.get(i);

            // update the gene tree node heights  TODO: will this introduce negative branch lengths?
            for (Node gNode : geneTree.getInternalNodes()) {
                final double gNodeHeight = gNode.getHeight();

                if (oldHeight <= gNodeHeight && gNodeHeight < upper) {
                    if (networkNode.isRoot() || isWithinChildBranch(parentNode, gNode, geneTree, upper)) {
                        // update the node height relative to 'upper'
                        final double gNewNodeHeight = upper - (upper - gNodeHeight) * (upper - newHeight) / (upper - oldHeight);
                        gNode.setHeight(gNewNodeHeight);
                        m++;
                    }
                } else if (lower < gNodeHeight && gNodeHeight < oldHeight) {
                    if (isWithinChildBranch(networkNode, gNode, geneTree, upper)) {
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

    // check if a gene tree node is within a child branch of the network node
    private boolean isWithinChildBranch(NetworkNode snNode, Node gNode, EmbeddedTree geneTree, double upper) {
        final int traversalNodeNr = snNode.getTraversalNumber();
        Integer withinBrNr;  Node ancNode = gNode;
        do {
            withinBrNr = geneTree.getEmbedding(ancNode.getNr(), traversalNodeNr);
            ancNode = ancNode.getParent();
        } while (withinBrNr < 0 && ancNode != null && ancNode.getHeight() < upper);

        return snNode.childBranchNumbers.contains(withinBrNr);
    }
}
