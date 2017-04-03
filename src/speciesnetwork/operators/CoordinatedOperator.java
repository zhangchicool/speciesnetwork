package speciesnetwork.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.IntegerParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Chi Zhang
 */

public abstract class CoordinatedOperator extends Operator {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Input.Validate.REQUIRED);
    public final Input<List<Tree>> geneTreesInput =
            new Input<>("geneTree", "Gene tree within the species network.", new ArrayList<>());
    public final Input<List<IntegerParameter>> embeddingsInput =
            new Input<>("embedding", "The embedding matrix.", new ArrayList<>());

    private int nLoci;

    @Override
    public void initAndValidate() {
        nLoci = geneTreesInput.get().size();  // if nLoci == 0 (no input), gene trees are not changed
        if (embeddingsInput.get().size() != nLoci) {
            throw new RuntimeException("Number of embeddings doesn't match number of gene trees!");
        }
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

        final List<Tree> geneTrees = geneTreesInput.get();
        final List<IntegerParameter> embeddings = embeddingsInput.get();
        int m = 0;  // # gene node heights changed relative to 'upper'
        int n = 0;  // # gene node heights changed relative to 'lower'
        for (int i = 0; i < nLoci; i++) {  // loop over all loci
            Tree gTree = geneTrees.get(i);
            IntegerParameter embedding = embeddings.get(i);

            // update the gene tree node heights  TODO: will this introduce negative branch lengths?
            for (Node gNode : gTree.getInternalNodes()) {
                final double gNodeHeight = gNode.getHeight();

                if (oldHeight <= gNodeHeight && gNodeHeight < upper) {
                    if (networkNode.isRoot() || isWithinChildBranch(parentNode, gNode, embedding, upper)) {
                        // update the node height relative to 'upper'
                        final double gNewNodeHeight = upper - (upper - gNodeHeight) * (upper - newHeight) / (upper - oldHeight);
                        gNode.setHeight(gNewNodeHeight);
                        m++;
                    }
                } else if (lower < gNodeHeight && gNodeHeight < oldHeight) {
                    if (isWithinChildBranch(networkNode, gNode, embedding, upper)) {
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
    private boolean isWithinChildBranch(NetworkNode snNode, Node gNode, IntegerParameter embedding, double upper) {
        final int traversalNodeNr = snNode.getTraversalNumber();
        Integer withinBrNr;  Node ancNode = gNode;
        do {
            withinBrNr = embedding.getMatrixValue(traversalNodeNr, ancNode.getNr());
            ancNode = ancNode.getParent();
        } while (withinBrNr < 0 && ancNode != null && ancNode.getHeight() < upper);

        return snNode.childBranchNumbers.contains(withinBrNr);
    }
}
