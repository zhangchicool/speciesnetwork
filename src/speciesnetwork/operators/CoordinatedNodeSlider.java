package speciesnetwork.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import speciesnetwork.NetworkNode;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Chi Zhang
 */

@Description("Randomly select an internal network node and move its height using a sliding window." +
        "Also update the affected gene tree node heights at each locus to maintain the compatibility.")
public class CoordinatedNodeSlider extends NodeSlider {
    public final Input<List<Tree>> geneTreesInput =
            new Input<>("geneTree", "Gene tree within the species network.", new ArrayList<>());
    public final Input<List<IntegerParameter>> embeddingsInput =
            new Input<>("embedding", "The embedding matrix.", new ArrayList<>());

    private double nLoci;

    @Override
    public void initAndValidate() {
        nLoci = geneTreesInput.get().size();
        if (embeddingsInput.get().size() != nLoci) {
            throw new RuntimeException("Number of embeddings doesn't match number of gene trees!");
        }
    }

    @Override
    public double proposal() {
        double logProposalRatio = super.proposal();

        // update the gene tree node heights accordingly
        logProposalRatio += updateRubberBand(snNode);

        return logProposalRatio;
    }

    /**
     * The RubberBand algorithm of Rannala & Yang, 2003 (Appendix Step 4)
     */
    private double updateRubberBand(NetworkNode snNode) {
        if (!snNode.isSpeciation())  // necessary only when speciation node
            return 0.0;

        double m = 0;  // # gene node heights changed relative to 'upper'
        double n = 0;  // # gene node heights changed relative to 'lower'
        final List<Tree> geneTrees = geneTreesInput.get();
        final List<IntegerParameter> embeddings = embeddingsInput.get();

        for (int i = 0; i < nLoci; i++) {  // loop over all loci
            Tree gTree = geneTrees.get(i);
            IntegerParameter embedding = embeddings.get(i);

            for (Node gNode : gTree.getInternalNodes()) {  // loop over each gene node
                final double gNodeHeight = gNode.getHeight();

                Integer snNextBrNr;
                if (lower < gNodeHeight && gNodeHeight < oldHeight) {
                    final int traversalNodeNr = snNode.getTraversalNumber();
                    Node ancNode = gNode;
                    do {
                        snNextBrNr = embedding.getMatrixValue(traversalNodeNr, ancNode.getNr());
                        ancNode = ancNode.getParent();
                    } while (snNextBrNr < 0 && ancNode.getHeight() < oldHeight);
                    // check if gNode is in populations represented by snNode's children
                    if (snNode.childBranchNumbers.contains(snNextBrNr)) {
                        // update the node height relative to 'lower'
                        final double gHeightNew = lower + (gNodeHeight - lower) * (newHeight - lower) / (oldHeight - lower);
                        gNode.setHeight(gHeightNew);
                        n++;
                    }
                }
                else if (oldHeight <= gNodeHeight && snNode.isRoot()) {
                    // also update the node height relative to 'lower' (as upper limit is infinity)
                    final double gHeightNew = lower + (gNodeHeight - lower) * (newHeight - lower) / (oldHeight - lower);
                    gNode.setHeight(gHeightNew);
                    n++;
                }
                else if (oldHeight <= gNodeHeight && gNodeHeight < upper) {
                    final Integer snBranchNr = snNode.gammaBranchNumber;
                    final int traversalNodeNr = snNode.getParentByBranch(snBranchNr).getTraversalNumber();
                    Node ancNode = gNode;
                    do {
                        snNextBrNr = embedding.getMatrixValue(traversalNodeNr, ancNode.getNr());
                        ancNode = ancNode.getParent();
                    } while (snNextBrNr < 0 && ancNode.getHeight() < upper);
                    // check if gNode is in the population represented by snNode
                    if (snBranchNr.equals(snNextBrNr)) {
                        // update the node height relative to 'upper'
                        final double gHeightNew = upper - (upper - gNodeHeight) * (upper - newHeight) / (upper - oldHeight);
                        gNode.setHeight(gHeightNew);
                        m++;
                    }
                }
            }
        }

        return m * Math.log((upper - newHeight)/(upper - oldHeight)) +
               n * Math.log((newHeight - lower)/(oldHeight - lower));
    }
}
