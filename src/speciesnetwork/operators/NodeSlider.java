package speciesnetwork.operators;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.util.Randomizer;

import speciesnetwork.Network;
import speciesnetwork.NetworkNode;

/**
 * @author Chi Zhang
 */

@Description("Randomly selects an internal network node and move its height using an uniform sliding window.")
public class NodeSlider extends Operator {
    public Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);
    public Input<List<RebuildEmbedding>> rebuildEmbeddingsInput = new Input<>("rebuildEmbedding",
            "Operator which rebuilds embedding of gene tree within species network.", new ArrayList<>());
    public Input<Double> windowSizeInput =
            new Input<>("windowSize", "The size of the sliding window, default 0.01.", 0.01);

    // empty constructor to facilitate construction by XML + initAndValidate
    public NodeSlider() {
    }

    @Override
    public void initAndValidate() {
    }

    /**
     * Propose a new network-node height from a uniform distribution.
     * If the new value is outside the boundary, the excess is reflected back into the interval.
     * The proposal ratio of this slider move is 1.0.
     * Then rebuild the embedding of the gene trees.
     */
    @Override
    public double proposal() {
        List<RebuildEmbedding> reembedOps = rebuildEmbeddingsInput.get();
        final int nLoci = reembedOps.size();

        // check the embedding in the current species network
        int oldChoices = 0;
        for (int i = 0; i < nLoci; i++) {
            RebuildEmbedding rebuildOperator = reembedOps.get(i);
            final int nChoices = rebuildOperator.rebuildEmbedding(false);
            if(nChoices < 0)
                throw new RuntimeException("Developer ERROR: current embedding invalid! geneTree " + i);
            oldChoices += nChoices;
        }

        // pick an internal node randomly
        Network speciesNetwork = speciesNetworkInput.get();
        List<NetworkNode> intNodes = speciesNetwork.getInternalNodes();
        NetworkNode snNode = intNodes.get(Randomizer.nextInt(intNodes.size()));

        // determine the lower and upper bounds
        NetworkNode leftParent = snNode.getLeftParent();
        NetworkNode rightParent = snNode.getRightParent();
        double upper = Double.MAX_VALUE;
        if (leftParent != null)
            upper = leftParent.getHeight();
        if (rightParent != null && upper > rightParent.getHeight())
            upper = rightParent.getHeight();
        NetworkNode leftChild = snNode.getLeftChild();
        NetworkNode rightChild = snNode.getRightChild();
        double lower = Double.MIN_NORMAL;
        if (leftChild != null)
            lower = leftChild.getHeight();
        if (rightChild != null && lower < rightChild.getHeight())
            lower = rightChild.getHeight();
        if (lower >= upper)
            throw new RuntimeException("Developer ERROR: upper bound must be larger than lower bound!");

        // propose a new height, reflect it back if it's outside the boundary
        final double windowSize = windowSizeInput.get();
        final double oldHeight = snNode.getHeight();
        double newHeight = oldHeight + (Randomizer.nextDouble() - 0.5) * windowSize;
        while (newHeight < lower || newHeight > upper) {
            if (newHeight < lower)
                newHeight = 2.0 * lower - newHeight;
            if (newHeight > upper)
                newHeight = 2.0 * upper - newHeight;
        }

        // update the new node height
        snNode.setHeight(newHeight);

        // update the embedding in the new species network
        int newChoices = 0;
        for (int i = 0; i < nLoci; i++) {
            RebuildEmbedding rebuildOperator = reembedOps.get(i);
            // rebuild the embedding
            final int nChoices = rebuildOperator.rebuildEmbedding(true);
            if(nChoices < 0)
                return Double.NEGATIVE_INFINITY;
            newChoices += nChoices;
        }

        return (newChoices - oldChoices) * Math.log(2);
    }
}
