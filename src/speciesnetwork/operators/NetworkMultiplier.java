package speciesnetwork.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;
import speciesnetwork.SanityChecks;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Chi Zhang
 */

@Description("Change network (and gene trees) internal node heights using a multiplier.")
public class NetworkMultiplier extends Operator {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);
    public final Input<List<Tree>> geneTreesInput =
            new Input<>("geneTree", "Gene tree within the species network.", new ArrayList<>());
    public final Input<Double> tuningInput =
            new Input<>("tuning", "A fine-tuning parameter (default is 0.3).", 0.3);

    private double tuning;

    @Override
    public void initAndValidate() {
        tuning = tuningInput.get();
    }

    @Override
    public double proposal() {
        final Network speciesNetwork = speciesNetworkInput.get();

        final double scaler = Math.exp (tuning * (Randomizer.nextDouble() - 0.5));

        // check proposed root height < origin height
        if (speciesNetwork.getRoot().getHeight() * scaler >= speciesNetwork.getOrigin().getHeight())
            return Double.NEGATIVE_INFINITY;

        // scale all internal network nodes
        speciesNetwork.startEditing(this);
        for (NetworkNode snNode : speciesNetwork.getInternalNodes()) {
            final double newHeight = scaler * snNode.getHeight();
            snNode.setHeight(newHeight);
        }
        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());
        double logProposalRatio =  Math.log(scaler) * speciesNetwork.getInternalNodes().length;

        // also scale all gene tree internal nodes
        for (Tree gTree : geneTreesInput.get()) {
            for (Node gNode : gTree.getInternalNodes()) {
                final double newHeight = scaler * gNode.getHeight();
                gNode.setHeight(newHeight);
            }
            logProposalRatio += Math.log(scaler) * gTree.getInternalNodeCount();
        }

        return logProposalRatio;
    }
}
