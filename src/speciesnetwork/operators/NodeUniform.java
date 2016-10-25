package speciesnetwork.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;
import speciesnetwork.SanityChecks;

/**
 * Randomly pick an internal network node.
 * Change its height uniformly between the lower and upper limit.
 * Time of origin is not changed.
 *
 * @author Chi Zhang
 */

@Description("Randomly selects an internal network node and move its height uniformly.")
public class NodeUniform extends Operator {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        final Network speciesNetwork = speciesNetworkInput.get();

        // pick an internal node randomly
        final NetworkNode[] internalNodes = speciesNetwork.getInternalNodes();
        final int randomIndex = Randomizer.nextInt(internalNodes.length);
        NetworkNode snNode = internalNodes[randomIndex];

        // determine the lower and upper bounds
        double upper = Double.MAX_VALUE;
        for (NetworkNode p: snNode.getParents()) {
            upper = Math.min(upper, p.getHeight());
        }

        double lower = 0.0;
        for (NetworkNode c: snNode.getChildren()) {
            lower = Math.max(lower, c.getHeight());
        }

        // propose a new height uniformly
        final double newHeight = Randomizer.nextDouble() * (upper - lower) + lower;
        speciesNetwork.startEditing(this);
        snNode.setHeight(newHeight);
        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());

        return 0.0;
    }
}
