package speciesnetwork.operators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Operator;
import beast.base.util.Randomizer;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;
import speciesnetwork.SanityChecks;

/**
 * @author Chi Zhang
 */

@Description("Randomly select an internal network node and move its height uniformly.")
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
        // double oldHeight = pickedNode.getHeight();
        double newHeight = Randomizer.nextDouble() * (upper - lower) + lower;
        speciesNetwork.startEditing(this);
        pickedNode.setHeight(newHeight);
        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());

        return 0.0;
    }
}
