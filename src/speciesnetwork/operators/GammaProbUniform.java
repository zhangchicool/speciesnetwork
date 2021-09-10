package speciesnetwork.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.util.Randomizer;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;

@Description("Changes the value of gamma by uniformly selecting a value in its range.")
public class GammaProbUniform extends Operator {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        final Network speciesNetwork = speciesNetworkInput.get();

        final int nReticulations = speciesNetwork.getReticulationNodeCount();
        if (nReticulations == 0)  // no reticulation
            return Double.NEGATIVE_INFINITY;

        speciesNetwork.startEditing(this);

        final int randomIndex = Randomizer.nextInt(nReticulations) + speciesNetwork.getReticulationOffset();
        final NetworkNode randomNode = speciesNetwork.getNode(randomIndex);

        final double newGamma = Randomizer.nextDouble();
        randomNode.setGammaProb(newGamma);

        return 0.0;
    }
}
