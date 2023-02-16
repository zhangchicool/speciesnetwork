package speciesnetwork.operators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Operator;
import beast.base.util.Randomizer;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;

@Description("Changes the value of gamma by applying a random walk to the logit of gamma.")
public class GammaProbRndWalk extends Operator {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);
    public final Input<Double> windowSizeInput =
            new Input<>("windowSize", "The size of the sliding window (default is 2).", 2.0);

    @Override
    public void initAndValidate() {
    }

    /**
     * y' = y + w*(u-0.5)
     * The proposal ratio on the transformed variable y is 1.
     * y = logit(r) = log(r/(1-r)) = log(r) - log(1-r)
     * dr/dy = d[exp(y)/(1+exp(y))] / dy = exp(y)/[(1+exp(y))^2]
     * So the proposal ratio for the original variable r is (dr'/dy') / (dr/dy) =
     * exp(w*(u-0.5)) * [(1+exp(y)) / (1+exp(y'))]^2
     */
    @Override
    public double proposal() {
        final Network speciesNetwork = speciesNetworkInput.get();
        final double windowSize = windowSizeInput.get();

        final int nReticulations = speciesNetwork.getReticulationNodeCount();
        if (nReticulations == 0)  // no reticulation
            return Double.NEGATIVE_INFINITY;

        speciesNetwork.startEditing(this);

        final int randomIndex = Randomizer.nextInt(nReticulations) + speciesNetwork.getReticulationOffset();
        final NetworkNode randomNode = speciesNetwork.getNode(randomIndex);

        final double currentGamma = randomNode.getGammaProb();                       // r
        final double currentLogOdds = Math.log(currentGamma / (1.0 - currentGamma)); // y
        final double logOddsShift = (Randomizer.nextDouble() - 0.5) * windowSize;    // w*(u-0.5)
        final double newLogOdds = currentLogOdds + logOddsShift;                     // y'

        final double logProposalRatio;
        final double newGamma;                                                       // r'
        if (newLogOdds >= 0.0) {  // to avoid numerical error, assuming logOddsShift is in (-1, 1)
            newGamma = 1.0 / (1.0 + Math.exp(-newLogOdds));
            logProposalRatio = logOddsShift +
                    2 * (Math.log(Math.exp(-currentLogOdds) + 1.0) - Math.log(Math.exp(-currentLogOdds) + Math.exp(logOddsShift)));
        } else {
            newGamma = Math.exp(newLogOdds) / (Math.exp(newLogOdds) + 1.0);
            logProposalRatio = logOddsShift +
                    2 * (Math.log(1.0 + Math.exp(currentLogOdds)) - Math.log(1.0 + Math.exp(newLogOdds)));
        }
        randomNode.setGammaProb(newGamma);

        return logProposalRatio;
    }
}
