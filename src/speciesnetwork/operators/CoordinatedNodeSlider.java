package speciesnetwork.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;
import speciesnetwork.SanityChecks;

/**
 * @author Chi Zhang
 */

@Description("Randomly select an internal network node and move its height using a sliding window." +
        "Also update the affected gene tree node heights at each locus to maintain the compatibility.")
public class CoordinatedNodeSlider extends CoordinatedOperator {
    public final Input<RealParameter> originInput =
            new Input<>("origin", "The time when the process started.", Validate.REQUIRED);
    public final Input<Double> windowSizeInput =
            new Input<>("windowSize", "Window size of the uniform proposal (default is 0.02).", 0.02);
    public final Input<Boolean> isNormalInput =
            new Input<>("isNormal", "Using normal proposal (default: uniform proposal).", false);
    public final Input<Double> sigmaInput =
            new Input<>("sigma", "Standard deviation of the normal proposal (default is 0.01).", 0.01);

    @Override
    public double proposal() {
        final Network speciesNetwork = speciesNetworkInput.get();

        // pick an internal node randomly, including origin
        final NetworkNode[] internalNodes = speciesNetwork.getInternalNodesWithOrigin();
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

        // propose a new height, reflect it back if it's outside the boundary
        double oldHeight = pickedNode.getHeight();
        double newHeight;
        if (isNormalInput.get()) {
            final double sigma = sigmaInput.get();
            newHeight = oldHeight + Randomizer.nextGaussian() * sigma;
        } else {
            final double windowSize = windowSizeInput.get();
            newHeight = oldHeight + (Randomizer.nextDouble() - 0.5) * windowSize;
        }
        while (newHeight < lower || newHeight > upper) {
            if (newHeight < lower)
                newHeight = 2.0 * lower - newHeight;
            if (newHeight > upper)
                newHeight = 2.0 * upper - newHeight;
        }

        // update the new node height
        if (pickedNode.isOrigin()) {
            final RealParameter originTime = originInput.get();
            if (outsideBounds(newHeight, originTime))
                return Double.NEGATIVE_INFINITY;

            originTime.setValue(newHeight);
        }
        speciesNetwork.startEditing(this);
        pickedNode.setHeight(newHeight);
        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());

        return updateRubberBand(pickedNode, oldHeight, newHeight, lower, upper);
    }

    private boolean outsideBounds(final double value, final RealParameter param) {
        final Double l = param.getLower();
        final Double h = param.getUpper();

        return (value < l || value > h);
    }
}
