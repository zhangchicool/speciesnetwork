package speciesnetwork;

import java.util.*;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.core.Distribution;
import beast.math.distributions.Beta;
import speciesnetwork.utils.NodeHeightComparator;

/**
 * Birth-Death-Hybridization model for the species network.
 * @author Chi Zhang
 */

@Description("Birth-death-hybridization model")
public class BirthDeathHybridization extends Distribution {
    public final Input<Network> networkInput =
            new Input<>("network", "The species network.", Validate.REQUIRED);
    public final Input<RealParameter> birthRateInput =
            new Input<>("birthRate", "Speciation rate, lambda.");
    public final Input<RealParameter> deathRateInput =
            new Input<>("deathRate", "Extinction rate, mu.");
    public final Input<RealParameter> hybridRateInput =
            new Input<>("hybridRate", "Hybridization rate, nu.");
    public final Input<RealParameter> netDiversification =
            new Input<>("netDiversification", "Net diversification rate: lambda-mu-nu.");
    public final Input<RealParameter> turnOverInput =
            new Input<>("turnOver", "Turn over rate: (mu+nu)/lambda.");
    public final Input<RealParameter> hybridProportion =
            new Input<>("hybridProportion", "Hybridization proportion: nu/(mu+nu).");
    public final Input<RealParameter> rhoProbInput =
            new Input<>("rho", "Sampling prob. of extant species, rho.");
    public final Input<RealParameter> betaShapeInput =
            new Input<>("betaShape", "Shape of the symmetric beta prior on gamma probs (default is 1).");

    private static Comparator<NetworkNode> hc = new NodeHeightComparator();

    private double lambda, mu, nu;
    private Beta betaPrior;
    final static double EPSILON = 1e-8;

    @Override
    public void initAndValidate() {
        // set up alpha and beta parameters (alpha = beta)
        betaPrior = new Beta();
        final RealParameter betaShape;
        if (betaShapeInput.get() == null)
            betaShape = new RealParameter("1.0");  // default
        else
            betaShape = betaShapeInput.get();
        betaPrior.alphaInput.setValue(betaShape, betaPrior);
        betaPrior.betaInput.setValue(betaShape, betaPrior);
    }

    private void updateParameters() {
        if (birthRateInput.get() != null && deathRateInput.get() != null && hybridRateInput.get() != null) {
            lambda = birthRateInput.get().getValue();
            mu = deathRateInput.get().getValue();
            nu = hybridRateInput.get().getValue();
        }
        else if (netDiversification.get() != null && turnOverInput.get() != null && hybridProportion.get() != null) {
            lambda = netDiversification.get().getValue() / (1 - turnOverInput.get().getValue());
            mu = lambda * turnOverInput.get().getValue() * (1 - hybridProportion.get().getValue());
            nu = lambda * turnOverInput.get().getValue() * hybridProportion.get().getValue();
        }
        else {
            throw new IllegalArgumentException("Either specify birthRate & deathRate & hybridizationRate " +
                    "OR specify netDiversification & turnOver & hybridizationProportion.");
        }
        if (lambda <= 0.0 || mu <= 0.0 || nu <= 0.0) {
            throw new RuntimeException("Speciation, extinction and hybridization rates must be positive!");
        }

        // assuming complete sampling, rho is unused
        // rho = rhoProbInput.get() == null ? 1.0 : rhoProbInput.get().getValue();
    }

    @Override
    public double calculateLogP() {
        final Network network = networkInput.get();

        // sort the network nodes according to their heights in ascending order
        List<NetworkNode> nodes = Arrays.asList(network.getAllNodes());
        nodes.sort(hc);

        // get current values of lambda and nu
        updateParameters();

        logP = 0.0;
        // calculate probability of the network
        for (int i = 1; i < nodes.size(); i++) {
            final NetworkNode node = nodes.get(i);
            final double nodeHeight = node.getHeight();
            final double nextHeight = nodes.get(i-1).getHeight();

            if (nodeHeight > EPSILON) {  // rule out extant species
                // number of branches in time interval (nodeHeight, nextHeight)
                final int nBranch = network.getBranchCount((nodeHeight + nextHeight) /2.0);
                final double totalRate = nBranch * (lambda + mu) + nu * nBranch * (nBranch -1) /2;
                logP += totalRate * (nextHeight - nodeHeight);

                // rate for a particular event at time nodeHeight
                if (node.isReticulation()) {
                    logP += Math.log(nu);
                    logP += betaPrior.logDensity(node.inheritProb);
                }
                else if (node.isSpeciation()) {
                    logP += Math.log(lambda);
                }
                else if (node.isLeaf()) {
                    logP += Math.log(mu);
                }
            }
         }

        return logP;
    }

    @Override
    protected boolean requiresRecalculation() {
        return networkInput.get().isDirty() ||
                (birthRateInput.get() != null && birthRateInput.get().somethingIsDirty()) ||
                (deathRateInput.get() != null && deathRateInput.get().somethingIsDirty()) ||
                (hybridRateInput.get() != null && hybridRateInput.get().somethingIsDirty()) ||
                (netDiversification.get() != null && netDiversification.get().somethingIsDirty()) ||
                (turnOverInput.get() != null && turnOverInput.get().somethingIsDirty()) ||
                (hybridProportion.get() != null && hybridProportion.get().somethingIsDirty());
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) { }
}
