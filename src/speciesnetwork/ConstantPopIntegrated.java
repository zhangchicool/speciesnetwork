package speciesnetwork;

import java.text.DecimalFormat;
import java.util.List;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;

/**
 * @author Huw Ogilvie
 */

public class ConstantPopIntegrated extends PopulationSizeModel {
    public final Input<RealParameter> invgammaShapeInput = new Input<>("alpha",
            "Shape of the inverse gamma distribution on population sizes.", Validate.REQUIRED);
    public final Input<RealParameter> invgammaScaleInput = new Input<>("beta",
            "Scale of the inverse gamma distribution on population sizes.", Validate.REQUIRED);
    public final Input<RealParameter> invgammaMeanInput = new Input<>("mean",
            "Mean of the inverse gamma distribution on population sizes.", Validate.XOR, invgammaScaleInput);

    @Override
    public void initAndValidate() {
    }

    @Override
    public double branchLogP(int speciesBranchNumber, double[] perGenePloidy,
                             List<Double[]> branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts) {
        final RealParameter invgammaShape = invgammaShapeInput.get();
        final RealParameter invgammaScale = invgammaScaleInput.get();
        final RealParameter invgammaMean = invgammaMeanInput.get();
        final double alpha = invgammaShape.getValue();
        final double beta;
        if (invgammaScale != null)
            beta = invgammaScale.getValue();
        else
            beta = invgammaMean.getValue() * (alpha - 1.0);
        final int nGenes = perGenePloidy.length;

        int branchQ = 0;
        double branchLogR = 0.0;
        double branchGamma = 0.0;

        for (int j = 0; j < nGenes; j++) {
            final int geneN = branchLineageCounts[j];
            final Double[] geneCoalescentTimes = branchCoalescentTimes.get(j);
            final int geneK = branchEventCounts[j];
            final double genePloidy = perGenePloidy[j]; 
            branchLogR -= geneK * Math.log(genePloidy);
            branchQ += geneK;

            double partialGamma = 0.0;
            for (int i = 0; i < geneK; i++) {
                partialGamma += (geneCoalescentTimes[i + 1] - geneCoalescentTimes[i])
                                * (geneN - i) * (geneN - i - 1.0) / 2.0;
            }
            if (geneN - geneK > 1) {
                partialGamma += (geneCoalescentTimes[geneK + 1] - geneCoalescentTimes[geneK])
                                * (geneN - geneK) * (geneN - geneK - 1.0) / 2.0;
            }
            branchGamma += partialGamma / genePloidy;
        }

        double logGammaRatio = 0.0;
        for (int i = 0; i < branchQ; i++) {
            logGammaRatio += Math.log(alpha + i);
        }

        return branchLogR + alpha * Math.log(beta) - (alpha + branchQ) * Math.log(beta + branchGamma) + logGammaRatio;
    }

    @Override
    public void initPopSizes(int nPopulation) {
        // do nothing
    }

    @Override
    public void initPopSizes(double popInitial) {
        // do nothing
    }

    @Override
    public void serialize(NetworkNode speciesNetworkNode, StringBuilder buf, DecimalFormat df) {
        // do nothing
    }

    /* private void debug(int speciesBranchNumber, double[] perGenePloidy,
                       List<Double[]> branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts) {
        for (int i = 0; i < branchCoalescentTimes.size(); i++) {
            StringBuffer sb = new StringBuffer();
            sb.append(speciesBranchNumber);
            sb.append(".");
            sb.append(i);
            sb.append(": ");
            sb.append(branchLineageCounts[i]);
            sb.append(" - ");
            sb.append(branchEventCounts[i]);
            sb.append(" | ");
            for (int j = 0; j < branchCoalescentTimes.get(i).length; j++) {
                sb.append(branchCoalescentTimes.get(i)[j]);
                sb.append(", ");
            }
            System.out.println(sb.toString());
        }
    } */
}
