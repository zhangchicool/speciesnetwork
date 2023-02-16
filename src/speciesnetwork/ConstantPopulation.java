package speciesnetwork;

import java.text.DecimalFormat;
import java.util.List;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;

/**
 * @author Huw Ogilvie
 */

public class ConstantPopulation extends PopulationSizeModel {
    public final Input<RealParameter> popSizesInput =
            new Input<>("popSizes", "Constant per-branch population sizes.", Validate.REQUIRED);

    @Override
    public void initAndValidate() {
    }

    @Override
    public double branchLogP(int speciesBranchNumber, double[] perGenePloidy,
                             List<Double[]> branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts) {
        final RealParameter popSizes = popSizesInput.get();
        final double popSize = popSizes.getValue(speciesBranchNumber);

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

        return branchLogR - (branchQ * Math.log(popSize)) - (branchGamma / popSize);
    }

    @Override
    public void initPopSizes(int nPopulation) {
        final RealParameter popSizes = popSizesInput.get();
        popSizes.setDimension(nPopulation);
    }

    @Override
    public void initPopSizes(double popInitial) {
        final RealParameter popSizes = popSizesInput.get();

        for (int i = 0; i < popSizes.getDimension(); i++) {
            popSizes.setValue(i, popInitial);
        }
    }

    @Override
    public void serialize(NetworkNode speciesNetworkNode, StringBuilder buf, DecimalFormat df) {
        // do nothing
    }
}
