package snetworktests;

import java.util.ArrayList;
import java.util.List;
import org.junit.Test;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import speciesnetwork.ConstantPopIntegrated;
import speciesnetwork.PopulationSizeModel;

public class ConstantPopIntegratedTest extends PopulationTestHelper {

    public ConstantPopIntegratedTest() {
        /*
        alpha  <- 5.0
        beta   <- 1.0
        gammaP <- 0.4
        logP <- 2*log(gammaP) + log(1-gammaP) + 2*log(1-gammaP) +
            log(1/2) + alpha*log(beta) - (alpha+1)*log(beta + 3*0.1/2 + 3*0.05/2 + 0.05/2) + log(gamma(alpha+1)/gamma(alpha)) +
            log(1/4) + alpha*log(beta) - (alpha+2)*log(beta + 0.07/2 + 0.1/2) + log(gamma(alpha+2)/gamma(alpha)) +
            log(1/2) + alpha*log(beta) - (alpha+1)*log(beta + 0.05/2 + 3*0.2/2) + log(gamma(alpha+1)/gamma(alpha)) +
            log(1/2) + alpha*log(beta) - (alpha+1)*log(beta + 3*0.05/2 + 0.25/2) + log(gamma(alpha+1)/gamma(alpha)) +
            alpha*log(beta) - alpha*log(beta + 0.2/2) + alpha*log(beta) - alpha*log(beta + 0.1/2) +
            log(1/4/8) + alpha*log(beta) - (alpha+5)*log(beta + 3*0.05/2 + 0.08/2 + 6*0.1/2 + 3*0.05/2 + 0.05/2) + log(gamma(alpha+5)/gamma(alpha))
         */
        expectedLogP = -2.010226796; // this should be the right answer (calculated by hand)

        nSpecies = 3;
        nBranches = 8;
        popSize = 0.1;
        ploidy = 2.0;

        newickSpeciesNetwork = "(((A:0.2,#H1[&gamma=0.4]:0.1)S1:0.3,((B:0.1)#H1:0.2,C:0.3)S2:0.2)R:0.1)";
        newickGeneTrees.add("(((a1:0.07,a2:0.07):0.48,(b1:0.25,b2:0.25):0.30):0.08,(b3:0.35,c1:0.35):0.28)");
        newickGeneTrees.add("((((a1:0.10,a2:0.10):0.50,(b1:0.05,b2:0.05):0.55):0.05,b3:0.65):0.05,c1:0.70)");
        
        embeddings = new ArrayList<>();
        final int[] embedding1 = {
                -1, -1, -1, -1,
                -1, -1, -1, -1,
                -1,  6, -1,  1,
                -1,  6, -1,  1,
                -1, -1,  7,  1,
                -1, -1,  2, -1,
                 4,  0, -1, -1,
                 4, -1, -1, -1,
                -1, -1, -1, -1,
                 5, -1, -1, -1,
                -1, -1, -1, -1};
        embeddings.add(embedding1);

        final int[] embedding2 = {
                -1, -1, -1, -1,
                -1, -1, -1, -1,
                -1, -1, -1, -1,
                -1, -1, -1, -1,
                 5, -1,  7,  1,
                 5, -1,  2, -1,
                 4,  0, -1, -1,
                 5, -1,  7,  1,
                -1, -1, -1, -1,
                -1, -1, -1, -1,
                -1, -1, -1, -1};
        embeddings.add(embedding2);
    }

    @Test
    public void testLogP() {
        super.testLogP();
    }

    @Override
    public TaxonSet generateSuperset() {
        List<Taxon> superSetList = new ArrayList<>();

        List<Taxon> taxonListA = new ArrayList<>();
        taxonListA.add(new Taxon("a1"));
        taxonListA.add(new Taxon("a2"));
        superSetList.add(new TaxonSet("A", taxonListA));

        List<Taxon> taxonListB = new ArrayList<>();
        taxonListB.add(new Taxon("b1"));
        taxonListB.add(new Taxon("b2"));
        taxonListB.add(new Taxon("b3"));
        superSetList.add(new TaxonSet("B", taxonListB));

        List<Taxon> taxonListC = new ArrayList<>();
        taxonListC.add(new Taxon("c1"));
        superSetList.add(new TaxonSet("C", taxonListC));

        return new TaxonSet(superSetList);
    }

    @Override
    public PopulationSizeModel generatePopulationModel() {
        final double alpha = 5.0, beta = 1.0;
        RealParameter alphaParameter = new RealParameter();
        RealParameter betaParameter = new RealParameter();
        alphaParameter.initByName("value", String.valueOf(alpha));
        betaParameter.initByName("value", String.valueOf(beta));

        state.initByName("stateNode", alphaParameter);
        state.initByName("stateNode", betaParameter);
        state.initialise();

        PopulationSizeModel populationModel = new ConstantPopIntegrated();
        populationModel.initByName("alpha", alphaParameter, "beta", betaParameter);
        
        return populationModel;
    }
}
