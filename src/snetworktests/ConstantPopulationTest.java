package snetworktests;

import java.util.ArrayList;
import java.util.List;
import org.junit.Test;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import speciesnetwork.ConstantPopulation;
import speciesnetwork.PopulationSizeModel;

public class ConstantPopulationTest extends PopulationTestHelper {

    public ConstantPopulationTest() {
        /*
        gammaP <- 0.4
        logGene1 <- - log(2*0.1) - 0.07/(2*0.1) - 0.1*3/(2*0.1) +
            2*log(1-gammaP) + log(gammaP) - 0.1/(2*0.1) +
            - log(2*0.1) - 0.05*3/(2*0.1) - 0.25/(2*0.1) - log(2*0.1) - 0.05/(2*0.1) +
            - 2*log(2*0.1) - 0.05*3/(2*0.1) - 0.08/(2*0.1)
        logGene2 <- - log(2*0.1) - 0.1/(2*0.1) +
            - log(2*0.1) - 0.05*3/(2*0.1) - 0.05/(2*0.1) +
            2*log(gammaP) - 0.2/(2*0.1) - 0.2*3/(2*0.1) +
            - 3*log(2*0.1) - 0.1*6/(2*0.1) - 0.05*3/(2*0.1) - 0.05/(2*0.1)
        logP <- logGene1 + logGene2
        */
        expectedLogP = -2.52067921; // -0.046217525 -2.474461685

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
        RealParameter popSizesParameter = new RealParameter();
        popSizesParameter.initByName("value", String.valueOf(popSize));

        state.initByName("stateNode", popSizesParameter);
        state.initialise();

        PopulationSizeModel populationModel = new ConstantPopulation();
        populationModel.initByName("popSizes", popSizesParameter);
        
        return populationModel;
    }
}
