package sb2tests;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import beast.core.State;
import beast.core.parameter.RealParameter;
import starbeast2.ConstantPopulationIO;
import starbeast2.MultispeciesCoalescent;

public class ConstantIOTest extends PopulationTestHelper {
    private final String newickSpeciesTree = "((s0:0.32057156677143211,s3:0.32057156677143211):1.2653250035015629,(s1:0.56540722294658641,s2:0.56540722294658641):1.0204893473264085)";
    private final String newickGeneTreeA = "((((s0_tip1:0.3416660303037105,s3_tip0:0.3416660303037105):0.024561190897159135,s0_tip0:0.36622722120086965):0.0643095990846464,s3_tip1:0.43053682028551604):1.4201019862262891,((s1_tip0:0.14473698225381706,s1_tip1:0.14473698225381706):0.5135479407233198,(s2_tip0:0.19897724687831703,s2_tip1:0.19897724687831703):0.4593076760988198):1.1923538835346683)";
    private final String newickGeneTreeB = "(((s0_tip0:0.04173231934154758,s0_tip1:0.04173231934154758):0.7845256741090114,(s3_tip0:0.09482581277282173,s3_tip1:0.09482581277282173):0.7314321806777372):0.8703651781500925,((s1_tip0:0.33170960882423645,s1_tip1:0.33170960882423645):0.29497523293318856,(s2_tip0:0.2908611340994834,s2_tip1:0.2908611340994834):0.3358237076579416):1.0699383298432266)";

    private final double popSize = 0.3;
    private final double alpha = 1.5;
    private final double beta = 2.5;
    private final double expectedLogP = -14.5233984762; // this should be the right answer (calculated by hand)

    private RealParameter alphaParameter;
    private RealParameter betaParameter;

    MultispeciesCoalescent msc;

    @Test
    public void testLogP() throws Exception {
        nSpecies = 4;
        nBranches = (nSpecies * 2) - 1;
        popSizes = new double[nBranches];

        initializeTrees(newickSpeciesTree, newickGeneTreeA, newickGeneTreeB);

        alphaParameter = new RealParameter();
        betaParameter = new RealParameter();

        alphaParameter.initByName("value", String.valueOf(alpha));
        betaParameter.initByName("value", String.valueOf(beta));

        // Create dummy state to allow statenode editing
        State state = new State();
        state.initByName("stateNode", alphaParameter);
        state.initByName("stateNode", betaParameter);
        state.initialise();

        populationModel = new ConstantPopulationIO();
        populationModel.initByName("alpha", alphaParameter, "beta", betaParameter);
        populationModel.initPopSizes(nBranches);
        populationModel.initPopSizes(popSize);

        msc = new MultispeciesCoalescent();
        msc.initByName("tree", speciesTree, "geneTree", geneTreeList, "taxonSuperSet", speciesSuperSet, "populationModel", populationModel);

        double calculatedLogP = msc.calculateLogP();
        assertEquals(expectedLogP, calculatedLogP, allowedError);
    }
}
