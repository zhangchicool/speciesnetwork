package snetworktests;

import org.junit.Test;
import static org.junit.Assert.assertEquals;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.TreeParser;
import speciesnetwork.NetworkParser;
import speciesnetwork.BirthHybridizationModel;

public class BirthHybridizationTest {
    TreeParser treeParser;
    NetworkParser networkParser;

    final String newickNetwork = "(((A:0.2,(B:0.1)#H1:0.1)S1:0.3,(#H1:0.2,C:0.3)S2:0.2)R:0.1)";
    final RealParameter birthRate = new RealParameter("30");
    final RealParameter hybridRate = new RealParameter("20");

    // lam  <- 30
    // nu   <- 20
    // logp <- -lam * 0.1 + log(lam) +
    //         -(2*lam + nu) * 0.2 + log(lam) +
    //         -(3*lam + 3*nu) * 0.1 + log(lam) +
    //         -(4*lam + 6*nu) * 0.1 + log(nu) +
    //         -(3*lam + 3*nu) * 0.1
    final double expectedLogP = -59.80067558;
    final double allowedError = 1e-6;

    public BirthHybridizationTest() {
        treeParser = new TreeParser();
        networkParser = new NetworkParser();

        try {
            treeParser.initByName("newick", newickNetwork, "IsLabelledNewick", true, "adjustTipHeights", false);
            networkParser.initByName("tree", treeParser);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Test
    public void testLogP() {
        // System.out.println(networkParser.toString());
        BirthHybridizationModel networkPrior = new BirthHybridizationModel();
        networkPrior.initByName("network", networkParser, "birthRate", birthRate, "hybridRate", hybridRate);

        final double calculatedLogP = networkPrior.calculateLogP();
        assertEquals(expectedLogP, calculatedLogP, allowedError);
    }
}
