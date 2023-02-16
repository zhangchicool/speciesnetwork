package snetworktests;

import java.util.ArrayList;
import java.util.List;
import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import beast.base.inference.State;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeParser;
import speciesnetwork.*;
import speciesnetwork.operators.RebuildEmbedding;

abstract class PopulationTestHelper {
    String newickSpeciesNetwork;
    List<String> newickGeneTrees = new ArrayList<>();

    TaxonSet speciesSuperset;
    TreeParser speciesTree;
    NetworkParser speciesNetwork;
    List<EmbeddedTree> geneTrees = new ArrayList<>();

    State state = null;
    MultispeciesCoalescent msc;
    int nSpecies;
    int nBranches;
    double popSize;
    double ploidy;
    double expectedLogP;
    List<int[]> embeddings;

    final double allowedError = 1e-6;

    abstract public TaxonSet generateSuperset();
    abstract public PopulationSizeModel generatePopulationModel();

    @Test
    public void testLogP() {
        speciesSuperset = generateSuperset();
        initializeSpeciesNetwork();
        initializeStateNodes();
        initializeGeneTrees(false);

        final PopulationSizeModel populationModel = generatePopulationModel();
        populationModel.initPopSizes(nBranches);
        populationModel.initPopSizes(popSize);

        msc = new MultispeciesCoalescent();
        msc.initByName("speciesNetwork", speciesNetwork, "geneTree", geneTrees, "populationModel", populationModel);

        double calculatedLogP = msc.calculateLogP();
        assertEquals(expectedLogP, calculatedLogP, allowedError);
    }

    private void initializeSpeciesNetwork() {
        speciesTree = new TreeParser();
        speciesTree.initByName("newick", newickSpeciesNetwork, "IsLabelledNewick", true, "adjustTipHeights", false);
        speciesNetwork = new NetworkParser();
        speciesNetwork.initByName("tree", speciesTree);
    }

    private void initializeStateNodes() {
        if (state == null) state = new State();
        state.initialise();
    }

    private void initializeGeneTrees(boolean reembed) {
        for (int i = 0; i < newickGeneTrees.size(); i++) {
            final String newick = newickGeneTrees.get(i);
            TreeParser treeParser = new TreeParser();
            treeParser.initByName("newick", newick, "IsLabelledNewick", true);
            EmbeddedTree embeddedTree = new EmbeddedTree(treeParser.getRoot());

            final int[] embedding = this.embeddings.get(i);
            final int nRow = treeParser.getNodeCount();
            final int nCol = embedding.length / nRow;
            embeddedTree.embedding.reset(nCol);
            for (int r = 0; r < nRow; r++) {
                for (int c = 0; c < nCol; c++)
                    embeddedTree.embedding.setDirection(r, c, embedding[r * nCol + c]);
            }

            geneTrees.add(embeddedTree);

        }
        if (reembed) { // rebuild the embedding
            RebuildEmbedding rebuildOperator = new RebuildEmbedding();
            rebuildOperator.initByName("speciesNetwork", speciesNetwork, "taxonset", speciesSuperset,
                    "geneTree", geneTrees);
            assertTrue(rebuildOperator.rebuildEmbedding());
        }
    }
}
