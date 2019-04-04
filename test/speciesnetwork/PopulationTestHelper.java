package speciesnetwork;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import beast.core.State;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import beast.util.TreeParser;
import speciesnetwork.operators.RebuildEmbedding;

abstract class PopulationTestHelper {
    String newickSpeciesNetwork;
    List<String> newickGeneTrees = new ArrayList<>();

    TaxonSet speciesSuperset;
    TreeParser speciesTree;
    NetworkParser speciesNetwork;
    List<GeneTreeInSpeciesNetwork> geneTreeWrappers = new ArrayList<>();

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
        msc.initByName("speciesNetwork", speciesNetwork, "geneTreeWithin", geneTreeWrappers, "populationModel", populationModel);

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

	protected TreeInterface treeFromRoot(Node root) {
		return new Tree(root);
	}

    protected void initializeGeneTrees(boolean reembed) {
		for (int i = 0; i < newickGeneTrees.size(); i++) {
            final String newick = newickGeneTrees.get(i);
            TreeParser treeParser = new TreeParser();
            treeParser.initByName("newick", newick, "IsLabelledNewick", true);
            TreeInterface embeddedTree = treeFromRoot(treeParser.getRoot());

            final int[] rawEmbedding = this.embeddings.get(i);
            final int nRow = treeParser.getNodeCount();
            final int nCol = rawEmbedding.length / nRow;
            Embedding embedding = new Embedding(nCol);
            for (int r = 0; r < nRow; r++) {
            	for (int c = 0; c < nCol; c++)
            		embedding.setDirection(r, c, rawEmbedding[r * nCol + c]);
            }

            GeneTreeInSpeciesNetwork geneTreeWrapper = new GeneTreeInSpeciesNetwork();
            geneTreeWrapper.initByName(
            		"geneTree", embeddedTree,
            		"embedding", embedding,
            		"ploidy", ploidy,
            		"speciesNetwork", speciesNetwork);
            geneTreeWrappers.add(geneTreeWrapper);
        }
        if (reembed) { // rebuild the embedding
            RebuildEmbedding rebuildOperator = new RebuildEmbedding();
            rebuildOperator.initByName("speciesNetwork", speciesNetwork, "taxonset", speciesSuperset,
                                       "geneTree", geneTreeWrappers);
            assertTrue(rebuildOperator.rebuildEmbedding());
        }
    }
}
