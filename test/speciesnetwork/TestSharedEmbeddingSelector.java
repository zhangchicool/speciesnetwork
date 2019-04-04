package speciesnetwork;

import java.util.ArrayList;
import java.util.List;

import beast.core.parameter.IntegerParameter;
import beast.evolution.tree.TreeInterface;

public class TestSharedEmbeddingSelector extends ConstantPopulationTest {
	private int ntrees = 0; 
	private int counter = 0; 
	private IntegerParameter range;
	private List<GeneTreeInSpeciesNetwork> trees = new ArrayList<GeneTreeInSpeciesNetwork>();

	@Override
	protected GeneTreeSelector geneTree(TreeInterface tree, Embedding embedding) {
        GeneTreeInSpeciesNetwork geneTreeWrapper = new GeneTreeInSpeciesNetwork();
        geneTreeWrapper.initByName(
        		"geneTree", tree,
        		"embedding", embedding,
        		"ploidy", ploidy,
				"taxa", generateSuperset(),
				"speciesNetwork", speciesNetwork);
        trees.set(counter, geneTreeWrapper);
        GeneTreeSelector selector = new GeneTreeSelector();
        selector.initByName(
        		"tree", trees,
        		"choice", range,
        		"index", counter
        		);
        ++counter;
		return selector;
	}

	@Override
    protected void initializeGeneTrees(boolean reembed) {
		ntrees = newickGeneTrees.size();
		Integer[] rawRange = new Integer[ntrees];
		for (int i=0; i<ntrees; ++i) {
			rawRange[i] = i;
			trees.add(null);
		}
		range = new IntegerParameter(rawRange);
		
		super.initializeGeneTrees(reembed);
	}
}
