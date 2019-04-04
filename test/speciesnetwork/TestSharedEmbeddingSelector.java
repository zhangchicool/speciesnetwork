package speciesnetwork;

import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import beast.core.parameter.IntegerParameter;
import beast.evolution.tree.Node;
import beast.util.TreeParser;
import speciesnetwork.EmbeddedTree;
import speciesnetwork.EmbeddedTreeInterface;
import speciesnetwork.EmbeddedTreeSelector;
import speciesnetwork.GeneTreeInSpeciesNetwork;
import speciesnetwork.operators.RebuildEmbedding;

public class TestSharedEmbeddingSelector extends ConstantPopulationTest {
	private int ntrees = 0; 
	private int counter = 0; 
	private IntegerParameter range;
	private List<EmbeddedTree> trees = new ArrayList<EmbeddedTree>();
	
	@Override
	protected EmbeddedTreeInterface treeFromRoot(Node root) {
        EmbeddedTreeSelector embeddedTree = new EmbeddedTreeSelector();
        trees.set(counter, new EmbeddedTree(root));
        embeddedTree.initByName(
        		"tree", trees,
        		"choice", range,
        		"index", counter
        		);
        ++counter;
		return embeddedTree;
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
