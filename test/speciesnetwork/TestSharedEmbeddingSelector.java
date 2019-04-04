package speciesnetwork;

import java.util.ArrayList;
import java.util.List;

import beast.core.parameter.IntegerParameter;
import beast.evolution.tree.Node;

public class TestSharedEmbeddingSelector extends ConstantPopulationTest {
	private int ntrees = 0; 
	private int counter = 0; 
	private IntegerParameter range;
	private List<EmbeddedTree> trees = new ArrayList<EmbeddedTree>();
	
	@Override
	protected GeneTreeSelector treeFromRoot(Node root) {
		GeneTreeSelector embeddedTree = new GeneTreeSelector();
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
