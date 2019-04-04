package speciesnetwork;

import beast.evolution.tree.Node;
import speciesnetwork.EmbeddedTree;
import speciesnetwork.EmbeddedTreeInterface;
import speciesnetwork.EmbeddedTreeSelector;

public class TestEmbeddingSelectors extends ConstantPopulationTest {
	@Override
	protected EmbeddedTreeInterface treeFromRoot(Node root) {
        EmbeddedTreeSelector embeddedTree = new EmbeddedTreeSelector();
        embeddedTree.initByName(
        		"tree", new EmbeddedTree(root),
        		"choice", "0"
        		);
		return embeddedTree;
	}
}
