package speciesnetwork;

import beast.evolution.tree.Node;

public class TestEmbeddingSelectors extends ConstantPopulationTest {
	@Override
	protected GeneTreeSelector treeFromRoot(Node root) {
        GeneTreeSelector embeddedTree = new GeneTreeSelector();
        embeddedTree.initByName(
        		"tree", new EmbeddedTree(root),
        		"choice", "0"
        		);
		return embeddedTree;
	}
}
