package speciesnetwork;

import beast.evolution.tree.TreeInterface;
public class TestEmbeddingSelectors extends ConstantPopulationTest {
	@Override
	protected GeneTreeSelector geneTree(TreeInterface tree, Embedding embedding) {
        GeneTreeInSpeciesNetwork geneTreeWrapper = new GeneTreeInSpeciesNetwork();
        geneTreeWrapper.initByName(
        		"geneTree", tree,
        		"embedding", embedding,
        		"ploidy", ploidy,
        		"speciesNetwork", speciesNetwork);
        GeneTreeSelector selector = new GeneTreeSelector();
        selector.initByName(
        		"tree", geneTreeWrapper,
        		"choice", "0"
        		);
		return selector;
	}
}
