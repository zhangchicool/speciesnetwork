package speciesnetwork;

import org.junit.Test;

import beast.evolution.tree.Tree;

public class TestEmbedding {
	@Test
	public void testPossibleEmbedding() {
		Tree tree = new Tree("((A:2,B:2):1,C:3):1;");
		EmbeddedTree eTree = new EmbeddedTree(tree.getRoot());
		Network nw = new NetworkParser();
		nw.initByName(
				"tree", new Tree("((A:2,(B:1)#H1:1):1,(C:1,#H1:0):2):1;"));
		GeneTreeInSpeciesNetwork gTree = new GeneTreeInSpeciesNetwork();
		gTree.initByName(
				"speciesNetwork", nw,
				"geneTree", eTree,
				"ploidy", 1.0);
		gTree.getSpeciesOccupancy();
	}
	
	@Test
	public void testImpossibleEmbedding() {
		Tree tree = new Tree("((A:0.5,B:0.5):2.5,C:3):1;");
		EmbeddedTree eTree = new EmbeddedTree(tree.getRoot());
		Network nw = new NetworkParser();
		nw.initByName(
				"tree", new Tree("((A:2,(B:1)#H1:1):1,(C:1,#H1:0):2):1;"));
		GeneTreeInSpeciesNetwork gTree = new GeneTreeInSpeciesNetwork();
		gTree.initByName(
				"speciesNetwork", nw,
				"geneTree", eTree,
				"ploidy", 1.0);
		gTree.getSpeciesOccupancy();
	}
}
