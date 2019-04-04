package speciesnetwork;

import java.util.Arrays;

import org.junit.Test;

import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Tree;

public class TestEmbedding {
	@Test
	public void testTrivialEmbedding() {
		Tree tree = new Tree("(a:2,b:2);");
		Network nw = new NetworkParser();
		nw.initByName(
				"tree", new Tree("(A:2,B:2);"));
		GeneTreeInSpeciesNetwork gTree = new GeneTreeInSpeciesNetwork();
		Taxon[] t1 = { new Taxon("a") };
		Taxon[] t2 = { new Taxon("b") };
		TaxonSet[] tx = {
				new TaxonSet("A", Arrays.asList(t1)),
				new TaxonSet("B", Arrays.asList(t2))			
				};
		TaxonSet taxa = new TaxonSet(Arrays.asList(tx));
		gTree.initByName(
				"speciesNetwork", nw,
				"geneTree", tree,
				"embedding", new Embedding(3, 3),
				"taxa", taxa,
				"ploidy", 1.0);
		gTree.rebuildEmbedding();
		gTree.getSpeciesOccupancy();
	}
	
	@Test
	public void testPossibleEmbedding() {
		Tree tree = new Tree("((a:2,b:2):1,c:3):1;");
		Network nw = new NetworkParser();
		nw.initByName(
				"tree", new Tree("((A:2,(B:1)#H1:1):1,(C:1,#H1:0):2):1;"));
		Taxon[] t1 = { new Taxon("a") };
		Taxon[] t2 = { new Taxon("b") };
		Taxon[] t3 = { new Taxon("c") };
		TaxonSet[] tx = {
				new TaxonSet("A", Arrays.asList(t1)),
				new TaxonSet("B", Arrays.asList(t2)),			
				new TaxonSet("C", Arrays.asList(t3))			
				};
		TaxonSet taxa = new TaxonSet(Arrays.asList(tx));
		GeneTreeInSpeciesNetwork gTree = new GeneTreeInSpeciesNetwork();
		gTree.initByName(
				"speciesNetwork", nw,
				"geneTree", tree,
				"embedding", new Embedding(5, 6),
				"taxa", taxa,
				"ploidy", 1.0);
		gTree.rebuildEmbedding();
		gTree.getSpeciesOccupancy();
	}
	
	@Test(expected = RuntimeException.class)
	public void testImpossibleEmbedding() {
		Tree tree = new Tree("((a:0.5,b:0.5):2.5,c:3):1;");
		Network nw = new NetworkParser();
		nw.initByName(
				"tree", new Tree("((A:2,(B:1)#H1:1):1,(C:1,#H1:0):2):1;"));
		Taxon[] t1 = { new Taxon("a") };
		Taxon[] t2 = { new Taxon("b") };
		Taxon[] t3 = { new Taxon("c") };
		TaxonSet[] tx = {
				new TaxonSet("A", Arrays.asList(t1)),
				new TaxonSet("B", Arrays.asList(t2)),			
				new TaxonSet("C", Arrays.asList(t3))			
				};
		TaxonSet taxa = new TaxonSet(Arrays.asList(tx));
		GeneTreeInSpeciesNetwork gTree = new GeneTreeInSpeciesNetwork();
		gTree.initByName(
				"speciesNetwork", nw,
				"geneTree", tree,
				"embedding", new Embedding(5, 6),
				"taxa", taxa,
				"ploidy", 1.0);
		gTree.rebuildEmbedding();
		gTree.getSpeciesOccupancy();
	}
}
