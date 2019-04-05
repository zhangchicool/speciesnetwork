package speciesnetwork;

import java.util.Arrays;

import org.junit.Assert;
import org.junit.Test;

import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Tree;

public class TestEmbedding {
	@Test
	public void testTrivialEmbedding() {
		// Test whether a trivial embedding of two individuals in two taxa, exactly
		// following them, gives valid species network occupancies
		Tree tree = new Tree("(a:2,b:2);");
		Network nw = new NetworkParser();
		nw.initByName("tree", new Tree("(A:2,B:2);"));
		GeneTreeInSpeciesNetwork gTree = new GeneTreeInSpeciesNetwork();
		TaxonSet[] tx = { new TaxonSet("A", Arrays.asList(new Taxon[] { new Taxon("a") })),
				new TaxonSet("B", Arrays.asList(new Taxon[] { new Taxon("b") })) };
		TaxonSet taxa = new TaxonSet(Arrays.asList(tx));
		gTree.initByName("speciesNetwork", nw, "geneTree", tree, "taxa", taxa, "ploidy", 1.0);
		gTree.rebuildEmbedding();
		double[][] x = gTree.getSpeciesOccupancy();
		Assert.assertArrayEquals(new double[] { 2.0, 0.0, 0.0 }, x[0], 1e-5);
		Assert.assertArrayEquals(new double[] { 0.0, 2.0, 0.0 }, x[1], 1e-5);
		Assert.assertArrayEquals(new double[] { 0.0, 0.0, Double.POSITIVE_INFINITY }, x[2], 1e-5);
	}

	@Test
	public void testPossibleEmbedding() {
		// Test whether a possible embedding into a network, exactly
		// following its branches, gives valid species network occupancies
		Tree tree = new Tree("((a:2,b:2):1,c:3):1;");
		Network nw = new NetworkParser();
		nw.initByName("tree", new Tree("((A:2,(B:1)#H1:1):1,(C:1.5,#H1:0.5):1.5):1;"));
		Taxon[] t1 = { new Taxon("a") };
		Taxon[] t2 = { new Taxon("b") };
		Taxon[] t3 = { new Taxon("c") };
		TaxonSet[] tx = { new TaxonSet("A", Arrays.asList(t1)), new TaxonSet("B", Arrays.asList(t2)),
				new TaxonSet("C", Arrays.asList(t3)) };
		TaxonSet taxa = new TaxonSet(Arrays.asList(tx));
		GeneTreeInSpeciesNetwork gTree = new GeneTreeInSpeciesNetwork();
		gTree.initByName("speciesNetwork", nw, "geneTree", tree, "taxa", taxa, "ploidy", 1.0);
		gTree.rebuildEmbedding();
		double[][] x = gTree.getSpeciesOccupancy();
		Assert.assertArrayEquals(new double[] { 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, x[0], 1e-5);
		Assert.assertArrayEquals(new double[] { 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 }, x[1], 1e-5);
		Assert.assertArrayEquals(new double[] { 0.0, 0.0, 1.5, 0.0, 0.0, 1.5, 0.0, 0.0 }, x[2], 1e-5);
		Assert.assertArrayEquals(new double[] { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 }, x[3], 1e-5);
		Assert.assertArrayEquals(new double[] { 0.0, 0.0, 0.0, Double.POSITIVE_INFINITY, 0.0, 0.0, 0.0, 0.0 }, x[4],
				1e-5);
	}

	@Test
	public void testOtherPossibleEmbedding() {
		// Test whether a possible embedding into a network, using a hybridization edge
		// and with gene tree splits before species tree splits, gives valid
		// species network occupancies
		Tree tree = new Tree("(a:4.0,(b:2.0,c:2.0):2.0):1;");
		Network nw = new NetworkParser();
		nw.initByName("tree", new Tree("((A:2,(B:1)#H1:1):1,(C:1.5,#H1:0.5):1.5):1;"));
		Taxon[] t1 = { new Taxon("a") };
		Taxon[] t2 = { new Taxon("b") };
		Taxon[] t3 = { new Taxon("c") };
		TaxonSet[] tx = { new TaxonSet("A", Arrays.asList(t1)), new TaxonSet("B", Arrays.asList(t2)),
				new TaxonSet("C", Arrays.asList(t3)) };
		TaxonSet taxa = new TaxonSet(Arrays.asList(tx));
		GeneTreeInSpeciesNetwork gTree = new GeneTreeInSpeciesNetwork();
		gTree.initByName("speciesNetwork", nw, "geneTree", tree, "taxa", taxa, "ploidy", 1.0);
		gTree.rebuildEmbedding();
		double[][] x = gTree.getSpeciesOccupancy();
		Assert.assertArrayEquals(new double[] { 2.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0 }, x[0], 1e-5);
		Assert.assertArrayEquals(new double[] { 0.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5 }, x[1], 1e-5);
		Assert.assertArrayEquals(new double[] { 0.0, 0.0, 1.5, 0.0, 0.0, 0.5, 0.0, 0.0 }, x[2], 1e-5);
		Assert.assertArrayEquals(new double[] { 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0 }, x[3], 1e-5);
		Assert.assertArrayEquals(new double[] { 0.0, 0.0, 0.0, Double.POSITIVE_INFINITY, 0.0, 0.0, 0.0, 0.0 }, x[4],
				1e-5);
	}

	@Test(expected = RuntimeException.class)
	public void testImpossibleEmbedding() {
		// Test whether an impossible embedding fails: The most recent common ancestor
		// of a and b is at 0.5, but the lineages to which a and b belong (A and B,
		// respectively) have already split by then. This should raise an exception.
		Tree tree = new Tree("((a:0.5,b:0.5):2.5,c:3):1;");
		Network nw = new NetworkParser();
		nw.initByName("tree", new Tree("((A:2,(B:1)#H1:1):1,(C:1.5,#H1:0.5):1.5):1;"));
		Taxon[] t1 = { new Taxon("a") };
		Taxon[] t2 = { new Taxon("b") };
		Taxon[] t3 = { new Taxon("c") };
		TaxonSet[] tx = { new TaxonSet("A", Arrays.asList(t1)), new TaxonSet("B", Arrays.asList(t2)),
				new TaxonSet("C", Arrays.asList(t3)) };
		TaxonSet taxa = new TaxonSet(Arrays.asList(tx));
		GeneTreeInSpeciesNetwork gTree = new GeneTreeInSpeciesNetwork();
		gTree.initByName("speciesNetwork", nw, "geneTree", tree, "taxa", taxa, "ploidy", 1.0);
		gTree.rebuildEmbedding();
		gTree.getSpeciesOccupancy();
	}
}
