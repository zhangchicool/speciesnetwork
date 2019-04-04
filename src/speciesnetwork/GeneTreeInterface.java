package speciesnetwork;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multiset;

import beast.evolution.tree.TreeInterface;

public interface GeneTreeInterface {
	public double getPloidy();
	public TreeInterface getTree();
	public Embedding getEmbedding();
	public void setEmbedding(Embedding newEmbedding);
	public double logGammaSum();
	public ListMultimap<Integer, Double> coalescentTimes();
	public Multiset<Integer> coalescentLineageCounts();
	public void rebuildEmbedding();
}