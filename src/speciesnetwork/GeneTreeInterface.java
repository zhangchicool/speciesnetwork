package speciesnetwork;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multiset;

import beast.core.Operator;

public interface GeneTreeInterface {

	public Embedding getEmbedding();

	public double logGammaSum();

	public ListMultimap<Integer, Double> coalescentTimes();

	public Multiset<Integer> coalescentLineageCounts();

	public void rebuildEmbedding(Operator operator);
}