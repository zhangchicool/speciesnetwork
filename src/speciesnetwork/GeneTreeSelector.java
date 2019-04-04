package speciesnetwork;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multiset;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.evolution.tree.TreeInterface;

public class GeneTreeSelector extends CalculationNode implements GeneTreeInterface {
	public final Input<List<GeneTreeInSpeciesNetwork>> trees = new Input<>("tree", "EmbeddedTree-list to select from",
			new ArrayList<>(), Validate.REQUIRED);
	public final Input<IntegerParameter> choice = new Input<>("choice", "EmbeddedTree-list to select from",
			Validate.REQUIRED);
	public final Input<Integer> index = new Input<>("index", "For re-use purposes, the relevant index of choice", 0);

	private GeneTreeInSpeciesNetwork tree() {
		return trees.get().get(choice.get().getNativeValue(index.get()));
	}

	@Override
	public TreeInterface getTree() {
		// TODO Auto-generated method stub
		return tree().getTree();
	}

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public Embedding getEmbedding() {
		return tree().getEmbedding();
	}

	@Override
	public void setEmbedding(Embedding newEmbedding) {
		tree().setEmbedding(newEmbedding);
	}

	@Override
	public double getPloidy() {
		return tree().getPloidy();
	}

	@Override
	public double logGammaSum() {
		return tree().logGammaSum();
	}

	@Override
	public ListMultimap<Integer, Double> coalescentTimes() {
		return tree().coalescentTimes();
	}

	@Override
	public Multiset<Integer> coalescentLineageCounts() {
		return tree().coalescentLineageCounts();
	}
}
