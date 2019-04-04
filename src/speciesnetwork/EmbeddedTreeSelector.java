package speciesnetwork;

import java.util.ArrayList;
import java.util.List;

import beast.core.BEASTInterface;
import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

public class EmbeddedTreeSelector extends CalculationNode implements EmbeddedTreeInterface {
	public final Input<List<BEASTInterface>> trees = new Input<>("tree", "EmbeddedTree-list to select from",
			new ArrayList<>(), Validate.REQUIRED, EmbeddedTreeInterface.class);
	public final Input<IntegerParameter> choice = new Input<>("choice", "EmbeddedTree-list to select from",
			Validate.REQUIRED);
	public final Input<Integer> index = new Input<>("index", "For re-use purposes, the relevant index of choice", 0);

	private EmbeddedTreeInterface tree() {
		return ((EmbeddedTreeInterface) (trees.get().get(choice.get().getNativeValue(index.get()))));
	}

	@Override
	public int getLeafNodeCount() {
		return tree().getLeafNodeCount();
	}

	@Override
	public int getInternalNodeCount() {
		return tree().getInternalNodeCount();
	}

	@Override
	public int getNodeCount() {
		return tree().getNodeCount();
	}

	@Override
	public Node getRoot() {
		return tree().getRoot();
	}

	@Override
	public Node getNode(int i) {
		return tree().getNode(i);
	}

	@Override
	public Node[] getNodesAsArray() {
		return tree().getNodesAsArray();
	}

	@Override
	public List<Node> getExternalNodes() {
		return tree().getExternalNodes();
	}

	@Override
	public List<Node> getInternalNodes() {
		return tree().getInternalNodes();
	}

	@Override
	public TaxonSet getTaxonset() {
		return tree().getTaxonset();
	}

	@Override
	public boolean somethingIsDirty() {
		return tree().somethingIsDirty();
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
	public void setMetaData(Node node, Double[] t, String pattern) {
		tree().setMetaData(node, t, pattern);
	}

	@Override
	public void initAndValidate() {
		trees.get().get(choice.get().getNativeValue(index.get())).initAndValidate();
	}

	@Override
	public void getMetaData(Node node, Double[] t, String pattern) {
		tree().getMetaData(node, t, pattern);
	}

	@Override
	public void makeCaterpillar(double rootHeight, double d, boolean b) {
		tree().makeCaterpillar(rootHeight, d, b);

	}

	@Override
	public void assignFromTree(Tree tempTree) {
		tree().assignFromTree(tempTree);
	}

}
