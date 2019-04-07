package speciesnetwork;

import java.util.ArrayList;
import java.util.List;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

public class TreeSelector extends CalculationNode implements TreeInterface {
	public Input<List<TreeInterface>> trees = new Input<>("tree", "The trees to select from", new ArrayList<>());
	public Input<IntegerParameter> index = new Input<>("index", "An integer parameter pointing to the tree currently chosen",
			Validate.REQUIRED);
	
	@Override
	public void initAndValidate() {
		if (index.get().lowerValueInput.get() != 0) {
			throw new RuntimeException("TreeSelector index lower must be 0");
		}
		if (index.get().upperValueInput.get() != trees.get().size() - 1) {
			throw new RuntimeException("TreeSelector index upper must match trees size");
		}
	}

	private TreeInterface tree() { return trees.get().get(index.get().getValue()); }
	
	@Override
	public int getLeafNodeCount() { return tree().getLeafNodeCount(); }

	@Override
	public int getInternalNodeCount() { return tree().getInternalNodeCount(); }

	@Override
	public int getNodeCount() { return tree().getNodeCount(); }

	@Override
	public Node getRoot() { return tree().getRoot(); }

	@Override
	public Node getNode(int i) { return tree().getNode(i); }

	@Override
	public Node[] getNodesAsArray() { return tree().getNodesAsArray(); }

	@Override
	public List<Node> getExternalNodes() { return tree().getExternalNodes(); }

	@Override
	public List<Node> getInternalNodes() { return tree().getInternalNodes(); }


	@Override
	public TaxonSet getTaxonset() { return tree().getTaxonset(); }


	@Override
	public void getMetaData(Node node, Double[] t, String pattern) { tree().getMetaData(node, t, pattern); }


	@Override
	public void setMetaData(Node node, Double[] t, String pattern) { tree().setMetaData(node, t, pattern); }

	@Override
	public boolean somethingIsDirty() {	return tree().somethingIsDirty() || index.get().somethingIsDirty() ; }

}
