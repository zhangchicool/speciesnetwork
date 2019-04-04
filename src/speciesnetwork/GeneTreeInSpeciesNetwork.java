package speciesnetwork;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multiset;
import com.google.common.collect.HashMultiset;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.util.Randomizer;

/**
 * @author Huw Ogilvie
 * @author Chi Zhang
 */

public class GeneTreeInSpeciesNetwork extends CalculationNode implements GeneTreeInterface {
	public final Input<Network> speciesNetworkInput = new Input<>("speciesNetwork",
			"Species network for embedding the gene tree.", Validate.REQUIRED);
	public final Input<TreeInterface> geneTreeInput = new Input<>("geneTree",
			"Gene tree embedded in the species network.", Validate.REQUIRED);
	public final Input<Embedding> embeddingInput = new Input<>("embedding",
			"Embedding of the gene tree in the species network", (Embedding) null);
	public final Input<Double> ploidyInput = new Input<>("ploidy", "Ploidy (copy number) for this gene (default is 2).",
			2.0);
	public final Input<TaxonSet> taxonSuperSetInput = new Input<>("taxa",
			"Taxon superset associating taxa with gene tree tips", Validate.REQUIRED);

	private boolean needsUpdate;

	// the coalescent times of this gene tree for all species branches
	protected ListMultimap<Integer, Double> coalescentTimes = ArrayListMultimap.create();
	protected ListMultimap<Integer, Double> storedCoalescentTimes = ArrayListMultimap.create();
	// the number of lineages at the tipward end of each species branch
	protected Multiset<Integer> coalescentLineageCounts = HashMultiset.create();
	protected Multiset<Integer> storedCoalescentLineageCounts = HashMultiset.create();

	protected double[][] speciesOccupancy;
	protected double[][] storedSpeciesOccupancy;
	protected double logGammaSum;
	protected double storedLogGammaSum;

	// heirs are the gene tree leaf tip numbers below each gene tree node or species
	// network node
	private Multimap<Node, Integer> geneNodeHeirs = HashMultimap.create();
	private Multimap<NetworkNode, Integer> speciesNodeHeirs = HashMultimap.create();
	private int geneNodeCount;
	private int traversalNodeCount;

	@Override
	public boolean requiresRecalculation() {
		needsUpdate = geneTreeInput.isDirty() || speciesNetworkInput.isDirty();
		return needsUpdate;
	}

	@Override
	public void store() {
		storedCoalescentTimes.clear();
		storedCoalescentLineageCounts.clear();

		storedCoalescentTimes.putAll(coalescentTimes);
		storedCoalescentLineageCounts.addAll(coalescentLineageCounts);

		storedSpeciesOccupancy = new double[speciesOccupancy.length][speciesOccupancy[0].length];
		System.arraycopy(speciesOccupancy, 0, storedSpeciesOccupancy, 0, speciesOccupancy.length);

		storedLogGammaSum = logGammaSum;

		super.store();
	}

	@Override
	public void restore() {
		ListMultimap<Integer, Double> tmpCoalescentTimes = coalescentTimes;
		Multiset<Integer> tmpCoalescentLineageCounts = coalescentLineageCounts;
		double[][] tmpSpeciesOccupancy = speciesOccupancy;
		double tmpLogGammaSum = logGammaSum;

		coalescentTimes = storedCoalescentTimes;
		coalescentLineageCounts = storedCoalescentLineageCounts;
		speciesOccupancy = storedSpeciesOccupancy;
		logGammaSum = storedLogGammaSum;

		storedCoalescentTimes = tmpCoalescentTimes;
		storedCoalescentLineageCounts = tmpCoalescentLineageCounts;
		storedSpeciesOccupancy = tmpSpeciesOccupancy;
		storedLogGammaSum = tmpLogGammaSum;

		super.restore();
	}

	public void initAndValidate() {
		needsUpdate = true;
	}

	protected void computeCoalescentTimes() {
		if (needsUpdate) {
			update();
		}
	}

	private void update() {
		Network speciesNetwork = speciesNetworkInput.get();
		logGammaSum = 0.0;

		final int geneTreeNodeCount = getTree().getNodeCount();
		final int speciesBranchCount = speciesNetwork.getBranchCount();
		speciesOccupancy = new double[geneTreeNodeCount][speciesBranchCount];

		// reset coalescent arrays as these values need to be recomputed after any
		// changes to the species or gene tree
		coalescentLineageCounts.clear();
		coalescentTimes.clear();

		final Node geneTreeRoot = getTree().getRoot();
		final NetworkNode speciesNetworkRoot = speciesNetwork.getRoot();
		final Integer speciesRootBranchNumber = speciesNetworkRoot.gammaBranchNumber;
		recurseCoalescentEvents(geneTreeRoot, speciesNetworkRoot, speciesRootBranchNumber, Double.POSITIVE_INFINITY);

		needsUpdate = false;
	}

	// forward in time recursion, unlike StarBEAST 2
	private void recurseCoalescentEvents(final Node geneTreeNode, final NetworkNode speciesNetworkNode,
			final Integer speciesBranchNumber, final double lastHeight) {
		final double geneNodeHeight = geneTreeNode.getHeight();
		final double speciesNodeHeight = speciesNetworkNode.getHeight();
		final int geneTreeNodeNumber = geneTreeNode.getNr();

		// check if coalescent node occurs in a descendant species network branch
		if (geneNodeHeight < speciesNodeHeight) {
			speciesOccupancy[geneTreeNodeNumber][speciesBranchNumber] += lastHeight - speciesNodeHeight;
			coalescentLineageCounts.add(speciesBranchNumber);
			if (speciesNetworkNode.isReticulation()) {
				final double gamma = speciesNetworkNode.inheritProb;
				if (speciesNetworkNode.gammaBranchNumber.equals(speciesBranchNumber)) {
					logGammaSum += Math.log(gamma);
				} else {
					logGammaSum += Math.log(1.0 - gamma);
				}
			}
			// traversal direction forward in time
			final int traversalNodeNumber = speciesNetworkNode.getTraversalNumber();
			final Integer nextSpeciesBranchNumber = embeddingInput.get().getDirection(geneTreeNodeNumber,
					traversalNodeNumber);
			assert (nextSpeciesBranchNumber >= 0);
			final NetworkNode nextSpeciesNode = speciesNetworkNode.getChildByBranch(nextSpeciesBranchNumber);
			assert nextSpeciesNode != null;
			recurseCoalescentEvents(geneTreeNode, nextSpeciesNode, nextSpeciesBranchNumber, speciesNodeHeight);
		} else if (geneTreeNode.isLeaf()) { // assumes tip node heights are always zero
			speciesOccupancy[geneTreeNodeNumber][speciesBranchNumber] += lastHeight;
			coalescentLineageCounts.add(speciesBranchNumber);
		} else {
			speciesOccupancy[geneTreeNodeNumber][speciesBranchNumber] += lastHeight - geneNodeHeight;
			coalescentTimes.put(speciesBranchNumber, geneNodeHeight);
			for (Node geneChildNode : geneTreeNode.getChildren()) {
				recurseCoalescentEvents(geneChildNode, speciesNetworkNode, speciesBranchNumber, geneNodeHeight);
			}
		}
	}

	public double[][] getSpeciesOccupancy() {
		if (needsUpdate)
			update();
		return speciesOccupancy;
	}

	/**
	 * @return the first tip node which is descendant of
	 * @param gTreeNode this can be in Tree.java as
	 *                  gTreeNode.getGeneTreeTipDescendant()
	 */
	public Node getGeneNodeDescendantTip(Node gTreeNode) {
		final List<Node> gTreeTips = getTree().getExternalNodes(); // tips
		for (Node tip : gTreeTips) {
			Node node = tip;
			while (node != null && !node.equals(gTreeNode)) {
				node = node.getParent();
			}
			if (node != null)
				return tip; // find you!
		}
		return null; // looped all the tips but nothing found
	}

	public void rebuildEmbedding() {
		final Network speciesNetwork = speciesNetworkInput.get();
		traversalNodeCount = speciesNetwork.getTraversalNodeCount();
		geneNodeCount = getTree().getNodeCount();

		getNodeHeirs();

		final Embedding newEmbedding = recurseRebuild(getTree().getRoot(), speciesNetwork.getRoot());
		if (newEmbedding == null) {
			throw new RuntimeException("No valid embedding found");
		}
		setEmbedding(newEmbedding);
	}

	private void getNodeHeirs() {
		// map of species network tip names to species network tip nodes
		final Map<String, NetworkNode> speciesNodeMap = new HashMap<>();
		for (NetworkNode speciesNode : speciesNetworkInput.get().getLeafNodes()) {
			final String speciesName = speciesNode.getLabel();
			speciesNodeMap.put(speciesName, speciesNode);
		}

		// map of gene tree tip names to species network tip nodes
		final Map<String, NetworkNode> geneTipMap = new HashMap<>();
		final TaxonSet taxonSuperSet = taxonSuperSetInput.get();
		for (Taxon species : taxonSuperSet.taxonsetInput.get()) {
			final String speciesName = species.getID();
			final NetworkNode speciesNode = speciesNodeMap.get(speciesName);
			final TaxonSet speciesTaxonSet = (TaxonSet) species;
			for (Taxon geneTip : speciesTaxonSet.taxonsetInput.get()) {
				final String gTipName = geneTip.getID();
				geneTipMap.put(gTipName, speciesNode);
			}
		}

		geneNodeHeirs.clear();
		speciesNodeHeirs.clear();
		for (final Node geneLeaf : getTree().getExternalNodes()) {
			final int gLeafNr = geneLeaf.getNr();
			final String gLeafName = geneLeaf.getID();
			final NetworkNode speciesLeaf = geneTipMap.get(gLeafName);
			// the heir for each gene leaf node is itself
			geneNodeHeirs.put(geneLeaf, gLeafNr);
			// the heirs for each species leaf node is the associated gene leaf nodes
			speciesNodeHeirs.put(speciesLeaf, gLeafNr);
		}

		recurseGeneHeirs(getTree().getRoot());
		for (final NetworkNode speciesLeaf : speciesNetworkInput.get().getLeafNodes()) {
			recurseSpeciesHeirs(speciesLeaf);
		}
	}

	private void recurseGeneHeirs(final Node gTreeNode) {
		for (Node child : gTreeNode.getChildren()) {
			recurseGeneHeirs(child);
			geneNodeHeirs.putAll(gTreeNode, geneNodeHeirs.get(child));
		}
	}

	private void recurseSpeciesHeirs(final NetworkNode sNetNode) {
		for (NetworkNode child : sNetNode.getChildren()) {
			speciesNodeHeirs.putAll(sNetNode, speciesNodeHeirs.get(child));
		}
		for (NetworkNode parent : sNetNode.getParents()) {
			recurseSpeciesHeirs(parent);
		}
	}

	// recursive, return value is a possible gene tree embedding, return null if no
	// valid embedding
	private Embedding recurseRebuild(final Node geneTreeNode, final NetworkNode speciesNetworkNode) {
		if (geneTreeNode.getHeight() < speciesNetworkNode.getHeight()) {
			// embed this gene tree node (lineage) in a descendant species network branch
			final int geneTreeNodeNr = geneTreeNode.getNr();
			final int traversalNodeNr = speciesNetworkNode.getTraversalNumber();
			final Collection<Integer> requiredHeirs = geneNodeHeirs.get(geneTreeNode);

			final Embedding[] altEmbeddings = new Embedding[2];
			double probSum = 0.0;
			int i = 0;
			for (Integer childBranchNr : speciesNetworkNode.childBranchNumbers) {
				final NetworkNode childSpeciesNode = speciesNetworkNode.getChildByBranch(childBranchNr);
				if (speciesNodeHeirs.get(childSpeciesNode).containsAll(requiredHeirs)) {
					altEmbeddings[i] = recurseRebuild(geneTreeNode, childSpeciesNode);
					if (altEmbeddings[i] == null)
						return null;
					altEmbeddings[i].setDirection(geneTreeNodeNr, traversalNodeNr, childBranchNr);

					if (childSpeciesNode.isReticulation()) {
						double childGamma;
						if (childSpeciesNode.gammaBranchNumber.equals(childBranchNr))
							childGamma = childSpeciesNode.getGammaProb();
						else
							childGamma = 1.0 - childSpeciesNode.getGammaProb();
						altEmbeddings[i].probability *= childGamma;
						altEmbeddings[i].probabilitySum *= childGamma;
					}

					probSum += altEmbeddings[i].probabilitySum;
					i++;
				}
			}
			if (i == 0 || probSum == 0.0)
				return null; // for a valid embedding, should never go here

			final double u = Randomizer.nextDouble() * probSum;
			if (u < altEmbeddings[0].probabilitySum) {
				altEmbeddings[0].probabilitySum = probSum;
				return altEmbeddings[0];
			} else {
				altEmbeddings[1].probabilitySum = probSum;
				return altEmbeddings[1];
			}
		} else if (geneTreeNode.isLeaf()) {
			return new Embedding(geneNodeCount, traversalNodeCount);
		} else {
			// embed both children of gene tree node in this species network branch
			Embedding embedding = null;
			for (Node childTreeNode : geneTreeNode.getChildren()) {
				Embedding childEmbedding = recurseRebuild(childTreeNode, speciesNetworkNode);
				if (childEmbedding == null) {
					return null;
				} else if (embedding == null) {
					embedding = childEmbedding;
				} else {
					embedding.mergeWith(childEmbedding);
				}
			}

			return embedding;
		}
	}

	@Override
	public TreeInterface getTree() {
		return geneTreeInput.get();
	}

	@Override
	public Embedding getEmbedding() {
		return embeddingInput.get();
	}

	@Override
	public void setEmbedding(Embedding newEmbedding) {
		embeddingInput.set(newEmbedding);
	}

	@Override
	public double getPloidy() {
		return ploidyInput.get();
	}

	@Override
	public double logGammaSum() {
		computeCoalescentTimes();
		return logGammaSum;
	}

	@Override
	public ListMultimap<Integer, Double> coalescentTimes() {
		return coalescentTimes;
	}

	@Override
	public Multiset<Integer> coalescentLineageCounts() {
		return coalescentLineageCounts;
	}
}
