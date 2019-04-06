package speciesnetwork;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multiset;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
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
	public final Input<TaxonSet> taxonSuperSetInput = new Input<>("taxa",
			"Taxon superset associating taxa with gene tree tips", Validate.REQUIRED);
	public Input<Embedding> embeddingInput = new Input<Embedding>("embedding",
			"Embedding between geneTree and speciesNetwork", (Embedding) null, Embedding.class);

	// the coalescent times of this gene tree for all species branches
	protected ListMultimap<Integer, Double> coalescentTimes = ArrayListMultimap.create();

	protected double[][] speciesOccupancy;
	protected double logGammaSum;

	// heirs are the gene tree leaf tip numbers below each gene tree node or species
	// network node
	private Multimap<NetworkNode, Integer> speciesNodeHeirs = HashMultimap.create();
	private Multimap<NetworkNode, Integer> storedSpeciesNodeHeirs = null;

	private final Multiset<Integer> coalescentLineageCounts = HashMultiset.create();
	private int traversalNodeCount;
	private int geneNodeCount;
  private final Map<String, NetworkNode>  geneTipMap = new HashMap<>();

	@Override
	public boolean requiresRecalculation() {
		boolean needsUpdate;
		needsUpdate = geneTreeInput.isDirty() || speciesNetworkInput.isDirty() || embeddingInput.isDirty();
		return needsUpdate;
	}

  @Override
	public void initAndValidate() {
		Network speciesNetwork = speciesNetworkInput.get();

		geneNodeCount = getTree().getNodeCount();
		traversalNodeCount = speciesNetwork.getTraversalNodeCount();
		Embedding e = embeddingInput.get();
		if (e == null) {
			embeddingInput.set(new Embedding(geneNodeCount, traversalNodeCount));			
		} else if (e.getDimension() == 0) {
			embeddingInput.get().assignFrom(new Embedding(geneNodeCount, traversalNodeCount));
		} else {
			if (e.geneNodeCount != geneNodeCount || e.traversalNodeCount != traversalNodeCount) {
				throw new RuntimeException("Wrong embedding shape");
			}
		}
    
    // map of species network tip names to species network tip nodes
		final Map<String, NetworkNode> speciesNodeMap = new HashMap<>();
		for (NetworkNode speciesNode : speciesNetworkInput.get().getLeafNodes()) {
			final String speciesName = speciesNode.getLabel();
			speciesNodeMap.put(speciesName, speciesNode);
		}

		// map of gene tree tip names to species network tip nodes
		geneTipMap.clear();
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
	}

	void update() {
		Network speciesNetwork = speciesNetworkInput.get();
		logGammaSum = 0.0;

		geneNodeCount = getTree().getNodeCount();
		traversalNodeCount = speciesNetwork.getTraversalNodeCount();
		final int speciesBranchCount = speciesNetwork.getBranchCount();
		speciesOccupancy = new double[geneNodeCount][speciesBranchCount];

		// reset coalescent arrays as these values need to be recomputed after any
		// changes to the species or gene tree
		coalescentLineageCounts.clear();
		coalescentTimes.clear();

		final Node geneTreeRoot = getTree();
		final NetworkNode speciesNetworkRoot = speciesNetwork.getRoot();
		final Integer speciesRootBranchNumber = speciesNetworkRoot.gammaBranchNumber;

    if (speciesNetwork.isDirty()) {
      getNodeHeirs();      
    }

		assert getEmbedding().geneNodeCount == geneNodeCount;
		assert getEmbedding().traversalNodeCount == traversalNodeCount;
		recurseCoalescentEvents(geneTreeRoot, speciesRootBranchNumber, Double.POSITIVE_INFINITY);
	}

	public Node getTree() {
		return geneTreeInput.get().getRoot();
	}

	// forward in time recursion
	private void recurseCoalescentEvents(final Node geneTreeNode, final Integer speciesBranchNumber,
			final double lastHeight) {
		Network network = speciesNetworkInput.get();
		NetworkNode speciesNetworkNode = network.getNode(network.getNodeNumber(speciesBranchNumber));
		final double speciesNodeHeight = speciesNetworkNode.getHeight();

		final int geneTreeNodeNumber = geneTreeNode.getNr();
		final double geneNodeHeight = geneTreeNode.getHeight();

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
			final Integer nextSpeciesBranchNumber = getEmbedding().getDirection(geneTreeNodeNumber,
					traversalNodeNumber);
			recurseCoalescentEvents(geneTreeNode, nextSpeciesBranchNumber, speciesNodeHeight);
		} else if (geneTreeNode.isLeaf()) {
			// TODO: assumes tip node heights are always zero
			speciesOccupancy[geneTreeNodeNumber][speciesBranchNumber] += lastHeight;
			coalescentLineageCounts.add(speciesBranchNumber);
		} else {
			speciesOccupancy[geneTreeNodeNumber][speciesBranchNumber] += lastHeight - geneNodeHeight;
			coalescentTimes.put(speciesBranchNumber, geneNodeHeight);
			for (Node geneChildNode : geneTreeNode.getChildren()) {
				recurseCoalescentEvents(geneChildNode, speciesBranchNumber, geneNodeHeight);
			}
		}
	}

  @Override
	public void rebuildEmbedding(Operator operator) {
		final Network speciesNetwork = speciesNetworkInput.get();
		traversalNodeCount = speciesNetwork.getTraversalNodeCount();
		geneNodeCount = getTree().getNodeCount();

    if (speciesNetwork.isDirty()) {
      getNodeHeirs();      
    }
		final Embedding newEmbedding = recurseRebuild(getTree(), speciesNetwork.getRoot());
		if (newEmbedding == null) {
			throw new RuntimeException("No valid embedding found");
		}
		newEmbedding.stored = getEmbedding(); 
		if (operator != null) {
			embeddingInput.get().startEditing(operator);
		}
		embeddingInput.get().assignFrom(newEmbedding);
		update();
	}

	void getNodeHeirs() {
		speciesNodeHeirs.clear();
		for (final Node geneLeaf : getTree().getAllLeafNodes()) {
			final int gLeafNr = geneLeaf.getNr();
			final String gLeafName = geneLeaf.getID();
			final NetworkNode speciesLeaf = geneTipMap.get(gLeafName);
			// the heirs for each species leaf node are the associated gene leaf nodes
			speciesNodeHeirs.put(speciesLeaf, gLeafNr);
		}

		recurseSpeciesHeirs(speciesNetworkInput.get().getRoot(), speciesNodeHeirs);
	}

	void recurseSpeciesHeirs(final NetworkNode sNetNode, Multimap<NetworkNode, Integer> geneHeirs) {
		if (!geneHeirs.get(sNetNode).isEmpty()) {
			return;
		}
		for (NetworkNode child : sNetNode.getChildren()) {
			recurseSpeciesHeirs(child, geneHeirs);
			geneHeirs.putAll(sNetNode, geneHeirs.get(child));
		}

	}

	// recursive, return value is a possible gene tree embedding, return null if no
	// valid embedding
	private Embedding recurseRebuild(final Node geneTreeNode, final NetworkNode speciesNetworkNode) {
		if (geneTreeNode.getHeight() < speciesNetworkNode.getHeight()) {
			// embed this gene tree node (lineage) in a descendant species network branch
			final int geneTreeNodeNr = geneTreeNode.getNr();
			final int traversalNodeNr = speciesNetworkNode.getTraversalNumber();

			double probSum = 0.0;
			List<Embedding> altEmbeddings = new ArrayList<>();
			for (Integer childBranchNr : speciesNetworkNode.childBranchNumbers) {
				final NetworkNode childSpeciesNode = speciesNetworkNode.getChildByBranch(childBranchNr);
				if (speciesSubNetworkContainsGeneSubTree(childSpeciesNode, geneTreeNode)) {
					Embedding em = recurseRebuild(geneTreeNode, childSpeciesNode);

					if (em == null)
						return null;
					em.setDirection(geneTreeNodeNr, traversalNodeNr, childBranchNr);

					if (childSpeciesNode.isReticulation()) {
						double childGamma;
						if (childSpeciesNode.gammaBranchNumber.equals(childBranchNr))
							childGamma = childSpeciesNode.getGammaProb();
						else
							childGamma = 1.0 - childSpeciesNode.getGammaProb();
						em.probability *= childGamma;
						em.probabilitySum *= childGamma;
					}

					altEmbeddings.add(em);
					probSum += em.probabilitySum;
				}
			}
			if (altEmbeddings.isEmpty() || probSum == 0.0)
				return null; // for a valid embedding, should never go here

			double u = Randomizer.nextDouble() * probSum;
			for (Embedding em : altEmbeddings) {
				if (u < em.probabilitySum) {
					em.probabilitySum = probSum;
					return em;
				} else {
					u -= em.probabilitySum;
				}
			}
			// There was some u left over after all the probability sums? This should never
			// happen.
			throw new RuntimeException("Random selection of sub-embedding failed.");
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

	private boolean speciesSubNetworkContainsGeneSubTree(NetworkNode childSpeciesNode, Node geneTreeNode) {
		Collection<Integer> speciesGeneHeirs = speciesNodeHeirs.get(childSpeciesNode);
		if (geneTreeNode.isLeaf()) {
      return speciesGeneHeirs.contains(geneTreeNode.getNr());
		}
		for (Node l : geneTreeNode.getAllLeafNodes()) {
			if (!speciesGeneHeirs.contains(l.getNr())) {
				return false;
			}
		}
		return true;
	}

	@Override
	public double logGammaSum() {
    if (requiresRecalculation()) {
      update();
    }
		return logGammaSum;
	}

	@Override
	public ListMultimap<Integer, Double> coalescentTimes() {
    if (requiresRecalculation()) {
      update();
    }
		return coalescentTimes;
	}

	@Override
	public Multiset<Integer> coalescentLineageCounts() {
    if (requiresRecalculation()) {
      update();
    }
		return coalescentLineageCounts;
	}

	@Override
	public Embedding getEmbedding() {
		return embeddingInput.get();
	}
	
  @Override
	protected void store() {
	  super.store();
    storedSpeciesNodeHeirs = speciesNodeHeirs;
    speciesNodeHeirs = HashMultimap.create();
	}
	
  @Override
	protected void restore() {
    super.restore();
    speciesNodeHeirs = storedSpeciesNodeHeirs;
	}
}
