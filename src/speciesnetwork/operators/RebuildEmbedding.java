package speciesnetwork.operators;

import java.util.*;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.StateNode;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import speciesnetwork.EmbeddedTree;
import speciesnetwork.Embedding;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;

/**
 * @author Huw Ogilvie
 * @author Chi Zhang
 */

@Description("Rebuild the embedding of a gene tree in the species network.")
public class RebuildEmbedding extends Operator {
    public final Input<Network> speciesNetworkInput = new Input<>("speciesNetwork",
            "The species network.", Validate.REQUIRED);
    public final Input<TaxonSet> taxonSuperSetInput = new Input<>("taxonset",
            "Super-set of taxon sets mapping lineages to species.", Validate.REQUIRED);
    public final Input<List<EmbeddedTree>> geneTreesInput = new Input<>("geneTree",
            "The gene tree within the species network.", new ArrayList<>());
    // operator input can be null so that the species network and gene trees are unchanged
    public final Input<Operator> operatorInput = new Input<>("operator",
            "Tree/Network operator to combine into RebuildEmbedding.");

    // heirs are the gene tree leaf tip numbers below each gene tree node or species network node
    private Multimap<Node, Integer> geneNodeHeirs;
    private Multimap<NetworkNode, Integer> speciesNodeHeirs;

    @Override
    public void initAndValidate() {
        // nLoci = geneTreesInput.get().size();
        geneNodeHeirs = HashMultimap.create();
        speciesNodeHeirs = HashMultimap.create();
    }

    @Override
    public double proposal() {
        final List<EmbeddedTree> geneTrees = geneTreesInput.get();
        final int nGeneTrees = geneTrees.size();

        // make the operation if possible
        double logHR = 0.0;
        if (operatorInput.get() != null) {
            logHR = operatorInput.get().proposal();
            if (logHR == Double.NEGATIVE_INFINITY)
                return Double.NEGATIVE_INFINITY;
        }

        // Tell BEAST that *all* gene trees will be edited
        // doing this for all trees avoids Trie combinatorial explosions
        for (int i = 0; i < nGeneTrees; i++) {
        	final EmbeddedTree geneTree = geneTrees.get(i);
            geneTree.startEditing(this);
            logHR += Math.log(geneTree.embedding.probability);
        }

        // then rebuild the embedding
        if (rebuildEmbedding() < 0)
            return Double.NEGATIVE_INFINITY;

        for (final EmbeddedTree geneTree: geneTrees)
            logHR -= Math.log(geneTree.embedding.probability);

        return logHR;
    }

    @Override
    public List<StateNode> listStateNodes() {
        List<StateNode> stateNodes = new ArrayList<>();

        if (operatorInput.get() != null)
            stateNodes.addAll(operatorInput.get().listStateNodes());
        stateNodes.addAll(super.listStateNodes());

        return stateNodes;
    }

    public int rebuildEmbedding() {
        final Network speciesNetwork = speciesNetworkInput.get();
        final List<EmbeddedTree> geneTrees = geneTreesInput.get();

        int embeddingSum = 0;
        for (EmbeddedTree geneTree : geneTrees) {
            getNodeHeirs(speciesNetwork, geneTree);

            final int geneNodeCount = geneTree.getNodeCount();
            final int traversalNodeCount = speciesNetwork.getTraversalNodeCount();
            final Embedding initialEmbedding = new Embedding(geneNodeCount);
            initialEmbedding.reset(traversalNodeCount, -1);
            final Set<Embedding> embeddingSet = new HashSet<>();
            embeddingSet.add(initialEmbedding);

            if (recurseRebuild(embeddingSet, geneTree.getRoot(), speciesNetwork.getRoot())) {
	            Embedding[] embeddingArray = new Embedding[embeddingSet.size()];
	            int i = 0;
	            double probabilitySum = 0.0;
	            for (Embedding e: embeddingSet) {
	            	probabilitySum += e.probability;
	            	embeddingArray[i] = e;
	            	i++;
	            }
	        	assert embeddingArray.length == 1 || speciesNetwork.getReticulationNodeCount() > 0;
	        	final Embedding randomEmbedding = embeddingArray[Randomizer.nextInt(embeddingArray.length)];
	        	randomEmbedding.probability = randomEmbedding.probability / probabilitySum; // normalize probability
	            geneTree.embedding = randomEmbedding;
	            embeddingSum += embeddingArray.length;
            } else {
            	return -1;
            }
        }

        System.out.println("Total number of embeddings: " + embeddingSum + ", reticulation count: " + speciesNetwork.getReticulationNodeCount());
        return 0;
    }

	private void getNodeHeirs(final Network speciesNetwork, final EmbeddedTree geneTree) {
        // map of species network tip names to species network tip nodes
        final Map<String, NetworkNode> speciesNodeMap = new HashMap<>();
        for (NetworkNode speciesNode: speciesNetwork.getLeafNodes()) {
            final String speciesName = speciesNode.getLabel();
            speciesNodeMap.put(speciesName, speciesNode);
        }

        // map of gene tree tip names to species network tip nodes
        final Map<String, NetworkNode> geneTipMap = new HashMap<>();
        final TaxonSet taxonSuperSet = taxonSuperSetInput.get();
        for (Taxon species: taxonSuperSet.taxonsetInput.get()) {
            final String speciesName = species.getID();
            final NetworkNode speciesNode = speciesNodeMap.get(speciesName);
            final TaxonSet speciesTaxonSet = (TaxonSet) species;
            for (Taxon geneTip: speciesTaxonSet.taxonsetInput.get()) {
                final String gTipName = geneTip.getID();
                geneTipMap.put(gTipName, speciesNode);
            }
        }

        geneNodeHeirs.clear();
        speciesNodeHeirs.clear();
        for (final Node geneLeaf: geneTree.getExternalNodes()) {
            final int gLeafNr = geneLeaf.getNr();
            final String gLeafName = geneLeaf.getID();
            final NetworkNode speciesLeaf = geneTipMap.get(gLeafName);
            // the heir for each gene leaf node is itself
            geneNodeHeirs.put(geneLeaf, gLeafNr);
            // the heirs for each species leaf node is the associated gene leaf nodes
            speciesNodeHeirs.put(speciesLeaf, gLeafNr);
        }

        recurseGeneHeirs(geneTree.getRoot());
        for (final NetworkNode speciesLeaf: speciesNetwork.getLeafNodes()) {
            recurseSpeciesHeirs(speciesLeaf);
        }
    }

    private void recurseGeneHeirs (final Node gTreeNode) {
        for (Node child : gTreeNode.getChildren()) {
            recurseGeneHeirs(child);
            geneNodeHeirs.putAll(gTreeNode, geneNodeHeirs.get(child));
        }
    }

    private void recurseSpeciesHeirs(final NetworkNode sNetNode) {
        for (NetworkNode child: sNetNode.getChildren()) {
            speciesNodeHeirs.putAll(sNetNode, speciesNodeHeirs.get(child));
        }
        for (NetworkNode parent: sNetNode.getParents()) {
            recurseSpeciesHeirs(parent);
        } 
    }

    // recursive, return value is the multiplication of gamma probabilities
    private boolean recurseRebuild(final Set<Embedding> embeddingSet, final Node geneTreeNode, final NetworkNode speciesNetworkNode) {   
        if (geneTreeNode.getHeight() < speciesNetworkNode.getHeight()) {
            // this coalescent node must be embedded in a descendant species network branch
            final int geneTreeNodeNr = geneTreeNode.getNr();
            final int traversalNodeNr = speciesNetworkNode.getTraversalNumber();
            final Collection<Integer> requiredHeirs = geneNodeHeirs.get(geneTreeNode);
        	boolean hasValidEmbedding = false;
            for (Integer childBranchNr: speciesNetworkNode.childBranchNumbers) {
                final NetworkNode childSpeciesNode = speciesNetworkNode.getChildByBranch(childBranchNr);
                if (speciesNodeHeirs.get(childSpeciesNode).containsAll(requiredHeirs)) {
		            if (hasValidEmbedding) {
		            	// already set embedding so duplicate all existing embeddings
		            	final Set<Embedding> alternativeSet = new HashSet<>();
		            	for (Embedding e: embeddingSet)
		            		alternativeSet.add(new Embedding(e));
		            	// add them to the set if the embedding is valid
		            	if (recurseRebuild(alternativeSet, geneTreeNode, childSpeciesNode))
		            		for (Embedding e: alternativeSet) {
		            			e.setDirection(geneTreeNodeNr, traversalNodeNr, childBranchNr);
		            			embeddingSet.add(e);
		            		};
		            } else if (recurseRebuild(embeddingSet, geneTreeNode, childSpeciesNode)) {
		            	hasValidEmbedding = true;
		            	// set the direction if the embedding is valid
		            	for (Embedding e: embeddingSet)
		            		e.setDirection(geneTreeNodeNr, traversalNodeNr, childBranchNr);
		            }
                }
            }

            return hasValidEmbedding; // has at least one possible embedding
        } else {
        	// embed both child nodes of the gene tree node
            for (Node childTreeNode : geneTreeNode.getChildren()) {
                if (!recurseRebuild(embeddingSet, childTreeNode, speciesNetworkNode))
            		return false;
            }

            return true;
        }
    }
}
