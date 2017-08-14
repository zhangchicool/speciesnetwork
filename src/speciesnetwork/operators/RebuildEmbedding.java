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
        // make the operation if possible
        double logHR = 0.0;
        if (operatorInput.get() != null) {
            logHR = operatorInput.get().proposal();
            if (logHR == Double.NEGATIVE_INFINITY)
                return Double.NEGATIVE_INFINITY;
        }

        // then rebuild the embedding
        if (rebuildEmbedding() < 0)
            return Double.NEGATIVE_INFINITY;

        /* TODO calculate HR for embedding */
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

        // Tell BEAST that *all* gene trees will be edited
        // doing this for all trees avoids Trie combinatorial explosions
        for (EmbeddedTree geneTree : geneTrees) {
            geneTree.startEditing(this);
        }

        for (EmbeddedTree geneTree : geneTrees) {
            getNodeHeirs(speciesNetwork, geneTree);

            final int geneNodeCount = geneTree.getNodeCount();
            final int traversalNodeCount = speciesNetwork.getTraversalNodeCount();
            final Embedding initialEmbedding = new Embedding(geneNodeCount);
            initialEmbedding.reset(traversalNodeCount, -1);

            final Set<Embedding> validEmbeddings = recurseRebuild(initialEmbedding, geneTree.getRoot(), speciesNetwork.getRoot());
            final int embeddingCount = validEmbeddings.size();
            if (embeddingCount == 0) return -1;
            final Embedding randomEmbedding = (Embedding) validEmbeddings.toArray()[Randomizer.nextInt(embeddingCount)];
            geneTree.embedding = randomEmbedding;
        }

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
    private Set<Embedding> recurseRebuild(
    		final Embedding embedding,
    		final Node geneTreeNode,
    		final NetworkNode speciesNetworkNode) {
    	final Set<Embedding> embeddingSet = new HashSet<>();
        if (geneTreeNode.getHeight() < speciesNetworkNode.getHeight()) {
            // this coalescent node must be embedded in a descendant species network branch
            final Collection<Integer> requiredHeirs = geneNodeHeirs.get(geneTreeNode);
            final List<Integer> compatibleSpeciesBranches = new ArrayList<>();
            for (Integer childBranchNr: speciesNetworkNode.childBranchNumbers) {
                final NetworkNode childSpeciesNode = speciesNetworkNode.getChildByBranch(childBranchNr);
                if (speciesNodeHeirs.get(childSpeciesNode).containsAll(requiredHeirs)) {
                    compatibleSpeciesBranches.add(childBranchNr);
                }
            }

            for (Integer nextSpeciesBranchNr: compatibleSpeciesBranches) {
	            final int traversalNodeNr = speciesNetworkNode.getTraversalNumber();
	            final int geneTreeNodeNr = geneTreeNode.getNr();
	            embedding.setDirection(geneTreeNodeNr, traversalNodeNr, nextSpeciesBranchNr);
	
	            final NetworkNode nextSpecies = speciesNetworkNode.getChildByBranch(nextSpeciesBranchNr);
	            embeddingSet.addAll(recurseRebuild(embedding, geneTreeNode, nextSpecies));
            }
        } else if (geneTreeNode.isLeaf()) {
            embeddingSet.add(embedding);
        } else {
            // embed both gene tree children
            for (Node childTreeNode : geneTreeNode.getChildren()) {
                embeddingSet.addAll(recurseRebuild(embedding, childTreeNode, speciesNetworkNode));
            }
        }

        return embeddingSet;
    }
}
