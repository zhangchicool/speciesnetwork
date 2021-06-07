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

    // heirs are the gene tree leaf numbers below each gene tree node or species network node
    private Multimap<Node, Integer> geneNodeHeirs;
    private Multimap<NetworkNode, Integer> speciesNodeHeirs;
    private int geneNodeCount;
    private int traversalNodeCount;

    @Override
    public void initAndValidate() {
        // nLoci = geneTreesInput.get().size();
        geneNodeHeirs = HashMultimap.create();
        speciesNodeHeirs = HashMultimap.create();
    }

    @Override
    public double proposal() {
        final List<EmbeddedTree> geneTrees = geneTreesInput.get();

        // make the operation if possible
        double operatorLogHR = 0.0;
        if (operatorInput.get() != null) {
            operatorLogHR = operatorInput.get().proposal();
            if (operatorLogHR == Double.NEGATIVE_INFINITY)
                return Double.NEGATIVE_INFINITY;
        }

        // Tell BEAST that *all* gene trees will be edited
        // doing this for all trees avoids Trie combinatorial explosions
        double embeddingLogHR = 0.0;
        for (final EmbeddedTree geneTree: geneTrees) {
            geneTree.startEditing(this);
            embeddingLogHR += Math.log(geneTree.embedding.probability) - Math.log(geneTree.embedding.probabilitySum);
        }

        // then rebuild the embedding
        if (!rebuildEmbedding())
            return Double.NEGATIVE_INFINITY;

        // finalize hastings ratio of rebuild embedding
        for (final EmbeddedTree geneTree: geneTrees) {
            embeddingLogHR -= Math.log(geneTree.embedding.probability) - Math.log(geneTree.embedding.probabilitySum);
        }
        
        return operatorLogHR + embeddingLogHR;
    }

    @Override
    public List<StateNode> listStateNodes() {
        List<StateNode> stateNodes = new ArrayList<>();

        final Operator op = operatorInput.get();
        if (op != null)
            stateNodes.addAll(op.listStateNodes());
        stateNodes.addAll(super.listStateNodes());

        return stateNodes;
    }

    public boolean rebuildEmbedding() {
        final List<EmbeddedTree> geneTrees = geneTreesInput.get();
        final Network speciesNetwork = speciesNetworkInput.get();
        traversalNodeCount = speciesNetwork.getInternalNodeCount();

        for (EmbeddedTree geneTree: geneTrees) {
            geneNodeCount = geneTree.getNodeCount();
            getNodeHeirs(speciesNetwork, geneTree);

            final Embedding newEmbedding = recurseRebuild(geneTree.getRoot(), speciesNetwork.getRoot());
            if (newEmbedding == null) return false;

            geneTree.embedding = newEmbedding;
        }

        return true;
    }

    private void getNodeHeirs(final Network speciesNetwork, final EmbeddedTree geneTree) {
        // map of species network tip names to species network tip nodes
        final Map<String, NetworkNode> speciesTipMap = new HashMap<>();
        for (NetworkNode speciesNode: speciesNetwork.getLeafNodes()) {
            final String speciesName = speciesNode.getLabel();
            speciesTipMap.put(speciesName, speciesNode);
        }

        // map of gene tree tip names to species network tip nodes
        final Map<String, NetworkNode> geneTipMap = new HashMap<>();
        final TaxonSet taxonSuperSet = taxonSuperSetInput.get();
        for (Taxon species: taxonSuperSet.taxonsetInput.get()) {
            final String speciesName = species.getID();
            final NetworkNode speciesNode = speciesTipMap.get(speciesName);
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
        // then map heirs for all internal nodes, recursively
        recurseGeneHeirs(geneTree.getRoot());
        recurseSpeciesHeirs(speciesNetwork.getRoot());
    }

    private void recurseGeneHeirs(final Node gTreeNode) {
        for (Node child : gTreeNode.getChildren()) {
            recurseGeneHeirs(child);
            geneNodeHeirs.putAll(gTreeNode, geneNodeHeirs.get(child));
        }
    }

    private void recurseSpeciesHeirs(final NetworkNode sNetNode) {
        for (NetworkNode child: sNetNode.getChildren()) {
            recurseSpeciesHeirs(child);
            speciesNodeHeirs.putAll(sNetNode, speciesNodeHeirs.get(child));
        }
    }

    // recursive, return value is a possible gene tree embedding, return null if no valid embedding
    private Embedding recurseRebuild(final Node geneTreeNode, final NetworkNode speciesNetworkNode) {
        if (geneTreeNode.isLeaf() && speciesNetworkNode.isLeaf()) {
            // reached the gene tree tip and species tip (height >= 0)
            return new Embedding(geneNodeCount, traversalNodeCount);
        }
        else if (geneTreeNode.getHeight() <= speciesNetworkNode.getHeight()) {
            // current gene tree node occurs in a descendant branch of current species node
            final int geneTreeNodeNr = geneTreeNode.getNr();
            final int traversalNodeNr = speciesNetworkNode.getTraversalNumber();
            final Collection<Integer> requiredHeirs = geneNodeHeirs.get(geneTreeNode);

            // there are at most two possible embeddings for this gene lineage
            final Embedding[] altEmbeddings = new Embedding[2];
            double probSum = 0.0;
            int i = 0;
            for (Integer childBranchNr: speciesNetworkNode.childBranchNumbers) {
                final NetworkNode childSpeciesNode = speciesNetworkNode.getChildByBranch(childBranchNr);
                if (speciesNodeHeirs.get(childSpeciesNode).containsAll(requiredHeirs)) {
                    // a valid embedding is possible, move on to find out
                    altEmbeddings[i] = recurseRebuild(geneTreeNode, childSpeciesNode);
                    if (altEmbeddings[i] == null) return null;

                    // the lineage passes through current species node and goes to this species child branch (forward in time)
                    altEmbeddings[i].setDirection(geneTreeNodeNr, traversalNodeNr, childBranchNr);

                    // deal with traversal probabilities
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
            if (i == 0 || probSum == 0.0) return null;  // for a valid embedding, should never go here

            // propose the embedding proportional to its probability
            final double u = Randomizer.nextDouble() * probSum;
            if (u < altEmbeddings[0].probabilitySum) {
                altEmbeddings[0].probabilitySum = probSum;
                return altEmbeddings[0];
            } else {
                altEmbeddings[1].probabilitySum = probSum;
                return altEmbeddings[1];
            }
        }
        else {
            // current gene tree node occurs above current species node
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
}
