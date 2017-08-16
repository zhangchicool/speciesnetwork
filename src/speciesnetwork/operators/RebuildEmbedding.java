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
        if (!rebuildEmbedding())
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

    public boolean rebuildEmbedding() {
        final Network speciesNetwork = speciesNetworkInput.get();
        final List<EmbeddedTree> geneTrees = geneTreesInput.get();

        for (EmbeddedTree geneTree: geneTrees) {
            getNodeHeirs(speciesNetwork, geneTree);

            final int traversalNodeCount = speciesNetwork.getTraversalNodeCount();
            final WeightedArrayList<Integer> wl = recurseRebuild(geneTree.getRoot(), speciesNetwork.getRoot());
            if (wl == null) return false;

            final Integer[] wa = (Integer[]) wl.toArray();
            geneTree.embedding.reset(traversalNodeCount, -1);
            assert wl.size() % 3 == 0;
            final int wlDepth = wl.size() / 3;
            for (int i = 0; i < wlDepth; i++) {
                geneTree.embedding.setDirection(wa[i * 3], wa[i * 3 + 1], wa[i * 3 + 2]);
            }
        }

        return true;
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
    private WeightedArrayList<Integer> recurseRebuild(final Node geneTreeNode, final NetworkNode speciesNetworkNode) {   
        if (geneTreeNode.getHeight() < speciesNetworkNode.getHeight()) {
            // this coalescent node must be embedded in a descendant species network branch
            final int geneTreeNodeNr = geneTreeNode.getNr();
            final int traversalNodeNr = speciesNetworkNode.getTraversalNumber();
            final Collection<Integer> requiredHeirs = geneNodeHeirs.get(geneTreeNode);
            final WeightedArrayList<Integer>[] chains = new WeightedArrayList[2];
            double probabilitySum = 0.0;
            int i = 0;
            for (Integer childBranchNr: speciesNetworkNode.childBranchNumbers) {
                final NetworkNode childSpeciesNode = speciesNetworkNode.getChildByBranch(childBranchNr);
                if (speciesNodeHeirs.get(childSpeciesNode).containsAll(requiredHeirs)) {
                	chains[i] = recurseRebuild(geneTreeNode, childSpeciesNode);
                	if (chains[i] == null) return null;
                	chains[i].add(geneTreeNodeNr);
                	chains[i].add(traversalNodeNr);
                	chains[i].add(childBranchNr);
                	probabilitySum += chains[i].probabilitySum;
                	i++;
                }
            }
            if (i == 0) return null;

            final double randomProb = Randomizer.nextDouble() * probabilitySum;
            if (randomProb < chains[0].probabilitySum) {
            	chains[0].probabilitySum = probabilitySum;
            	return chains[0];
            } else {
            	chains[1].probabilitySum = probabilitySum;
            	return chains[1];
            }
        } else if (geneTreeNode.isLeaf()) {
        	WeightedArrayList<Integer> wl = new WeightedArrayList<>();
        	wl.probability = 1.0;
        	wl.probabilitySum = 1.0;
        	return wl;
        } else {
        	// embed both child nodes of the gene tree node
        	WeightedArrayList<Integer> wl = new WeightedArrayList<>();
        	wl.probability = 1.0;
        	wl.probabilitySum = 1.0;
            for (Node childTreeNode : geneTreeNode.getChildren()) {
                 final WeightedArrayList<Integer> childList = recurseRebuild(childTreeNode, speciesNetworkNode);
                 if (childList == null) return null;
                 wl.addAll(childList);
                 wl.probability *= childList.probability;
                 wl.probabilitySum *= childList.probabilitySum;
            }

            return wl;
        }
    }
}
