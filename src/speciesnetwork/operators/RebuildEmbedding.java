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
    private int nLoci;

    @Override
    public void initAndValidate() {
        nLoci = geneTreesInput.get().size();
        geneNodeHeirs = HashMultimap.create();
        speciesNodeHeirs = HashMultimap.create();
    }

    @Override
    public double proposal() {
    	if (nLoci == 0) {
            throw new RuntimeException(String.format("ERROR: no gene trees specified for the operator %s!", getID()));
    	}

        // count the number of alternative traversing choices for the current state
        final int oldChoices = initializeEmbedding(false);
        if (oldChoices < 0)
            throw new RuntimeException("Developer ERROR: current embedding is invalid!");
        final int oldEmbedCount = countEmbedding();

        // make the operation if possible
        double logHR = 0.0;
        if (operatorInput.get() != null) {
            logHR = operatorInput.get().proposal();
            if (logHR == Double.NEGATIVE_INFINITY)
                return Double.NEGATIVE_INFINITY;
        }
        // Update calculation nodes as subsequent operators may depend on state nodes made dirty by this operation.
        // if (!operator.listStateNodes().isEmpty())  // copied from JointOperator
        //   operator.listStateNodes().get(0).getState().checkCalculationNodesDirtiness();

        // then rebuild the embedding AND
        // count the number of alternative traversing choices for the new state
        final int newEmbedCount = countEmbedding();
        final int newChoices = initializeEmbedding(true);
        if (newChoices < 0)
            return Double.NEGATIVE_INFINITY;

        return logHR + (newChoices - oldChoices) * Math.log(2) - Math.log(newEmbedCount) + Math.log(oldEmbedCount);
    }

    @Override
    public List<StateNode> listStateNodes() {
        List<StateNode> stateNodes = new ArrayList<>();

        if (operatorInput.get() != null)
            stateNodes.addAll(operatorInput.get().listStateNodes());
        stateNodes.addAll(super.listStateNodes());

        return stateNodes;
    }

    public int initializeEmbedding(boolean rebuild) {
        final Network speciesNetwork = speciesNetworkInput.get();
        final List<EmbeddedTree> geneTrees = geneTreesInput.get();

        // Tell BEAST that *all* gene trees will be edited
        // doing this for all trees avoids Trie combinatorial explosions
        if (rebuild) {
	        for (EmbeddedTree geneTree : geneTrees) {
	            geneTree.startEditing(this);
	        }
        }

        int nChoices = 0;
        for (int i = 0; i < nLoci; i++) {
            EmbeddedTree geneTree = geneTrees.get(i);

            getNodeHeirs(speciesNetwork, geneTree);

            final int n;
            if (rebuild) {
            	final int traversalNodeCount = speciesNetwork.getTraversalNodeCount();
                geneTree.resetEmbedding(traversalNodeCount, -1);
                n = recurseRebuild(geneTree.getRoot(), speciesNetwork.getRoot(), geneTree, true);
            } else {
                n = recurseRebuild(geneTree.getRoot(), speciesNetwork.getRoot(), geneTree, false);
            }
            // geneTree.printEmbedding();
            if (n < 0)
                return -(i + 1);
            nChoices += n;
        }

        return nChoices;
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

    private int recurseRebuild(final Node geneTreeNode, final NetworkNode speciesNetworkNode,
                               EmbeddedTree geneTree, boolean rebuild) {
        int nChoices = 0;

        // this coalescence node must be embedded in a descendant species network branch
        if (geneTreeNode.getHeight() < speciesNetworkNode.getHeight()) {
            final int traversalNodeNumber = speciesNetworkNode.getTraversalNumber();
            final int geneTreeNodeNumber = geneTreeNode.getNr();
            final Collection<Integer> requiredHeirs = geneNodeHeirs.get(geneTreeNode);
            final List<Integer> compatibleSpeciesBranches = new ArrayList<>();

            for (Integer branchNumber: speciesNetworkNode.childBranchNumbers) {
                final NetworkNode childSpeciesNode = speciesNetworkNode.getChildByBranch(branchNumber);
                if (speciesNodeHeirs.get(childSpeciesNode).containsAll(requiredHeirs)) {
                    compatibleSpeciesBranches.add(branchNumber);
                }
            }
            if (compatibleSpeciesBranches.size() == 0) {
                return -1; // for a valid embedding, should never go here
            } else if (compatibleSpeciesBranches.size() > 1) {
                nChoices++;
            }

            Integer nextSpeciesBranchNumber;
            if (rebuild) {
                final int nextBranchIndex = Randomizer.nextInt(compatibleSpeciesBranches.size());
                nextSpeciesBranchNumber = compatibleSpeciesBranches.get(nextBranchIndex);
                geneTree.setEmbedding(geneTreeNodeNumber, traversalNodeNumber, nextSpeciesBranchNumber);
            } else {
                nextSpeciesBranchNumber = geneTree.getEmbedding(geneTreeNodeNumber, traversalNodeNumber);
            }
            assert (nextSpeciesBranchNumber >= 0);
            final NetworkNode nextSpecies = speciesNetworkNode.getChildByBranch(nextSpeciesBranchNumber);
            final int moreChoices = recurseRebuild(geneTreeNode, nextSpecies, geneTree, rebuild);
            if (moreChoices < 0) return -1;

            return nChoices + moreChoices;
        } else if (geneTreeNode.isLeaf()) {
            return 0;
        } else {
            // embed both gene tree children
            final int leftChoices = recurseRebuild(geneTreeNode.getLeft(), speciesNetworkNode, geneTree, rebuild);
            if (leftChoices < 0) return -1;
            final int rightChoices = recurseRebuild(geneTreeNode.getRight(), speciesNetworkNode, geneTree, rebuild);
            if (rightChoices < 0) return -1;

            return nChoices + leftChoices + rightChoices;
        }
    }

    public int countEmbedding (){
        final Network speciesNetwork = speciesNetworkInput.get();
        final List<EmbeddedTree> geneTrees = geneTreesInput.get();

        int nEmbeddings = 0;
        for (int i = 0; i < nLoci; i++) {
            EmbeddedTree geneTree = geneTrees.get(i);
            getNodeHeirs(speciesNetwork, geneTree);

            final int n = countEmbedding(geneTree.getRoot(), speciesNetwork.getRoot());
            if (n < 0)
                return -(i + 1);  // return which gene tree goes wrong

            nEmbeddings += n;
        }
        return nEmbeddings;
    }

    // recursive
    private int countEmbedding(final Node geneTreeNode, final NetworkNode speciesNetworkNode) {
        if (geneTreeNode.getHeight() < speciesNetworkNode.getHeight()) {
            int nEmbeddings = 0;
            for (Integer childBranchNr: speciesNetworkNode.childBranchNumbers) {
                final NetworkNode childSpeciesNode = speciesNetworkNode.getChildByBranch(childBranchNr);
                final Collection<Integer> requiredHeirs = geneNodeHeirs.get(geneTreeNode);
                if (speciesNodeHeirs.get(childSpeciesNode).containsAll(requiredHeirs)) {
                    final NetworkNode nextSpecies = speciesNetworkNode.getChildByBranch(childBranchNr);
                    final int nNext = countEmbedding(geneTreeNode, nextSpecies);
                    if (nNext < 0)
                        return -1;
                    nEmbeddings += nNext;
                }
            }
            if (nEmbeddings == 0)
                return -1;  // for a valid embedding, should never go here
            return nEmbeddings;
        } else if (geneTreeNode.isLeaf()) {
            return 1;
        } else {
            int nEmbeddings = 1;
            for (Node childTreeNode : geneTreeNode.getChildren()) {
                nEmbeddings *= countEmbedding(childTreeNode, speciesNetworkNode);
                if (nEmbeddings < 0)
                    return -1;
            }
            return nEmbeddings;
        }
    }
}
