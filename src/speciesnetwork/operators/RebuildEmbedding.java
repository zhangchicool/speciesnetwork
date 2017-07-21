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

    @Override
    public void initAndValidate() {
        // nLoci = geneTreesInput.get().size();
        geneNodeHeirs = HashMultimap.create();
        speciesNodeHeirs = HashMultimap.create();
    }

    @Override
    public double proposal() {
        final List<EmbeddedTree> geneTrees = geneTreesInput.get();

        // count the number of alternative traversing choices for the current state
        int oldChoices = 0;
        for (EmbeddedTree geneTree : geneTrees) {
            oldChoices += geneTree.choicesCount;
        }

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

        // then rebuild the embedding
        if (rebuildEmbedding() < 0)
            return Double.NEGATIVE_INFINITY;

        // count the number of alternative traversing choices for the new state
        int newChoices = 0;
        for (EmbeddedTree geneTree : geneTrees) {
            newChoices += geneTree.choicesCount;
        }

        return logHR + (newChoices - oldChoices) * Math.log(2);
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

        for (EmbeddedTree geneTree : geneTrees) {
            // Tell BEAST that *all* gene trees will be edited
            // doing this for all trees avoids Trie combinatorial explosions
            geneTree.startEditing(this);

            getNodeHeirs(speciesNetwork, geneTree);

            final int traversalNodeCount = speciesNetwork.getTraversalNodeCount();
            geneTree.resetEmbedding(traversalNodeCount, -1);
            final int nChoices = recurseRebuild(geneTree.getRoot(), speciesNetwork.getRoot(), geneTree);
            if (nChoices < 0)
                return -1;  // no valid embedding
            geneTree.choicesCount = nChoices;

            // final int nEmbeddings = countEmbedding(geneTree.getRoot(), speciesNetwork.getRoot());
            // geneTree.embeddingCount = nEmbeddings;
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

    // recursive
    private int recurseRebuild(final Node geneTreeNode, final NetworkNode speciesNetworkNode, EmbeddedTree geneTree) {
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
            int nAlternative = compatibleSpeciesBranches.size() - 1;
            if (nAlternative < 0)
                return -1;  // for a valid embedding, should never go here

            final int nextBranchIndex = Randomizer.nextInt(nAlternative + 1);
            final Integer nextSpeciesBranchNr = compatibleSpeciesBranches.get(nextBranchIndex);
            final int traversalNodeNr = speciesNetworkNode.getTraversalNumber();
            final int geneTreeNodeNr = geneTreeNode.getNr();
            geneTree.setEmbedding(geneTreeNodeNr, traversalNodeNr, nextSpeciesBranchNr);
            // else nextSpeciesBranchNumber = geneTree.getEmbedding(geneTreeNodeNumber, traversalNodeNumber);
            //assert (nextSpeciesBranchNumber >= 0);

            final NetworkNode nextSpecies = speciesNetworkNode.getChildByBranch(nextSpeciesBranchNr);
            final int nNext = recurseRebuild(geneTreeNode, nextSpecies, geneTree);
            if (nNext < 0)
                return -1;
            return nAlternative + nNext;
        } else if (geneTreeNode.isLeaf()) {
            return 0;
        } else {
            int nAlternative = 0;
            // embed both gene tree children
            for (Node childTreeNode : geneTreeNode.getChildren()) {
                final int nNext = recurseRebuild(childTreeNode, speciesNetworkNode, geneTree);
                if (nNext < 0)
                    return -1;
                nAlternative += nNext;
            }
            return nAlternative;
        }
    }

    // recursive
    private int countEmbedding(final Node geneTreeNode, final NetworkNode speciesNetworkNode) {
        if (geneTreeNode.getHeight() < speciesNetworkNode.getHeight()) {
            int nEmbeddings = 0;
            for (Integer childBranchNr: speciesNetworkNode.childBranchNumbers) {
                final NetworkNode childSpeciesNode = speciesNetworkNode.getChildByBranch(childBranchNr);
                final Collection<Integer> requiredHeirs = geneNodeHeirs.get(geneTreeNode);
                if (speciesNodeHeirs.get(childSpeciesNode).containsAll(requiredHeirs)) {
                    final int nNext = countEmbedding(geneTreeNode, childSpeciesNode);
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
