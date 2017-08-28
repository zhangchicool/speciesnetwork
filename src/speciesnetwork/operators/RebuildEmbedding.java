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
        for (final EmbeddedTree geneTree : geneTrees) {
            geneTree.startEditing(this);
            logHR += Math.log(geneTree.embedding.probability);
        }

        // then rebuild the embedding
        if (!rebuildEmbedding())
            return Double.NEGATIVE_INFINITY;

        for (final EmbeddedTree geneTree : geneTrees)
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

        final NetworkNode networkRoot = speciesNetwork.getRoot();
        for (EmbeddedTree geneTree : geneTrees) {
            getNodeHeirs(speciesNetwork, geneTree);

            final int geneNodeCount = geneTree.getNodeCount();
            final int networkBranchCount = speciesNetwork.getBranchCount();
            final int traversalNodeCount = speciesNetwork.getTraversalNodeCount();
            double[][] traverseProb = new double[geneNodeCount][networkBranchCount];

            // traverse and store the probabilities of each lineage traversing the network branches
            if (!traverseEmbedding(geneTree.getRoot(), networkRoot, networkRoot.gammaBranchNumber, traverseProb))
                return false;

            // propose a new embedding by traversing each lineage proportional to its traverse probability
            geneTree.embedding.reset(traversalNodeCount);
            geneTree.embedding.probability = 1.0;
            if (!recurseRebuild(geneTree.getRoot(), networkRoot, geneTree, traverseProb))
                return false;
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

    private void recurseGeneHeirs(final Node gTreeNode) {
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

    private boolean recurseRebuild(final Node geneTreeNode, final NetworkNode speciesNetworkNode,
                                   EmbeddedTree geneTree, double[][] traverseProb) {
        final int geneLineageNr = geneTreeNode.getNr();

        if (geneTreeNode.getHeight() < speciesNetworkNode.getHeight()) {
            // embed this gene tree node (lineage) in a descendant species network branch
            final Collection<Integer> requiredHeirs = geneNodeHeirs.get(geneTreeNode);
            final List<Integer> candidateSpeciesBranchNrs = new ArrayList<>();

            double probSum = 0.0;
            for (Integer childBranchNr : speciesNetworkNode.childBranchNumbers) {
                final NetworkNode childSpeciesNode = speciesNetworkNode.getChildByBranch(childBranchNr);
                if (speciesNodeHeirs.get(childSpeciesNode).containsAll(requiredHeirs)) {
                    candidateSpeciesBranchNrs.add(childBranchNr);
                    probSum += traverseProb[geneLineageNr][childBranchNr];
                }
            }
            if (candidateSpeciesBranchNrs.size() == 0)
                return false;  // for a valid embedding, should never go here
            assert (probSum > 0.0);

            // propose a traverse direction proportional to its probability
            final double u = Randomizer.nextDouble();
            Integer nextSpeciesBranchNr = candidateSpeciesBranchNrs.get(0);
            if (u > traverseProb[geneLineageNr][nextSpeciesBranchNr] / probSum)
                nextSpeciesBranchNr = candidateSpeciesBranchNrs.get(1);

            geneTree.embedding.probability *= traverseProb[geneLineageNr][nextSpeciesBranchNr] / probSum;

            final int traversalNodeNr = speciesNetworkNode.getTraversalNumber();
            geneTree.embedding.setDirection(geneLineageNr, traversalNodeNr, nextSpeciesBranchNr);

            final NetworkNode nextSpeciesNode = speciesNetworkNode.getChildByBranch(nextSpeciesBranchNr);
            return recurseRebuild(geneTreeNode, nextSpeciesNode, geneTree, traverseProb);
        }
        else {
            // embed both children of gene tree node in this species network branch
            for (Node childTreeNode : geneTreeNode.getChildren()) {
                if (!recurseRebuild(childTreeNode, speciesNetworkNode, geneTree, traverseProb))
                return false;
            }
            return true;
        }
    }

    private boolean traverseEmbedding(final Node geneTreeNode, final NetworkNode speciesNetworkNode,
                                      final Integer parentSpBranchNr, double[][] traverseProb) {
        final int geneLineageNr = geneTreeNode.getNr();

        if (geneTreeNode.getHeight() < speciesNetworkNode.getHeight()) {
            traverseProb[geneLineageNr][parentSpBranchNr] = 0.0;
            for (Integer childSpBranchNr : speciesNetworkNode.childBranchNumbers) {
                final NetworkNode childSpeciesNode = speciesNetworkNode.getChildByBranch(childSpBranchNr);
                final Collection<Integer> requiredHeirs = geneNodeHeirs.get(geneTreeNode);
                if (speciesNodeHeirs.get(childSpeciesNode).containsAll(requiredHeirs)) {
                    if(!traverseEmbedding(geneTreeNode, childSpeciesNode, childSpBranchNr, traverseProb))
                        return false;

                    traverseProb[geneLineageNr][parentSpBranchNr] += traverseProb[geneLineageNr][childSpBranchNr];
                }
            }
            if (traverseProb[geneLineageNr][parentSpBranchNr] == 0.0)
                return false;  // for a valid embedding, should never go here

            if (speciesNetworkNode.isReticulation()) {
                if (speciesNetworkNode.gammaBranchNumber.equals(parentSpBranchNr))
                    traverseProb[geneLineageNr][parentSpBranchNr] *= speciesNetworkNode.getGammaProb();
                else
                    traverseProb[geneLineageNr][parentSpBranchNr] *= 1.0 - speciesNetworkNode.getGammaProb();
            }
        }
        else {
            traverseProb[geneLineageNr][parentSpBranchNr] = 1.0;
            for (Node childTreeNode : geneTreeNode.getChildren()) {
                if(!traverseEmbedding(childTreeNode, speciesNetworkNode, parentSpBranchNr, traverseProb))
                    return false;

                final int childLineageNr = childTreeNode.getNr();
                traverseProb[geneLineageNr][parentSpBranchNr] *= traverseProb[childLineageNr][parentSpBranchNr];
            }
        }

        return true;
    }
}
