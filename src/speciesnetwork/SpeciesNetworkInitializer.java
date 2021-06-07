package speciesnetwork;

import java.util.*;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import speciesnetwork.operators.RebuildEmbedding;
import speciesnetwork.simulator.CoalescentSimulator;

/**
 * @author Huw Ogilvie
 * @author Chi Zhang
 */

@Description("Set a starting point for a species-network analysis.")
public class SpeciesNetworkInitializer extends Tree implements StateNodeInitialiser {

    private enum Method {
        RANDOM("random"), // random starting trees (caterpillar)
        POINT("point"),   // point estimate from the sequences
        USER("user");     // user defined starting state

        Method(final String name) {
            this.ename = name;
        }

        @Override
        public String toString() {
            return ename;
        }

        private final String ename;
    }
    public final Input<Method> initMethod = new Input<>("method",
            "Initializing method (default random)", Method.RANDOM, Method.values());
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "Species network to initialize.", Validate.REQUIRED);
    public final Input<RealParameter> originInput =
            new Input<>("origin", "The time when the process started.", Validate.REQUIRED);
    public final Input<List<EmbeddedTree>> geneTreesInput =
            new Input<>("geneTree", "Gene tree to initialize.", new ArrayList<>());
    public final Input<TaxonSet> taxonSuperSetInput = new Input<>("taxonset",
            "Super-set of taxon sets mapping lineages to species.", Validate.REQUIRED);
    public final Input<RebuildEmbedding> rebuildEmbeddingInput = new Input<>("rebuildEmbedding",
            "Operator which rebuilds embedding of gene trees within species network.", Validate.REQUIRED);
    public final Input<CoalescentSimulator> coalSimulatorInput =
            new Input<>("coalescentSimulator", "Simulate gene trees to initialize.");

    // @Override
    // public void initAndValidate() { super.initAndValidate(); }

    @Override
    public void initStateNodes() {
        Log.info.println("Initializing species network and gene trees.");

        final Method method = initMethod.get();
        switch( method ) {
            case POINT:
                pointInit();
                break;
            case RANDOM:
                randomInit();
                break;
            case USER:
                userStates();
                break;
        }

        final double tOrigin = originInput.get().getValue();
        final double tMRCA = speciesNetworkInput.get().getRoot().getHeight();
        if (tOrigin < tMRCA)
            throw new IllegalArgumentException("Time of origin (" + tOrigin + ") < time of MRCA (" + tMRCA + ")!");

        // initialize embedding for all gene trees
        if (!rebuildEmbeddingInput.get().rebuildEmbedding())
            throw new RuntimeException("Failed to initialize gene tree embedding!");
    }

    private void randomInit() {
        // initialize caterpillar species tree, if no user defined network
        final Network sNetwork = speciesNetworkInput.get();
        // scale the species network according to the time of origin
        final double tOrigin = originInput.get().getValue();
        sNetwork.scale(tOrigin/sNetwork.getOrigin().getHeight());

        if (coalSimulatorInput.get() == null) {
            // map of species network tip names to species network tip nodes
            final Map<String, NetworkNode> speciesTipMap = new HashMap<>();
            for (NetworkNode speciesNode: sNetwork.getLeafNodes()) {
                final String speciesName = speciesNode.getLabel();
                speciesTipMap.put(speciesName, speciesNode);
            }

            // map of gene tip names to species network tip nodes
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

            final double rootHeight = sNetwork.getRoot().getHeight();
            final List<EmbeddedTree> geneTrees = geneTreesInput.get();
            for (final EmbeddedTree gtree : geneTrees) {
                // initialize caterpillar gene tree
                gtree.makeCaterpillar(rootHeight, rootHeight / gtree.getInternalNodeCount(), true);

                // adjust the heights of gene tree tips to be equal to the height of corresponding species tip
                for (Node geneLeaf: gtree.getExternalNodes()) {
                    final String gLeafName = geneLeaf.getID();
                    final NetworkNode speciesLeaf = geneTipMap.get(gLeafName);
                    geneLeaf.setHeight(speciesLeaf.getHeight());
                }
            }
        }
        else {
            // initialize simulated gene trees
            coalSimulatorInput.get().simulate();
        }
    }

    /**
     * build gene trees and species network from alignments
     */
    private void pointInit() {
        throw new RuntimeException("Initialization from point estimate is not yet supported!");
    }

    /**
     * initialize gene trees and species network from user defined starting state
     */
    private void userStates() {
        // nothing to do here at this moment
    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(speciesNetworkInput.get());
        stateNodes.addAll(geneTreesInput.get());
    }
}
