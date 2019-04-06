package speciesnetwork.operators;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.StateNode;
import beast.evolution.alignment.TaxonSet;
import speciesnetwork.GeneTreeInterface;
import speciesnetwork.Network;

/**
 * @author Huw Ogilvie
 * @author Chi Zhang
 */

@Description("Rebuild the embedding of a gene tree in the species network.")
public class RebuildEmbedding extends Operator {
    public final Input<List<GeneTreeInterface>> geneTreesInput = new Input<>("geneTree",
            "The gene tree within the species network.", new ArrayList<>());
    // operator input can be null so that the species network and gene trees are unchanged
    public final Input<Operator> operatorInput = new Input<>("operator",
            "Tree/Network operator to combine into RebuildEmbedding.");


    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        final List<GeneTreeInterface> geneTrees = geneTreesInput.get();

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
        for (GeneTreeInterface geneTree: geneTrees) {
            embeddingLogHR += Math.log(geneTree.getEmbedding().probability) -
            		Math.log(geneTree.getEmbedding().probabilitySum);
            
            // then rebuild the embedding
            try {
            	geneTree.rebuildEmbedding(this);
            } catch (RuntimeException e) {
                return Double.NEGATIVE_INFINITY;
            }
        }

        // finalize hastings ratio of rebuild embedding
        for (final GeneTreeInterface geneTree: geneTrees) {
        	embeddingLogHR -= Math.log(geneTree.getEmbedding().probability) -
        			Math.log(geneTree.getEmbedding().probabilitySum);
        }
        
        return operatorLogHR + embeddingLogHR;
    }

    @Override
    public List<StateNode> listStateNodes() {
        List<StateNode> stateNodes = new ArrayList<>();

        if (operatorInput.get() != null)
            stateNodes.addAll(operatorInput.get().listStateNodes());
		for (GeneTreeInterface geneTree: geneTreesInput.get()) {
            stateNodes.add(geneTree.getEmbedding());        	
        }
        stateNodes.addAll(super.listStateNodes());

        return stateNodes;
    }
}
