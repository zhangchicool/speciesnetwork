package speciesnetwork.operators;

import java.util.*;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.StateNode;

/**
 * @author Chi Zhang
 */

@Description("Combine a gene tree operator with RebuildEmbedding.")
public class JointReembedding extends Operator {
    public Input<Operator> treeOperatorInput = new Input<>("operator",
            "Tree operator to combine into RebuildEmbedding.", Validate.REQUIRED);
    public Input<List<RebuildEmbedding>> rebuildEmbeddingInput = new Input<>("rebuildEmbedding",
            "Operator which rebuilds embedding of gene tree within species network.", new ArrayList<>());

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        Operator treeOp = treeOperatorInput.get();
        List<RebuildEmbedding> reembedOps = rebuildEmbeddingInput.get();

        // count the number of alternative traversing choices for the current state (n0)
        int oldChoices = 0;
        for (RebuildEmbedding reembedOp: reembedOps) {
            final int nChoices = reembedOp.getNumberOfChoices();
            if (nChoices < 0)
                throw new RuntimeException("Developer ERROR: current embedding invalid!");
            oldChoices += nChoices;
        }

        // first make the operation
        double logHR = treeOp.proposal();
        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;
        // Update calculation nodes as subsequent operators may depend on state nodes made dirty by this operation.
        if (!treeOp.listStateNodes().isEmpty())  // copied from JointOperator
            treeOp.listStateNodes().get(0).getState().checkCalculationNodesDirtiness();

        // then rebuild the embedding
        int newChoices = 0;
        for (RebuildEmbedding reembedOp: reembedOps) {
            final int nChoices = reembedOp.initializeEmbedding();
            if (nChoices < 0)
                return Double.NEGATIVE_INFINITY;
            newChoices += nChoices;
            if (!reembedOp.listStateNodes().isEmpty()) // copied from JointOperator
                reembedOp.listStateNodes().get(0).getState().checkCalculationNodesDirtiness();
        }

        return logHR + (newChoices - oldChoices) * Math.log(2);
    }

    @Override
    public List<StateNode> listStateNodes() {
        List<StateNode> stateNodeList = new ArrayList<>();
        List<RebuildEmbedding> reembedOps = rebuildEmbeddingInput.get();

        stateNodeList.addAll(treeOperatorInput.get().listStateNodes());
        for (RebuildEmbedding reembedOp: reembedOps)
            stateNodeList.addAll(reembedOp.listStateNodes());

        return stateNodeList;
    }
}
