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
    public final Input<Operator> operatorInput = new Input<>("operator",
            "Tree/Network operator to combine into RebuildEmbedding.", Validate.REQUIRED);
    public final Input<List<RebuildEmbedding>> rebuildEmbeddingInput = new Input<>("rebuildEmbedding",
            "Operator which rebuilds embedding of gene tree within species network.", new ArrayList<>());

    @Override
    public void initAndValidate() {
        if (rebuildEmbeddingInput.get().size() == 0)
            throw new RuntimeException("No RebuildEmbedding operator!");
    }

    @Override
    public double proposal() {
        Operator operator = operatorInput.get();
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
        double logHR = operator.proposal();
        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;
        // Update calculation nodes as subsequent operators may depend on state nodes made dirty by this operation.
        // if (!operator.listStateNodes().isEmpty())  // copied from JointOperator
        //   operator.listStateNodes().get(0).getState().checkCalculationNodesDirtiness();

        // then rebuild the embedding
        // and count the number of alternative traversing choices for the new state (n1)
        int newChoices = 0;
        for (RebuildEmbedding reembedOp: reembedOps) {
            final int nChoices = reembedOp.initializeEmbedding();
            if (nChoices < 0)
                return Double.NEGATIVE_INFINITY;
            newChoices += nChoices;
        }

        return logHR + (newChoices - oldChoices) * Math.log(2);
    }

    @Override
    public List<StateNode> listStateNodes() {
        List<StateNode> stateNodes = new ArrayList<>();

        stateNodes.addAll(operatorInput.get().listStateNodes());
        List<RebuildEmbedding> reembedOps = rebuildEmbeddingInput.get();
        for (RebuildEmbedding reembedOp: reembedOps)
            stateNodes.addAll(reembedOp.listStateNodes());

        return stateNodes;
    }
}
