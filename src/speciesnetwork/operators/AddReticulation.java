package speciesnetwork.operators;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.util.Randomizer;

import speciesnetwork.Network;
import speciesnetwork.NetworkNode;
import speciesnetwork.SanityChecks;

/**
 * This proposal adds a reticulation event by connecting two existing branches (with length l1 and l2) with a new branch.
 * The same branch can be picked twice (and forms a loop to that branch). The cutting proportion of each picked branch by
 * the connecting point, w1 and w2 ~ Uniform(0,1). Let l11 = l1 * w1, l12 = l1 * (1-w1), l21 = l2 * w2, l22 = l2 * (1-w2)
 * If the root branch is picked (with length l1 unknown), let l11 = w1 ~ exp(b), and l12 = unknown.
 * The direction of the new branch is determined by the two connecting points, the higher is speciation node, and the
 * lower is reticulation node. The gamma prob r = w3 ~ Uniform(0,1).
 * The Jacobian is l1 * l2. (li = 1 if the branch is a root branch. If the root branch is picked twice, the Jacobian = 1)
 *
 * The AddReticulation and DeleteReticulation are chosen with equal prob. If there is no reticulation in the network,
 * the DeleteReticulation move is aborted.
 * Let k be the number of branches in the current network. The probability of adding this branch is (1/k)(1/k)
 * Let m be the number of reticulation branches in the proposed network. The probability of selecting the same branch to
 * remove is (1/m).
 * The Hastings ratio is (1/m) / (1/k)(1/k)(g1)(g2)(g3)) = k^2 / (m * g1 * g2 * g3).
 * g1 = b * exp(-b * l11) if root, otherwise g1 = 1 (uniform density); g2 = b * exp(-b * l21) if root, otherwise g2 = 1.
 * g3 = 1.
 *
 * @author Chi Zhang
 */

@Description("Add a reticulation branch to the species network.")
public class AddReticulation extends Operator {
    public Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);
    public Input<List<RebuildEmbedding>> rebuildEmbeddingInput = new Input<>("rebuildEmbedding",
            "Operator which rebuilds embedding of gene tree within species network.", new ArrayList<>());

    private final double lambda = 1.0;  // rate of exponential distribution

    // empty constructor to facilitate construction by XML + initAndValidate
    public AddReticulation() {
    }

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        final Network speciesNetwork = speciesNetworkInput.get();
        final List<RebuildEmbedding> reembedOps = rebuildEmbeddingInput.get();

        SanityChecks.checkNetworkSanity(speciesNetwork.getRoot());

        // count the number of alternative traversing choices for the current state (n0)
        int oldChoices = 0;
        for (RebuildEmbedding reembedOp: reembedOps) {
            final int nChoices = reembedOp.getNumberOfChoices();
            if (nChoices < 0)
                throw new RuntimeException("Developer ERROR: current embedding invalid!");
            oldChoices += nChoices;
        }

        // number of branches in the current network
        final int nBranches = speciesNetwork.getBranchCount();  // k

        // pick two branches randomly, including the root branch
        final int pickedBranchNr1 = Randomizer.nextInt(nBranches);
        final int pickedBranchNr2 = Randomizer.nextInt(nBranches);  // allow picking the same branch

        // get the nodes associated with each branch
        final int pickedNodeNr1 = speciesNetwork.getNodeNumber(pickedBranchNr1);
        NetworkNode pickedNode1 = speciesNetwork.getNode(pickedNodeNr1);
        final int pickedNodeNr2 = speciesNetwork.getNodeNumber(pickedBranchNr2);
        NetworkNode pickedNode2 = speciesNetwork.getNode(pickedNodeNr2);
        NetworkNode pickedParent1 = pickedNode1.getParentByBranch(pickedBranchNr1);  // null if pickedNode1 is root
        NetworkNode pickedParent2 = pickedNode2.getParentByBranch(pickedBranchNr2);  // null if pickedNode2 is root

        // propose the attaching position at each branch
        final double l1, l2, l11, l21;
        double logProposalRatio = 0.0;

        if (pickedNode1.isRoot()) {
            l1 = 1;
            l11 = Randomizer.nextExponential(lambda);
            logProposalRatio += lambda * l11 - Math.log(lambda);
        } else {
            l1 = pickedParent1.getHeight() - pickedNode1.getHeight();
            l11 = l1 * Randomizer.nextDouble();
        }
        if (pickedNode2.isRoot()) {
            l2 = 1;
            l21 = Randomizer.nextExponential(lambda);
            logProposalRatio += lambda * l21 - Math.log(lambda);
        } else {
            l2 = pickedParent2.getHeight() - pickedNode2.getHeight();
            l21 = l2 * Randomizer.nextDouble();
        }
        logProposalRatio += Math.log(l1) + Math.log(l2);  // the Jacobian

        // final double middleNodeHeight1 = pickedNode1.getHeight() + l11;
        // final double middleNodeHeight2 = pickedNode2.getHeight() + l21;

        // start moving
        speciesNetwork.startEditing(this);

        // create two new nodes
        NetworkNode middleNode1 = new NetworkNode(speciesNetwork);
        NetworkNode middleNode2 = new NetworkNode(speciesNetwork);

        // set height
        middleNode1.setHeight(pickedNode1.getHeight() + l11);
        middleNode2.setHeight(pickedNode2.getHeight() + l21);

        // add a branch joining the two middle nodes (picked branches)
        if (middleNode1.getHeight() > middleNode2.getHeight()) {
            speciesNetwork.addReticulationBranch(middleNode1, middleNode2, pickedBranchNr1, pickedBranchNr2);
            middleNode2.setGamma(Randomizer.nextDouble());
        } else {
            speciesNetwork.addReticulationBranch(middleNode2, middleNode1, pickedBranchNr2, pickedBranchNr1);
            middleNode1.setGamma(Randomizer.nextDouble());
        }

        // number of reticulation branches in the proposed network
        final int nReticulationBranches = 2 * speciesNetwork.getReticulationNodeCount();  // m
        logProposalRatio += 2 * Math.log(nBranches) - Math.log(nReticulationBranches);

        SanityChecks.checkNetworkSanity(speciesNetwork.getRoot());

        // update the embedding in the new species network
        int newChoices = 0;
        for (RebuildEmbedding reembedOp: reembedOps) {
            final int nChoices = reembedOp.initializeEmbedding();
            if (nChoices < 0)
                return Double.NEGATIVE_INFINITY;
            newChoices += nChoices;
            // System.out.println(String.format("Gene tree %d: %d choices", i, nChoices));
            if (!reembedOp.listStateNodes().isEmpty()) // copied from JointOperator
                reembedOp.listStateNodes().get(0).getState().checkCalculationNodesDirtiness();
        }
        logProposalRatio += (newChoices - oldChoices) * Math.log(2);

        return logProposalRatio;
    }
}
