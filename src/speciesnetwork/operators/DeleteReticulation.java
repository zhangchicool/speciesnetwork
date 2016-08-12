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
 * This proposal delete a reticulation branch from the species network. If there is no reticulation, this is aborted.
 * The two branches at each connecting point are joined, resulting branches with length l1 and l2 respectively.
 * If the deleted node is root, its child becomes the new root. The gamma prob r is removed.
 * The Jacobian is 1/(l1*l2). (li = 1 if the deleted node is root. If the root branch is picked twice, the Jacobian = 1)
 *
 * The AddReticulation and DeleteReticulation are chosen with equal prob.
 * Let m' be the number of reticulation branches in the current network. The probability of selecting the this branch to
 * remove is (1/m').
 * Let k' be the number of branches in the proposed network. The probability of adding this branch is (1/k')(1/k')

 * The Hastings ratio is (1/k')(1/k')(g1)(g2)(g3) / (1/m') = m' * g1 * g2 * g3 / k'^2.
 * g1 = b * exp(-b * l11) if root, otherwise g1 = 1 (uniform density); g2 = b * exp(-b * l21) if root, otherwise g2 = 1.
 * g3 = 1. See also AddReticulation.
 *
 * @author Chi Zhang
 */

@Description("Delete a reticulation branch from the species network.")
public class DeleteReticulation extends Operator {
    public Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);
    public Input<List<RebuildEmbedding>> rebuildEmbeddingInput = new Input<>("rebuildEmbedding",
            "Operator which rebuilds embedding of gene tree within species network.", new ArrayList<>());

    private final double lambda = 1.0;  // rate of exponential distribution

    // empty constructor to facilitate construction by XML + initAndValidate
    public DeleteReticulation() {
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

        final int nHybridNodes = speciesNetwork.getReticulationNodeCount();

        //number of reticulation branches in the current network
        final int nReticulationBranches = 2 * nHybridNodes;  // m'

        // pick a reticulation branch randomly
        final int hybridNodeNr = Randomizer.nextInt(nHybridNodes) + speciesNetwork.getReticulationOffset();
        NetworkNode hybridNode = speciesNetwork.getNode(hybridNodeNr);
        final int hybridBranchNr;
        if (Randomizer.nextBoolean()) {
            hybridBranchNr = hybridNode.gammaBranchNumber;
        } else {
            hybridBranchNr = hybridNode.gammaBranchNumber + 1;
        }
        // branch hybridBranchNr is connecting hybridNode and pickedParent
        NetworkNode pickedParent = hybridNode.getParentByBranch(hybridBranchNr);

        // get the child node and another parent node of hybridNode
        final int childBranchNr1 = hybridNode.childBranchNumbers.get(0);
        NetworkNode childNode1 = hybridNode.getChildByBranch(childBranchNr1);

        NetworkNode parentNode1;
        if (hybridBranchNr == hybridNode.gammaBranchNumber) {
            parentNode1 = hybridNode.getParentByBranch(hybridBranchNr + 1);
        } else {
            parentNode1 = hybridNode.getParentByBranch(hybridBranchNr);
        }

        // get the parent node and another child node of pickedParent
        final int childBranchNr2;
        if (hybridBranchNr == pickedParent.childBranchNumbers.get(0)) {
            childBranchNr2 = pickedParent.childBranchNumbers.get(1);
        } else {
            childBranchNr2 = pickedParent.childBranchNumbers.get(0);
        }
        NetworkNode childNode2 = pickedParent.getChildByBranch(childBranchNr2);

        final int parentBranchNr2 = pickedParent.gammaBranchNumber;
        NetworkNode parentNode2 = pickedParent.getParentByBranch(parentBranchNr2);  // null if pickedParent is root

        final double l1, l2, l11, l21;
        double logProposalRatio = 0.0;

        // work out the Jacobian
        if (hybridNode == childNode2 && pickedParent == parentNode1) {
            // the two attaching points are on the same branch
            if (pickedParent.isRoot()) {
                l1 = l2 = 1;
                l11 = hybridNode.getHeight() - childNode1.getHeight();
                l21 = pickedParent.getHeight() - childNode1.getHeight();
                logProposalRatio += 2 * Math.log(lambda) - lambda * (l11 + l21);
            } else {
                l1 = l2 = parentNode2.getHeight() - childNode1.getHeight();
            }
        } else {
            // the two attaching points are on different branches
            l1 = parentNode1.getHeight() - childNode1.getHeight();
            if (pickedParent.isRoot()) {
                l2 = 1;
                l21 = pickedParent.getHeight() - childNode2.getHeight();
                logProposalRatio += Math.log(lambda) - lambda * l21;
            } else {
                l2 = parentNode2.getHeight() - childNode2.getHeight();
            }
        }
        logProposalRatio += - Math.log(l1) - Math.log(l2);  // the Jacobian

        // start moving
        speciesNetwork.startEditing(this);

        // delete the reticulation branch
        speciesNetwork.deleteReticulationBranch(hybridBranchNr);

        // number of branches in the proposed network
        final int nBranches = speciesNetwork.getBranchCount();  // k'
        logProposalRatio += Math.log(nReticulationBranches) - 2 * Math.log(nBranches);

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
