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
 * This proposal deletes a reticulation branch from the species network. If there is no reticulation, or if the branch
 * is connecting two reticulation nodes, this is aborted. The two branches at each connecting point are joined,
 * resulting branches with length l1 and l2 respectively. The gamma prob r is removed.
 * The Jacobian is 1 / (l1 * l2).
 *
 * The AddReticulation and DeleteReticulation are chosen with equal prob.
 * Let m be the number of reticulation branches in the current network. The probability of selecting the this branch to
 * remove is 1/m.
 * Let k be the number of branches in the proposed network. The probability of adding this branch is (1/k)(1/k)
 * The Hastings ratio is (1/k)(1/k)(1)(1)(1) / (1/m) = m / k^2.
 *
 * See also AddReticulation.
 *
 * @author Chi Zhang
 */

@Description("Delete a reticulation branch from the species network.")
public class DeleteReticulation extends Operator {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);

    // empty constructor to facilitate construction by XML + initAndValidate
    public DeleteReticulation() {
    }

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        final Network speciesNetwork = speciesNetworkInput.get();
        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());

        final int nHybridNodes = speciesNetwork.getReticulationNodeCount();
        if (nHybridNodes == 0)  // there is no reticulation branch to delete
            return Double.NEGATIVE_INFINITY;
        //number of reticulation branches in the current network
        final int nReticulationBranches = 2 * nHybridNodes;  // m'

        // pick a reticulation branch randomly
        final Integer hybridBranchNr = Randomizer.nextInt(nReticulationBranches) + speciesNetwork.getReticulationOffset();
        final int hybridNodeNr = speciesNetwork.getNodeNumber(hybridBranchNr);
        // branch hybridBranchNr is connecting hybridNode and pickedParent
        NetworkNode hybridNode = speciesNetwork.getNode(hybridNodeNr);
        NetworkNode pickedParent = hybridNode.getParentByBranch(hybridBranchNr);
        if (pickedParent.isReticulation())  // cannot delete a branch connecting two reticulation nodes
            return Double.NEGATIVE_INFINITY;

        // get the child node and another parent node of hybridNode
        final Integer childBranchNr1 = hybridNode.childBranchNumbers.get(0);
        NetworkNode childNode1 = hybridNode.getChildByBranch(childBranchNr1);
        final Integer parentBranchNr1;
        if (hybridNode.gammaBranchNumber.equals(hybridBranchNr)) {
            parentBranchNr1 = hybridNode.gammaBranchNumber + 1;
        } else {
            parentBranchNr1 = hybridNode.gammaBranchNumber;
        }
        NetworkNode parentNode1 = hybridNode.getParentByBranch(parentBranchNr1);

        // get the parent node and another child node of pickedParent
        final Integer childBranchNr2;
        if (pickedParent.childBranchNumbers.get(0).equals(hybridBranchNr)) {
            childBranchNr2 = pickedParent.childBranchNumbers.get(1);
        } else {
            childBranchNr2 = pickedParent.childBranchNumbers.get(0);
        }
        NetworkNode childNode2 = pickedParent.getChildByBranch(childBranchNr2);
        final Integer parentBranchNr2 = pickedParent.gammaBranchNumber;
        NetworkNode parentNode2 = pickedParent.getParentByBranch(parentBranchNr2);

        double logProposalRatio = 0.0;

        // work out the Jacobian
        final double l1, l2;
        if (hybridNode == childNode2 && pickedParent == parentNode1) {
            // the two attaching points are on the same branch
            l1 = l2 = parentNode2.getHeight() - childNode1.getHeight();
        } else {
            // the two attaching points are on different branches
            l1 = parentNode1.getHeight() - childNode1.getHeight();
            l2 = parentNode2.getHeight() - childNode2.getHeight();
        }
        logProposalRatio -= Math.log(l1) + Math.log(l2);  // the Jacobian

        // start moving
        speciesNetwork.startEditing(this);

        // delete the reticulation branch
        speciesNetwork.deleteReticulationBranch(hybridBranchNr);

        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());

        // number of branches in the proposed network
        final int nBranches = speciesNetwork.getBranchCount();  // k'
        logProposalRatio += Math.log(nReticulationBranches) - 2 * Math.log(nBranches);

        return logProposalRatio;
    }
}
