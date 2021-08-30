package speciesnetwork.operators;

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
        // number of reticulation branches in the current network
        final int nReticulationBranches = 2 * nHybridNodes;  // m'

        // pick a reticulation branch randomly
        final Integer hybridBranchNr = Randomizer.nextInt(nReticulationBranches) + speciesNetwork.getReticulationOffset();
        final int hybridNodeNr = speciesNetwork.getNodeNumber(hybridBranchNr);
        // branch with hybridBranchNr is connecting hybridNode and parentNode
        NetworkNode hybridNode = speciesNetwork.getNode(hybridNodeNr);
        NetworkNode parentNode = hybridNode.getParentByBranch(hybridBranchNr);
        if (parentNode.isReticulation())  // cannot delete a branch connecting two reticulation nodes
            return Double.NEGATIVE_INFINITY;

        // get the parent node and another child node of parentNode
        final Integer pNParentBranchNr = parentNode.gammaBranchNumber;
        NetworkNode pNParentNode = parentNode.getParentByBranch(pNParentBranchNr);
        final Integer pNChildBranchNr;
        if (parentNode.childBranchNumbers.get(0).equals(hybridBranchNr))
            pNChildBranchNr = parentNode.childBranchNumbers.get(1);
        else
            pNChildBranchNr = parentNode.childBranchNumbers.get(0);
        NetworkNode pNChildNode = parentNode.getChildByBranch(pNChildBranchNr);

        // get the child node and another parent node of hybridNode
        final Integer hNChildBranchNr = hybridNode.childBranchNumbers.get(0);
        NetworkNode hNChildNode = hybridNode.getChildByBranch(hNChildBranchNr);
        final Integer hNParentBranchNr;
        if (hybridNode.gammaBranchNumber.equals(hybridBranchNr))
            hNParentBranchNr = hybridNode.gammaBranchNumber + 1;
        else
            hNParentBranchNr = hybridNode.gammaBranchNumber;
        NetworkNode hNParentNode = hybridNode.getParentByBranch(hNParentBranchNr);

        // work out the Jacobian
        final double l1, l2;
        if (parentNode == hNParentNode && hybridNode == pNChildNode) {
            // the two attaching points are on the same branch
            l1 = l2 = pNParentNode.getHeight() - hNChildNode.getHeight();
        } else {
            // the two attaching points are on different branches
            l1 = hNParentNode.getHeight() - hNChildNode.getHeight();
            l2 = pNParentNode.getHeight() - pNChildNode.getHeight();
        }
        double logProposalRatio = - Math.log(l1) - Math.log(l2);

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
