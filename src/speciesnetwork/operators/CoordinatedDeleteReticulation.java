package speciesnetwork.operators;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.util.Randomizer;
import speciesnetwork.EmbeddedTree;
import speciesnetwork.MultispeciesCoalescent;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;
import speciesnetwork.SanityChecks;
import speciesnetwork.simulator.CoalescentSimulator;

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

@Description("Delete a reticulation branch from the species network." +
             "Also update the affected gene trees to maintain compatibility.")
public class CoordinatedDeleteReticulation extends Operator {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);
    public final Input<List<EmbeddedTree>> geneTreesInput = new Input<>("geneTree",
            "The gene tree within the species network.", new ArrayList<>());
    public final Input<MultispeciesCoalescent> MSNCInput =
            new Input<>("MSNC", "The multispecies network coalescent.", Validate.REQUIRED);
    public final Input<CoalescentSimulator> coalSimulatorInput = new Input<>("coalescentSimulator",
            "Simulate gene trees given the species network.", Validate.REQUIRED);

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

        // calculate coalescent prob. of current gene trees in current species network
        final MultispeciesCoalescent MSNC = MSNCInput.get();
        logProposalRatio += MSNC.coalescentProb();

        // start moving species network
        speciesNetwork.startEditing(this);

        // delete the reticulation branch
        speciesNetwork.deleteReticulationBranch(hybridBranchNr);

        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());

        // number of branches in the proposed network
        final int nBranches = speciesNetwork.getBranchCount();  // k'
        logProposalRatio += Math.log(nReticulationBranches) - 2 * Math.log(nBranches);

        // update gene trees (simulate random gene trees under MSNC)
        final List<EmbeddedTree> geneTrees = geneTreesInput.get();
        for (EmbeddedTree geneTree : geneTrees) {
            geneTree.startEditing(this);  // *all* gene trees will be edited
        }
        CoalescentSimulator geneTreesSimulator = coalSimulatorInput.get();
        geneTreesSimulator.simulate();
        geneTreesSimulator.popSizesInput.get().getCurrentEditable(this);  // hack to let state&node properly stored

        // calculate coalescent prob. of new gene trees in new species network
        // do NOT call MSNC.calculateLogP(); doing that would update 'logP' which should not be changed at this stage
        logProposalRatio -= MSNC.coalescentProb();

        return logProposalRatio;
    }
}
