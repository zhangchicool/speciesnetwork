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
 * Pick an internal network node randomly.
 * If bifurcation node: move the end of either the two parent branches.
 * If reticulation node: move the top of either the two child branches.
 * Propose a destination branch randomly from all possible branches, and attach the picked branch to the new position.
 *
 * @author Chi Zhang
 */

@Description("Relocate the source of an edge starting with speciation node, " +
             "or the destination of an edge ending with hybridization node.")
public class RelocateBranch extends Operator {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);
    public final Input<Boolean> isWideInput =
            new Input<>("isWide", "If true, change the node height (default is false).", false);

    private boolean isWide;

    // empty constructor to facilitate construction by XML + initAndValidate
    public RelocateBranch() {
    }

    @Override
    public void initAndValidate() {
        isWide = isWideInput.get();
    }

    @Override
    public double proposal() {
        final Network speciesNetwork = speciesNetworkInput.get();
        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());

        // number of branches in the network
        final int branchCount = speciesNetwork.getBranchCount();

        // pick an internal node randomly
        final NetworkNode[] internalNodes = speciesNetwork.getInternalNodes();
        final int rIndex = Randomizer.nextInt(internalNodes.length);
        NetworkNode pN = internalNodes[rIndex];

        // start moving
        speciesNetwork.startEditing(this);

        double logProposalRatio = 0.0;

        if (pN.isReticulation()) {
            // move the end of either the two parent branches
            final Integer pickedBranchNr, pNpNPBranchNr;
            if (Randomizer.nextBoolean()) {
                pickedBranchNr = pN.gammaBranchNumber;
                pNpNPBranchNr = pN.gammaBranchNumber + 1;
            } else {
                pickedBranchNr = pN.gammaBranchNumber + 1;
                pNpNPBranchNr = pN.gammaBranchNumber;
            }
            // determine nodes around
            NetworkNode pP = pN.getParentByBranch(pickedBranchNr);
            NetworkNode pNP = pN.getParentByBranch(pNpNPBranchNr);
            final Integer pNpCBranchNr = pN.childBranchNumbers.get(0);
            NetworkNode pC = pN.getChildByBranch(pNpCBranchNr);

            // pick a candidate branch randomly
            Integer attachBranchNr;
            do {
                attachBranchNr = Randomizer.nextInt(branchCount);
            } while (attachBranchNr.equals(pickedBranchNr) || attachBranchNr.equals(pNpNPBranchNr));
            final int aCNodeNr = speciesNetwork.getNodeNumber(attachBranchNr);
            NetworkNode aC = speciesNetwork.getNode(aCNodeNr);
            NetworkNode aP = aC.getParentByBranch(attachBranchNr);

            final double upper = aP.getHeight();
            final double lower = aC.getHeight();

            // propose an attachment height
            final double newHeight = lower + (upper - lower) * Randomizer.nextDouble();

            // deal with the node relationships
            if (pN.getHeight() < pP.getHeight()) {
                //
                
//                otherParent.childBranchNumbers.remove(otherParentBrNr);
//                otherParent.childBranchNumbers.add(childBranchNr);
//                attachParent.childBranchNumbers.remove(attachBranchNr);
//                attachParent.childBranchNumbers.add(otherParentBrNr);
//                pickedNode.childBranchNumbers.remove(childBranchNr);
//                pickedNode.childBranchNumbers.add(attachBranchNr);
//                otherParent.updateRelationships();
//                pickedChild.updateRelationships();
//                attachParent.updateRelationships();
//                attachChild.updateRelationships();
//                pickedNode.updateRelationships();
            } else {
                //


            }

        } else {}

        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());

        return logProposalRatio;
    }
}
