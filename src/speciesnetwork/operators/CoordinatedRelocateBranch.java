package speciesnetwork.operators;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.util.Randomizer;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;
import speciesnetwork.utils.SanityChecks;

/**
 * Pick an internal network node randomly.
 * If bifurcation node: move the end of either the two parent branches.
 * If reticulation node: move the top of either the two child branches.
 * Propose a destination branch randomly from all possible branches, and attach the picked branch to the new position.
 * Also update the affected gene tree lineages at each locus to maintain the compatibility.
 *
 * @author Chi Zhang
 */

@Description("Relocate the source of an edge starting with speciation node, " +
             "or the destination of an edge ending with hybridization node." +
             "Also update the affected gene tree lineages to maintain compatibility.")
public class CoordinatedRelocateBranch extends CoordinatedOperator {

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        final Network speciesNetwork = speciesNetworkInput.get();
        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());

        // pick an internal node randomly
        final NetworkNode[] internalNodes = speciesNetwork.getInternalNodes();
        int rIndex = Randomizer.nextInt(internalNodes.length);
        final NetworkNode pN = internalNodes[rIndex];

        // start moving
        speciesNetwork.startEditing(this);

        final double logProposalRatio;

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
            // determine nodes around pN                                  // pN: picked node
            final NetworkNode pP = pN.getParentByBranch(pickedBranchNr);  // pP: parent of pN by picked branch
            final NetworkNode pNP = pN.getParentByBranch(pNpNPBranchNr);  // pNP: another parent of pN
            final Integer pNpCBranchNr = pN.childBranchNumbers.get(0);
            final NetworkNode pC = pN.getChildByBranch(pNpCBranchNr);     // pC: (only) child of pN
            // upper and lower bounds for the backward move
            final double bounds = pNP.getHeight() - pC.getHeight();

            // join pC and pNP with a single branch
            pNP.childBranchNumbers.remove(pNpNPBranchNr);
            pNP.childBranchNumbers.add(pNpCBranchNr);
            pN.childBranchNumbers.remove(pNpCBranchNr);
            pNP.updateRelationships();
            pC.updateRelationships();

            // look for all the candidate branches to attach to
            // do not change the direction of picked branch at the moment
            List<Integer> candidateBrNrs = new ArrayList<>();
            for (NetworkNode node : speciesNetwork.getAllNodesExceptOrigin()) {
                if (node != pN && node.getHeight() < pP.getHeight()) {
                    candidateBrNrs.add(node.gammaBranchNumber);
                    if (node.isReticulation())
                        candidateBrNrs.add(node.gammaBranchNumber + 1);
                }
            }
            if (candidateBrNrs.size() == 0)
                return Double.NEGATIVE_INFINITY;

            // pick a candidate branch randomly
            rIndex = Randomizer.nextInt(candidateBrNrs.size());
            final Integer attachBranchNr = candidateBrNrs.get(rIndex);
            final int aCNodeNr = speciesNetwork.getNodeNumber(attachBranchNr);
            final NetworkNode aC = speciesNetwork.getNode(aCNodeNr);      // aC: child at attaching branch
            final NetworkNode aP = aC.getParentByBranch(attachBranchNr);  // aP: parent at attaching branch

            final double upper = aP.getHeight();
            final double lower = aC.getHeight();
            // propose an attachment height
            final double newHeight = lower + (upper - lower) * Randomizer.nextDouble();

            // relocate the picked node and branch
            if (newHeight < pP.getHeight()) {
                // the reticulation direction is not changed
                pN.setHeight(newHeight);
                // separate aC and aP by pN
                aP.childBranchNumbers.remove(attachBranchNr);
                aP.childBranchNumbers.add(pNpNPBranchNr);
                pN.childBranchNumbers.add(attachBranchNr);
                aP.updateRelationships();
                pN.updateRelationships();
                aC.updateRelationships();
            }
            else {
                // cannot move in this case
                return Double.NEGATIVE_INFINITY;
            }

            logProposalRatio = Math.log(upper - lower) - Math.log(bounds);

            // TODO: update gene trees
        }
        else {
            // move the top of either the two child branches
            final Integer pickedBranchNr, pNpNCBranchNr;
            if (Randomizer.nextBoolean()) {
                pickedBranchNr = pN.childBranchNumbers.get(0);
                pNpNCBranchNr = pN.childBranchNumbers.get(1);
            } else {
                pickedBranchNr = pN.childBranchNumbers.get(1);
                pNpNCBranchNr = pN.childBranchNumbers.get(0);
            }
            // determine nodes around pN                                  // pN: picked node
            final NetworkNode pC = pN.getChildByBranch(pickedBranchNr);   // pC: child of pN by picked branch
            final NetworkNode pNC = pN.getChildByBranch(pNpNCBranchNr);   // pNC: another child of pN
            final Integer pNpPBranchNr = pN.gammaBranchNumber;
            final NetworkNode pP = pN.getParentByBranch(pNpPBranchNr);    // pP: (only) parent of pN
            // upper and lower bounds for the backward move
            final double bounds = pP.getHeight() - pNC.getHeight();

            // join pP and pNC with a single branch
            pP.childBranchNumbers.remove(pNpPBranchNr);
            pP.childBranchNumbers.add(pNpNCBranchNr);
            pN.childBranchNumbers.remove(pNpNCBranchNr);
            pP.updateRelationships();
            pNC.updateRelationships();

            // look for all the candidate branches to attach to
            // do not change the direction of picked branch at the moment
            List<Integer> candidateBrNrs = new ArrayList<>();
            for (NetworkNode node : speciesNetwork.getInternalNodesWithOrigin()) {
                if (node.getHeight() > pC.getHeight()) {
                    for (Integer childBrNr : node.childBranchNumbers) {
                        if (!childBrNr.equals(pickedBranchNr) && !childBrNr.equals(pNpPBranchNr))
                            candidateBrNrs.add(childBrNr);
                    }
                }
            }
            if (candidateBrNrs.size() == 0)
                return Double.NEGATIVE_INFINITY;

            // pick a candidate branch randomly
            rIndex = Randomizer.nextInt(candidateBrNrs.size());
            final Integer attachBranchNr = candidateBrNrs.get(rIndex);
            final int aCNodeNr = speciesNetwork.getNodeNumber(attachBranchNr);
            final NetworkNode aC = speciesNetwork.getNode(aCNodeNr);      // aC: child at attaching branch
            final NetworkNode aP = aC.getParentByBranch(attachBranchNr);  // aP: parent at attaching branch

            final double upper = aP.getHeight();
            final double lower = aC.getHeight();
            // propose an attachment height
            final double newHeight = lower + (upper - lower) * Randomizer.nextDouble();

            // relocate the picked node and branch
            if (newHeight > pC.getHeight()) {
                pN.setHeight(newHeight);
                // separate aC and aP by pN
                aP.childBranchNumbers.remove(attachBranchNr);
                aP.childBranchNumbers.add(pNpPBranchNr);
                pN.childBranchNumbers.add(attachBranchNr);
                // speciesNetwork.updateRelationships();
                aP.updateRelationships();
                pN.updateRelationships();
                aC.updateRelationships();
            }
            else {
                // cannot move in this case
                return Double.NEGATIVE_INFINITY;
            }

            logProposalRatio = Math.log(upper - lower) - Math.log(bounds);

            // TODO: update gene trees
        }

        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());

        return logProposalRatio;
    }
}
