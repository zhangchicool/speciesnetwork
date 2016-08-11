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
 * TODO: narrow move that doesn't change the height
 * @author Chi Zhang
 */

@Description("Relocate the source of an edge starting with speciation node, " +
             "or the destination of an edge ending with hybridization node.")
public class EdgeRelocator extends Operator {
    public Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);
    public Input<List<RebuildEmbedding>> rebuildEmbeddingInput = new Input<>("rebuildEmbedding",
            "Operator which rebuilds embedding of gene tree within species network.", new ArrayList<>());
    public Input<Boolean> isNarrowInput =
            new Input<>("isNarrow", "If true, do not change the node height.", false);

    private final double lambda = 1.0;  // rate of exponential distribution

    // empty constructor to facilitate construction by XML + initAndValidate
    public EdgeRelocator() {
    }

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        Network speciesNetwork = speciesNetworkInput.get();
        SanityChecks.checkNetworkSanity(speciesNetwork.getRoot());

        // pick an internal node randomly
        final NetworkNode[] internalNodes = speciesNetwork.getInternalNodes();
        int randomIndex = Randomizer.nextInt(internalNodes.length);
        NetworkNode pickedNode = internalNodes[randomIndex];
        final int pickedNodeNr = pickedNode.getNr();

        // start moving
        speciesNetwork.startEditing(this);

        double logProposalRatio = 0.0;
        if (pickedNode.isReticulation()) {
            // move the end of either the two parent branches
            final int pickedBranchNr, alterBranchNr;
            if (Randomizer.nextBoolean()) {
                pickedBranchNr = pickedNode.gammaBranchNumber;
                alterBranchNr = pickedNode.gammaBranchNumber + 1;
            } else {
                pickedBranchNr = pickedNode.gammaBranchNumber + 1;
                alterBranchNr = pickedNode.gammaBranchNumber;
            }

            NetworkNode pickedParent = pickedNode.getParentByBranch(pickedBranchNr);
            final double upperLimit = pickedParent.getHeight();  // upper bound

            final int pickedChildBrNr = pickedNode.childBranchNumbers.get(0);
            NetworkNode pickedChild = pickedNode.getChildByBranch(pickedChildBrNr);

            NetworkNode alterParent = pickedNode.getParentByBranch(alterBranchNr);
            if (upperLimit > alterParent.getHeight())
                logProposalRatio -= Math.log(alterParent.getHeight() - pickedChild.getHeight());
            else
                logProposalRatio -= Math.log(upperLimit - pickedChild.getHeight());

            // look for all the candidate branches to attach to
            List<Integer> candidateBranchNumbers = new ArrayList<>();
            for (NetworkNode node: speciesNetwork.getAllNodes()) {
                // do not attach to the original position
                if (node != pickedNode && node != pickedChild && node.getHeight() < upperLimit) {
                    candidateBranchNumbers.add(node.gammaBranchNumber);
                    if (node.isReticulation())
                        candidateBranchNumbers.add(node.gammaBranchNumber + 1);
                }
            }

            // pick a candidate branch randomly
            if (candidateBranchNumbers.size() == 0)
                return Double.NEGATIVE_INFINITY;
            randomIndex = Randomizer.nextInt(candidateBranchNumbers.size());
            final int attachBranchNr = candidateBranchNumbers.get(randomIndex);

            final int attachChildNr = speciesNetwork.getNodeNumber(attachBranchNr);
            NetworkNode attachChild = speciesNetwork.getNode(attachChildNr);
            NetworkNode attachParent = attachChild.getParentByBranch(attachBranchNr);

            // propose an attachment height
            final double upper;
            if (upperLimit > attachParent.getHeight())
                upper = attachParent.getHeight();
            else
                upper = upperLimit;
            final double lower = attachChild.getHeight();
            final double newHeight = lower + (upper - lower) * Randomizer.nextDouble();
            logProposalRatio += Math.log(upper - lower);

            // set new height
            pickedNode.setHeight(newHeight);

            // deal with the node relationships
            pickedNode.deleteChild(pickedChild);
            pickedNode.deleteParent(alterParent);
            pickedNode.addChild(attachChild);
            pickedNode.addParent(attachParent);

            pickedChild.deleteParent(pickedNode);
            pickedChild.addParent(alterParent);

            alterParent.deleteChild(pickedNode);
            alterParent.addChild(pickedChild);

            attachChild.deleteParent(attachParent);
            attachChild.addParent(pickedNode);

            attachParent.deleteChild(attachChild);
            attachParent.addChild(pickedNode);

        } else {
            // move the top of either the two child branches
            final int pickedBranchNr, alterBranchNr;
            if (Randomizer.nextBoolean()) {
                pickedBranchNr = pickedNode.childBranchNumbers.get(0);
                alterBranchNr = pickedNode.childBranchNumbers.get(1);
            } else {
                pickedBranchNr = pickedNode.childBranchNumbers.get(1);
                alterBranchNr = pickedNode.childBranchNumbers.get(0);
            }

            NetworkNode pickedChild = pickedNode.getChildByBranch(pickedBranchNr);
            final double lowerLimit = pickedChild.getHeight();  // lower bound

            NetworkNode alterChild = pickedNode.getChildByBranch(alterBranchNr);

            final int pickedParentBrNr = pickedNode.gammaBranchNumber;
            NetworkNode pickedParent = pickedNode.getParentByBranch(pickedParentBrNr);
            if (pickedParent != null) {
                if (lowerLimit < alterChild.getHeight())
                    logProposalRatio -= Math.log(pickedParent.getHeight() - alterChild.getHeight());
                else
                    logProposalRatio -= Math.log(pickedParent.getHeight() - lowerLimit);
            }
            else {  // pickedParent is null when pickedNode is root
                if (lowerLimit < alterChild.getHeight())
                    logProposalRatio += Math.log(lambda) - lambda * (pickedNode.getHeight() - alterChild.getHeight());
                else
                    return Double.NEGATIVE_INFINITY;
            }

            // look for all the candidate branches to attach to
            List<Integer> candidateBranchNumbers = new ArrayList<>();
            for (NetworkNode node: speciesNetwork.getAllNodes()) {
                // do not attach to the original position
                if (node != pickedNode && node != alterChild) {
                    NetworkNode pNode = node.getParentByBranch(node.gammaBranchNumber);
                    if (pNode.getHeight() > lowerLimit && pNode != pickedNode)
                        candidateBranchNumbers.add(node.gammaBranchNumber);
                    if (node.isReticulation()) {
                        pNode = node.getParentByBranch(node.gammaBranchNumber + 1);
                        if (pNode.getHeight() > lowerLimit && pNode != pickedNode)
                            candidateBranchNumbers.add(node.gammaBranchNumber + 1);
                    }
                }
            }

            // pick a candidate branch randomly
            if (candidateBranchNumbers.size() == 0)
                return Double.NEGATIVE_INFINITY;
            randomIndex = Randomizer.nextInt(candidateBranchNumbers.size());
            final int attachBranchNr = candidateBranchNumbers.get(randomIndex);

            final int attachChildNr = speciesNetwork.getNodeNumber(attachBranchNr);
            NetworkNode attachChild = speciesNetwork.getNode(attachChildNr);
            NetworkNode attachParent = attachChild.getParentByBranch(attachBranchNr);

            // propose an attachment height
            final double newHeight;
            if (attachParent != null) {
                final double lower;
                if (lowerLimit < attachChild.getHeight())
                    lower = attachChild.getHeight();
                else
                    lower = lowerLimit;
                final double upper = attachParent.getHeight();
                newHeight = lower + (upper - lower) * Randomizer.nextDouble();
                logProposalRatio += Math.log(upper - lower);
            } else {  // attachParent is null when attachChild is root
                final double length = Randomizer.nextExponential(lambda);
                newHeight = attachChild.getHeight() + length;
                logProposalRatio -= Math.log(lambda) - lambda * (length);
            }

            // set new height
            pickedNode.setHeight(newHeight);

            // deal with the node relationships
            pickedNode.deleteChild(alterChild);
            pickedNode.deleteParent(pickedParent);
            pickedNode.addChild(attachChild);
            pickedNode.addParent(attachParent);

            alterChild.deleteParent(pickedNode);
            alterChild.addParent(pickedParent);

            if (pickedParent != null) {
                pickedParent.deleteChild(pickedNode);
                pickedParent.addChild(alterChild);
            }

            attachChild.deleteParent(attachParent);
            attachChild.addParent(pickedNode);

            if (attachParent != null) {
                attachParent.deleteChild(attachChild);
                attachParent.addChild(pickedNode);
            } else {
                // need to switch root
                speciesNetwork.swapRoot(pickedNodeNr);
            }
        }

        return logProposalRatio;
    }
}
