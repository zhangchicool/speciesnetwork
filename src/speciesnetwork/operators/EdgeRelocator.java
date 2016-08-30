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
    private final boolean isNarrow = isNarrowInput.get();

    // empty constructor to facilitate construction by XML + initAndValidate
    public EdgeRelocator() {
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
            final int pickedParentBrNr, otherParentBrNr;
            if (Randomizer.nextBoolean()) {
                pickedParentBrNr = pickedNode.gammaBranchNumber;
                otherParentBrNr = pickedNode.gammaBranchNumber + 1;
            } else {
                pickedParentBrNr = pickedNode.gammaBranchNumber + 1;
                otherParentBrNr = pickedNode.gammaBranchNumber;
            }

            NetworkNode pickedParent = pickedNode.getParentByBranch(pickedParentBrNr);
            final double upperLimit = pickedParent.getHeight();  // upper bound

            final int pickedChildBrNr = pickedNode.childBranchNumbers.get(0);
            NetworkNode pickedChild = pickedNode.getChildByBranch(pickedChildBrNr);

            NetworkNode otherParent = pickedNode.getParentByBranch(otherParentBrNr);

            // look for all the candidate branches to attach to
            List<Integer> candidateBranchNumbers = new ArrayList<>();
            if (isNarrow) {
                for (NetworkNode node : speciesNetwork.getAllNodes()) {
                    // do not attach to the original position
                    if (node != pickedNode && node != pickedChild && node.getHeight() < upperLimit) {
                        NetworkNode gParent = node.getParentByBranch(node.gammaBranchNumber);
                        NetworkNode oParent = node.getParentByBranch(node.gammaBranchNumber + 1);
                        if (gParent.getHeight() > pickedNode.getHeight())
                            candidateBranchNumbers.add(node.gammaBranchNumber);
                        if (node.isReticulation() && oParent.getHeight() > pickedNode.getHeight())
                            candidateBranchNumbers.add(node.gammaBranchNumber + 1);
                    }
                }
            } else {
                for (NetworkNode node : speciesNetwork.getAllNodes()) {
                    // do not attach to the original position
                    if (node != pickedNode && node != pickedChild && node.getHeight() < upperLimit) {
                        candidateBranchNumbers.add(node.gammaBranchNumber);
                        if (node.isReticulation())
                            candidateBranchNumbers.add(node.gammaBranchNumber + 1);
                    }
                }

                if (upperLimit > otherParent.getHeight())
                    logProposalRatio -= Math.log(otherParent.getHeight() - pickedChild.getHeight());
                else
                    logProposalRatio -= Math.log(upperLimit - pickedChild.getHeight());
            }

            // pick a candidate branch randomly
            if (candidateBranchNumbers.size() == 0)
                return Double.NEGATIVE_INFINITY;
            randomIndex = Randomizer.nextInt(candidateBranchNumbers.size());
            final int attachBranchNr = candidateBranchNumbers.get(randomIndex);

            final int attachChildNr = speciesNetwork.getNodeNumber(attachBranchNr);
            NetworkNode attachChild = speciesNetwork.getNode(attachChildNr);
            NetworkNode attachParent = attachChild.getParentByBranch(attachBranchNr);

            if (!isNarrow) {
                final double upper;
                if (upperLimit > attachParent.getHeight())
                    upper = attachParent.getHeight();
                else
                    upper = upperLimit;
                final double lower = attachChild.getHeight();
                // propose an attachment height
                final double newHeight = lower + (upper - lower) * Randomizer.nextDouble();
                // set new height
                pickedNode.setHeight(newHeight);

                logProposalRatio += Math.log(upper - lower);
            }

            // deal with the node relationships
            otherParent.childBranchNumbers.remove(otherParentBrNr);
            otherParent.childBranchNumbers.add(pickedChildBrNr);
            attachParent.childBranchNumbers.remove(attachBranchNr);
            attachParent.childBranchNumbers.add(otherParentBrNr);
            pickedNode.childBranchNumbers.remove(pickedChildBrNr);
            pickedNode.childBranchNumbers.add(attachBranchNr);

            otherParent.updateChildren();
            pickedChild.updateParents();
            attachParent.updateChildren();
            attachChild.updateParents();
            pickedNode.updateParents();
            pickedNode.updateChildren();

        } else {
            // move the top of either the two child branches
            final int pickedChildBrNr, otherChildBrNr;
            if (Randomizer.nextBoolean()) {
                pickedChildBrNr = pickedNode.childBranchNumbers.get(0);
                otherChildBrNr = pickedNode.childBranchNumbers.get(1);
            } else {
                pickedChildBrNr = pickedNode.childBranchNumbers.get(1);
                otherChildBrNr = pickedNode.childBranchNumbers.get(0);
            }

            NetworkNode pickedChild = pickedNode.getChildByBranch(pickedChildBrNr);
            final double lowerLimit = pickedChild.getHeight();  // lower bound

            NetworkNode otherChild = pickedNode.getChildByBranch(otherChildBrNr);

            final int pickedParentBrNr = pickedNode.gammaBranchNumber;
            NetworkNode pickedParent = pickedNode.getParentByBranch(pickedParentBrNr);

            // look for all the candidate branches to attach to
            List<Integer> candidateBranchNumbers = new ArrayList<>();
            if (isNarrow) {
                for (NetworkNode node : speciesNetwork.getAllNodes()) {
                    // do not attach to the original position
                    if (node != pickedNode && node.getHeight() > lowerLimit) {
                        for (int childBrNr : node.childBranchNumbers) {
                            NetworkNode gChild = node.getChildByBranch(childBrNr);
                            if (gChild != pickedNode && gChild.getHeight() < pickedNode.getHeight())
                                candidateBranchNumbers.add(childBrNr);
                        }
                    }
                }
            } else {
                for (NetworkNode node : speciesNetwork.getAllNodes()) {
                    // do not attach to the original position
                    if (node != pickedNode && node.getHeight() > lowerLimit) {
                        for (int childBrNr : node.childBranchNumbers) {
                            if (node.getChildByBranch(childBrNr) != pickedNode)
                                candidateBranchNumbers.add(childBrNr);
                        }
                    }
                }

                if (pickedParent != null) {
                    if (lowerLimit < otherChild.getHeight())
                        logProposalRatio -= Math.log(pickedParent.getHeight() - otherChild.getHeight());
                    else
                        logProposalRatio -= Math.log(pickedParent.getHeight() - lowerLimit);
                } else {  // pickedParent is null when pickedNode is root
                    if (lowerLimit < otherChild.getHeight())
                        logProposalRatio += Math.log(lambda) - lambda * (pickedNode.getHeight() - otherChild.getHeight());
                    else
                        return Double.NEGATIVE_INFINITY;
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

            if (!isNarrow) {
                final double newHeight;
                // propose an attachment height
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
            }

            // deal with the node relationships
            if (pickedParent != null) {
                pickedParent.childBranchNumbers.remove(pickedParentBrNr);
                pickedParent.childBranchNumbers.add(otherChildBrNr);
                pickedParent.updateChildren();
            }
            if (attachParent != null) {
                attachParent.childBranchNumbers.remove(attachBranchNr);
                attachParent.childBranchNumbers.add(pickedParentBrNr);
                attachParent.updateChildren();
            }
            pickedNode.childBranchNumbers.remove(otherChildBrNr);
            pickedNode.childBranchNumbers.add(attachBranchNr);

            otherChild.updateParents();
            attachChild.updateParents();
            pickedNode.updateParents();
            pickedNode.updateChildren();

            // swap root
            if (pickedParent == null && attachParent != null) {
                speciesNetwork.swapRoot(attachParent.getNr());
            }
            else if (pickedParent != null && attachParent == null) {
                speciesNetwork.swapRoot(pickedParent.getNr());
            }
        }

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
