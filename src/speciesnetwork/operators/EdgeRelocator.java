package speciesnetwork.operators;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.IntegerParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import speciesnetwork.Network;
import speciesnetwork.NetworkNode;

/**
 * @author Chi Zhang
 */

@Description("Relocate the source of an edge starting with speciation node, " +
        "or the destination of an edge ending with hybridization node.")
public class EdgeRelocator extends Operator {
    public Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);
    public Input<List<Tree>> geneTreesInput =
            new Input<>("geneTree", "list of gene trees embedded in species network", new ArrayList<>());
    public Input<List<IntegerParameter>> embeddingsInput =
            new Input<>("embedding", "The matrices to embed the gene trees in the species network.", new ArrayList<>());
    public Input<TaxonSet> taxonSuperSetInput =
            new Input<>("taxonSuperset", "Super-set of taxon sets mapping lineages to species.", Validate.REQUIRED);

    // empty constructor to facilitate construction by XML + initAndValidate
    public EdgeRelocator() {
    }

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        Network speciesNetwork = speciesNetworkInput.get();

        // pick an internal node randomly
        final NetworkNode[] internalNodes = speciesNetwork.getInternalNodes();
        int randomIndex = Randomizer.nextInt(internalNodes.length);
        NetworkNode pickedNode = internalNodes[randomIndex];
        final int pickedNodeNr = pickedNode.getNr();

        final double newInterval, oldInterval;  // for calculating proposal ratio

        // start moving
        speciesNetwork.startEditing(this);

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
            if (pickedParent == null)
                throw new RuntimeException("Developer ERROR: picked node has no parent at this branch!");
            final double upperLimit = pickedParent.getHeight();  // upper bound

            NetworkNode pickedChild = null;
            for(Integer i: pickedNode.childBranchNumbers) {
                pickedChild = pickedNode.getChildByBranch(i);
            }
            if (pickedChild == null)
                throw new RuntimeException("Developer ERROR: picked node has no child!");

            NetworkNode alterParent = pickedNode.getParentByBranch(alterBranchNr);
            if (upperLimit > alterParent.getHeight())
                oldInterval = alterParent.getHeight() - pickedChild.getHeight();
            else
                oldInterval = upperLimit - pickedChild.getHeight();

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
            newInterval = upper - lower;

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


        }

        return 0.0; //Math.log(newInterval/oldInterval);
    }
}
