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
 * This proposal adds a reticulation event by connecting two existing branches (with length l1 and l2) with a new branch.
 * The same branch can be picked twice (and forms a loop to that branch). The cutting proportion of each picked branch by
 * the connecting point, w1 and w2 ~ Uniform(0,1). Let l11 = l1 * w1, l12 = l1 * (1-w1), l21 = l2 * w2, l22 = l2 * (1-w2)
 * The direction of the new branch is determined by the two connecting points, the higher is speciation node, and the
 * lower is reticulation node. The gamma prob r = w3 ~ Uniform(0,1).
 * The Jacobian is l1 * l2.
 *
 * The AddReticulation and DeleteReticulation are chosen with equal prob. If there is no reticulation in the network,
 * the DeleteReticulation move is aborted.
 * Let k be the number of branches in the current network. The probability of adding this branch is (1/k)(1/k)
 * Let m be the number of reticulation branches in the proposed network. The probability of selecting the same branch to
 * remove is (1/m).
 * The Hastings ratio is (1/m) / [(1/k)(1/k)(g1)(g2)(g3)] = k^2 / m, with g1 = g2 = g3 = 1 (uniform density).
 *
 * See also DeleteReticulation.
 *
 * @author Chi Zhang
 */

@Deprecated  // not working properly
@Description("Add a reticulation branch to the species network." +
             "Also update the affected gene trees to maintain compatibility.")
public class CoordinatedAddReticulation extends Operator {
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

        // number of branches in the current network
        final int nBranches = speciesNetwork.getBranchCount();  // k

        // pick two branches randomly, including the root branch
        final Integer pickedBranchNr1 = Randomizer.nextInt(nBranches);
        final Integer pickedBranchNr2 = Randomizer.nextInt(nBranches);  // allow picking the same branch

        // get the nodes associated with each branch
        final int pickedNodeNr1 = speciesNetwork.getNodeNumber(pickedBranchNr1);
        NetworkNode pickedNode1 = speciesNetwork.getNode(pickedNodeNr1);
        final int pickedNodeNr2 = speciesNetwork.getNodeNumber(pickedBranchNr2);
        NetworkNode pickedNode2 = speciesNetwork.getNode(pickedNodeNr2);
        NetworkNode pickedParent1 = pickedNode1.getParentByBranch(pickedBranchNr1);
        NetworkNode pickedParent2 = pickedNode2.getParentByBranch(pickedBranchNr2);

        // propose the attaching position at each branch
        final double l1, l2, l11, l21;
        l1 = pickedParent1.getHeight() - pickedNode1.getHeight();
        l11 = l1 * Randomizer.nextDouble();
        l2 = pickedParent2.getHeight() - pickedNode2.getHeight();
        l21 = l2 * Randomizer.nextDouble();

        double logProposalRatio = Math.log(l1) + Math.log(l2);  // the Jacobian

        // calculate coalescent prob. of current gene trees in current species network
        final MultispeciesCoalescent MSNC = MSNCInput.get();
        logProposalRatio += MSNC.coalescentProb();

        // start moving species network
        speciesNetwork.startEditing(this);

        // create two new nodes
        NetworkNode middleNode1 = new NetworkNode(speciesNetwork);
        NetworkNode middleNode2 = new NetworkNode(speciesNetwork);
        // set height
        middleNode1.setHeight(pickedNode1.getHeight() + l11);
        middleNode2.setHeight(pickedNode2.getHeight() + l21);

        // add a branch joining the two middle nodes (picked branches)
        if (middleNode1.getHeight() < middleNode2.getHeight()) {
            speciesNetwork.addReticulationBranch(middleNode1, middleNode2, pickedBranchNr1, pickedBranchNr2);
            middleNode1.setGammaProb(Randomizer.nextDouble());
        } else {
            speciesNetwork.addReticulationBranch(middleNode2, middleNode1, pickedBranchNr2, pickedBranchNr1);
            middleNode2.setGammaProb(Randomizer.nextDouble());
        }

        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());

        // number of reticulation branches in the proposed network
        final int nReticulationBranches = 2 * speciesNetwork.getReticulationNodeCount();  // m
        logProposalRatio += 2 * Math.log(nBranches) - Math.log(nReticulationBranches);

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
