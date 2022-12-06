package speciesnetwork.operators;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Operator;
import beast.base.util.Randomizer;
import speciesnetwork.EmbeddedTree;
import speciesnetwork.MultispeciesCoalescent;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;
import speciesnetwork.SanityChecks;
import speciesnetwork.simulator.CoalescentSimulator;

/**
 * Pick an internal network node randomly.
 * If bifurcation node: move the end of either the two parent branches.
 * If reticulation node: move the top of either the two child branches.
 * Propose a destination branch randomly from all possible branches, and attach the picked branch to the new position.
 * Also update the affected gene tree lineages at each locus to maintain the compatibility.
 *
 * @author Chi Zhang
 */

@Deprecated  // not working properly
@Description("Relocate the source of an edge starting with speciation node, " +
             "or the destination of an edge ending with hybridization node." +
             "Also update the affected gene trees to maintain compatibility.")
public class CoordinatedRelocateBranch extends Operator {
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

        // calculate coalescent prob. of current gene trees in current species network
        final MultispeciesCoalescent MSNC = MSNCInput.get();
        double logProposalRatio = MSNC.coalescentProb();

        // start moving species network
        speciesNetwork.startEditing(this);

        // pick an internal node randomly
        final NetworkNode[] internalNodes = speciesNetwork.getInternalNodes();
        int rIndex = Randomizer.nextInt(internalNodes.length);
        final NetworkNode pN = internalNodes[rIndex];

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
            List<Integer> candidateBrNrs = new ArrayList<>();
            for (NetworkNode node : speciesNetwork.getAllNodesExceptOrigin()) {
                if (node != pN && (pP.isSpeciation() || node.getHeight() < pP.getHeight())) {
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
            else if (pP.isReticulation()) {
                // cannot move in this case
                return Double.NEGATIVE_INFINITY;
            }
            else {
                // pP is a bifurcation node and the reticulation direction is changed
                pN.setHeight(pP.getHeight());
                pP.setHeight(newHeight);
                // determine nodes around pP
                final Integer pPpPPBranchNr = pP.gammaBranchNumber;
                final NetworkNode pPP = pP.getParentByBranch(pPpPPBranchNr); // pPP: parent of pP
                final Integer pPpPCBranchNr;
                if (pP.childBranchNumbers.get(0).equals(pickedBranchNr))
                    pPpPCBranchNr = pP.childBranchNumbers.get(1);
                else
                    pPpPCBranchNr = pP.childBranchNumbers.get(0);
                final NetworkNode pPC = pP.getChildByBranch(pPpPCBranchNr);  // pPC: another child of pP

                pP.childBranchNumbers.remove(pPpPCBranchNr);
                pN.childBranchNumbers.add(pPpPCBranchNr);
                if (pP == aC) {
                    // the attaching branch is pP-pPP
                    pP.childBranchNumbers.add(pNpNPBranchNr);
                } else {
                    pPP.childBranchNumbers.remove(pPpPPBranchNr);
                    pPP.childBranchNumbers.add(pNpNPBranchNr);
                    // separate aC and aP by pP
                    aP.childBranchNumbers.remove(attachBranchNr);
                    aP.childBranchNumbers.add(pPpPPBranchNr);
                    pP.childBranchNumbers.add(attachBranchNr);
                    pPP.updateRelationships();
                    aP.updateRelationships();
                    aC.updateRelationships();
                }
                // speciesNetwork.updateRelationships();
                pP.updateRelationships();
                pN.updateRelationships();
                pPC.updateRelationships();
            }

            logProposalRatio += Math.log(upper - lower) - Math.log(bounds);
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
            List<Integer> candidateBrNrs = new ArrayList<>();
            for (NetworkNode node : speciesNetwork.getInternalNodesWithOrigin()) {
                if (pC.isReticulation() || node.getHeight() > pC.getHeight()) {
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
            else if (!pC.isReticulation()) {
                // cannot move in this case
                return Double.NEGATIVE_INFINITY;
            }
            else {
                // pC is a reticulation node and the reticulation direction is changed
                pN.setHeight(pC.getHeight());
                pC.setHeight(newHeight);
                // determine nodes around pC
                final Integer pCpCCBranchNr = pC.childBranchNumbers.get(0);
                final NetworkNode pCC = pC.getChildByBranch(pCpCCBranchNr);  // pCC: child of pC
                final Integer pCpCPBranchNr;
                if (pC.gammaBranchNumber.equals(pickedBranchNr))
                    pCpCPBranchNr = pC.gammaBranchNumber + 1;
                else
                    pCpCPBranchNr = pC.gammaBranchNumber;
                final NetworkNode pCP = pC.getParentByBranch(pCpCPBranchNr); // pCP: another parent of pC

                pCP.childBranchNumbers.remove(pCpCPBranchNr);
                pCP.childBranchNumbers.add(pNpPBranchNr);
                if (aP == pC) {
                    // the attaching branch is pC-pCC
                    pN.childBranchNumbers.add(pCpCPBranchNr);
                } else {
                    pN.childBranchNumbers.add(pCpCCBranchNr);
                    // separate aC and aP by pC
                    aP.childBranchNumbers.remove(attachBranchNr);
                    aP.childBranchNumbers.add(pCpCPBranchNr);
                    pC.childBranchNumbers.remove(pCpCCBranchNr);
                    pC.childBranchNumbers.add(attachBranchNr);
                    // speciesNetwork.updateRelationships();
                    pCC.updateRelationships();
                    aP.updateRelationships();
                    aC.updateRelationships();
                }
                pCP.updateRelationships();
                pN.updateRelationships();
                pC.updateRelationships();
            }

            logProposalRatio += Math.log(upper - lower) - Math.log(bounds);
        }

        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());

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
