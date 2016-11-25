package speciesnetwork.simulator;

import java.io.*;
import java.util.List;
import java.util.ArrayList;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Runnable;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.util.Randomizer;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;
import speciesnetwork.tools.SummarizeNetwork;

/**
 * @author Chi Zhang
 */

@Description("Simulate a species network under the pure birth and hybridization process.")
public class YuleHybridSimulator extends Runnable {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "Species network to be simulated.", Validate.REQUIRED);
    public final Input<RealParameter> originInput =
            new Input<>("origin", "The time when the process started.", Validate.REQUIRED);
    public final Input<RealParameter> tmrcaInput =
            new Input<>("tmrca", "The time of most recent common ancestor (root).", Validate.XOR, originInput);

    public final Input<RealParameter> birthRateInput =
            new Input<>("birthRate", "Speciation rate, lambda.", Validate.REQUIRED);
    public final Input<RealParameter> hybridRateInput =
            new Input<>("hybridRate", "Hybridization rate, nu.", Validate.REQUIRED);
    public final Input<Integer> nHybridizationInput = new Input<>("nHybrid", "Number of Hybridization.");

    public final Input<String> outputFileNameInput =
            new Input<>("outputFileName", "If provided, write to this file rather than to standard out.");
    public final Input<Integer> iterationsInput =
            new Input<>("iterations", "Number of iterations to simulate (default is 1).");

    public final Input<SummarizeNetwork> summaryInput =
            new Input<>("summary", "Print summary statistics of the networks.");

    @Override
    public void initAndValidate() {
    }

    @Override
    public void run() throws IOException {
        String outputFileName = outputFileNameInput.get();

        PrintStream out;  // where to print
        if (outputFileName == null) {
            out = System.out;
        } else {
            String msg = "Writing";
            if (new File(outputFileName).exists())
                msg = "Warning: Overwriting";
            System.err.println(msg + " file " + outputFileName);
            out = new PrintStream(outputFileName);
        }
        // print header
        out.println("#NEXUS\n\nBegin trees;");

        final int nrOfIterations;
        if (iterationsInput.get() == null)
            nrOfIterations = 1;
        else
            nrOfIterations = iterationsInput.get();
        for (int iteration = 1; iteration <= nrOfIterations; iteration++) {
            out.println("tree SIM_" + iteration + " =" + simulate().toString() + ";");
        }
        out.println("End;");

        if (summaryInput.get() != null) {
            summaryInput.get().run();
        }
    }

    protected Network simulate() {
        Network speciesNetwork = speciesNetworkInput.get();
        TaxonSet speciesTaxa = speciesNetwork.taxonSetInput.get();
        List<String> speciesNames = new ArrayList<>();

        final int numTips;  // number of extant species (-1 for no such condition)
        if (speciesTaxa != null) {
            speciesNames = speciesTaxa.asStringList();
            numTips = speciesNames.size();
        } else {
            numTips = -1;
        }
        final int numHybrid; // number of hybridization (-1 for no such condition)
        if (nHybridizationInput.get() != null)
            numHybrid = nHybridizationInput.get();
        else
            numHybrid = -1;

        if (numTips < 0) {
            simulate(speciesNetwork);
        }
        else if (numHybrid < 0) {  // condition on numTips
            do {  // simulate until we get the desired condition (caution!)
               simulate(speciesNetwork);
            } while (speciesNetwork.getLeafNodeCount() != numTips);
            // set the tip labels to match the taxa labels
            for (int i = 0; i < numTips; i++) {
                NetworkNode leaf = speciesNetwork.getNode(i);
                leaf.setLabel(speciesNames.get(i));
            }
        }
        else {  // condition on numTips and numHybrid, and no bubble in this case
            do {  // simulate until we get the desired condition (caution!)
                simulate(speciesNetwork);
            } while (speciesNetwork.getLeafNodeCount() != numTips ||
                    speciesNetwork.getReticulationNodeCount() != numHybrid || speciesNetwork.hasBubble());
            // set the tip labels to match the taxa labels
            for (int i = 0; i < numTips; i++) {
                NetworkNode leaf = speciesNetwork.getNode(i);
                leaf.setLabel(speciesNames.get(i));
            }
        }

        return speciesNetwork;
    }

    private void simulate(Network speciesNetwork) {
        final double lambda = birthRateInput.get().getValue();
        final double nu = hybridRateInput.get().getValue();

        final double timeOrigin;
        if (originInput.get() != null)
            timeOrigin = originInput.get().getValue();
        else
            timeOrigin = tmrcaInput.get().getValue() + 1.0 / lambda;

        // set the initial states
        speciesNetwork.makeDummy();
        final NetworkNode origin = speciesNetwork.getOrigin();
        origin.setHeight(timeOrigin);
        final NetworkNode root = new NetworkNode(speciesNetwork);
        speciesNetwork.addSpeciationNode(root);
        origin.getChildren().add(root);
        root.getParents().add(origin);

        final List<NetworkNode> networkNodeList = new ArrayList<>();
        networkNodeList.add(root);
        boolean atRoot = true;

        double currentTime = timeOrigin;
        // start the simulation
        while (currentTime > 0.0) {
            final int k = networkNodeList.size();  // current number of branches
            final double totalRate = k*lambda + 0.5*k*(k-1)*nu;

            final double waitingTime;
            if (atRoot && originInput.get() == null) {
                waitingTime = 1.0 / lambda;
                atRoot = false;
            } else {
                waitingTime = Randomizer.nextExponential(totalRate);
            }
            currentTime -= waitingTime;

            if (currentTime > 0.0) {
                int rnd;
                if (Randomizer.nextDouble() <= k*lambda/totalRate) {
                    // speciation event, pick a random branch to split
                    rnd = Randomizer.nextInt(k);
                    final NetworkNode pNode = networkNodeList.get(rnd);
                    networkNodeList.remove(pNode);
                    pNode.setHeight(currentTime);

                    final NetworkNode cNode1 = new NetworkNode(speciesNetwork);
                    final NetworkNode cNode2 = new NetworkNode(speciesNetwork);
                    speciesNetwork.addSpeciationNode(cNode1);
                    speciesNetwork.addSpeciationNode(cNode2);
                    networkNodeList.add(cNode1);
                    networkNodeList.add(cNode2);

                    pNode.getChildren().add(cNode1);
                    pNode.getChildren().add(cNode2);
                    cNode1.getParents().add(pNode);
                    cNode2.getParents().add(pNode);
                } else {
                    // hybridization event, pick two branches to join
                    rnd = Randomizer.nextInt(k);
                    final NetworkNode pNode1 = networkNodeList.get(rnd);
                    networkNodeList.remove(pNode1);
                    rnd = Randomizer.nextInt(k-1);
                    final NetworkNode pNode2 = networkNodeList.get(rnd);
                    networkNodeList.remove(pNode2);
                    speciesNetwork.deleteNode(pNode1);
                    speciesNetwork.deleteNode(pNode2);
                    speciesNetwork.addReticulationNode(pNode1);  // use pNode1 as the hybrid node
                    pNode1.setHeight(currentTime);
                    pNode1.setGammaProb(Randomizer.nextDouble());

                    final NetworkNode cNode = new NetworkNode(speciesNetwork);
                    speciesNetwork.addSpeciationNode(cNode);
                    networkNodeList.add(cNode);

                    pNode1.getParents().addAll(pNode2.getParents());
                    for (NetworkNode parent: pNode2.getParents()) {
                        parent.getChildren().remove(pNode2);
                        parent.getChildren().add(pNode1);
                    }
                    pNode1.getChildren().add(cNode);
                    cNode.getParents().add(pNode1);
                }
            }
        }

        // reached the present, set node labels
        for (int i = 0; i < networkNodeList.size(); i++) {
            NetworkNode node = networkNodeList.get(i);
            speciesNetwork.deleteNode(node);
            speciesNetwork.addLeafNode(node);
            node.setLabel("T" + (i+1));
        }
        speciesNetwork.resetInternalNodeLabels();

        // deal with branch numbers and relationships
        for (int i = 0; i < speciesNetwork.getNodeCount(); i++) {
            NetworkNode node = speciesNetwork.getNode(i);
            node.setNr(i);
            node.gammaBranchNumber = speciesNetwork.getBranchNumber(i);
        }
        speciesNetwork.resetAllVisited();
        setChildBranchNrs(speciesNetwork.getOrigin());
    }

    private Integer setChildBranchNrs(NetworkNode node) {
        if (node.isVisited()) {
            return node.gammaBranchNumber + 1;
        } else {
            for (NetworkNode child: node.getChildren()) {
                node.childBranchNumbers.add(setChildBranchNrs(child));
            }
            node.setVisited(true);
            return node.gammaBranchNumber;
        }
    }
}
