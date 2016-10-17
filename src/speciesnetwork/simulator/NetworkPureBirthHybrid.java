package speciesnetwork.simulator;

import java.io.*;
import java.util.*;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Runnable;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.util.Randomizer;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;

/**
 * @author Chi Zhang
 */

@Description("Simulate a species network under the pure birth and hybridization process.")
public class NetworkPureBirthHybrid extends Runnable {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "Species network to be simulated.", Validate.REQUIRED);
    public final Input<RealParameter> originInput =
            new Input<>("origin", "The time when the process started.", Validate.REQUIRED);
    public final Input<RealParameter> birthRateInput =
            new Input<>("birthRate", "Speciation rate, lambda.", Validate.REQUIRED);
    public final Input<RealParameter> hybridRateInput =
            new Input<>("hybridRate", "Hybridization rate, nu.", Validate.REQUIRED);

    // number of extant species to condition on (-1 for no such condition)
    private int numTips = -1;

    @Override
    public void initAndValidate() {
        Network speciesNetwork = speciesNetworkInput.get();
        TaxonSet species = speciesNetwork.taxonSetInput.get();
        if (species != null)
            numTips = species.asStringList().size();
    }

    @Override
    public void run() {
        simulate();
    }

    public Network simulate() {
        Network speciesNetwork = speciesNetworkInput.get();
        NetworkNode[] networkTips = speciesNetwork.getLeafNodes();
        final double timeOrigin = originInput.get().getValue();
        final double lambda = birthRateInput.get().getValue();
        final double nu = hybridRateInput.get().getValue();

        // set the initial states
        speciesNetwork.makeDummy();

        final List<NetworkNode> networkNodeList = new ArrayList<>();
        networkNodeList.add(speciesNetwork.getRoot());
        double currentTime = timeOrigin;
        // start the simulation
        while (currentTime > 0.0) {
            final int k = networkNodeList.size();  // current number of branches
            final double totalRate = k*lambda + 0.5*k*(k-1)*nu;
            final double waitingTime = Randomizer.nextExponential(totalRate);
            currentTime -= waitingTime;

            if (currentTime > 0.0) {
                int rnd;
                if (Randomizer.nextDouble() <= k*lambda/totalRate) {
                    // speciation event, pick a random branch to split
                    rnd = Randomizer.nextInt(k);
                    final NetworkNode pNode = networkNodeList.get(rnd);
                    networkNodeList.remove(pNode);

                    final NetworkNode cNode1 = new NetworkNode(speciesNetwork);
                    final NetworkNode cNode2 = new NetworkNode(speciesNetwork);
                    speciesNetwork.addNode(cNode1);
                    speciesNetwork.addNode(cNode2);
                    networkNodeList.add(cNode1);
                    networkNodeList.add(cNode2);

                    pNode.setHeight(currentTime);
                    pNode.getChildren().add(cNode1);
                    pNode.getChildren().add(cNode2);
                    cNode1.getParents().add(pNode);
                    cNode2.getParents().add(pNode);
                } else {
                    // hybridization event, pick two branches to join
                    rnd = Randomizer.nextInt(k);
                    final NetworkNode pNode1 = networkNodeList.get(rnd);
                    rnd = Randomizer.nextInt(k-1);
                    final NetworkNode pNode2 = networkNodeList.get(rnd);
                    networkNodeList.remove(pNode1);
                    networkNodeList.remove(pNode2);
                    final NetworkNode hNode = speciesNetwork.joinNode(pNode1, pNode2);

                    final NetworkNode cNode = new NetworkNode(speciesNetwork);
                    speciesNetwork.addNode(cNode);
                    networkNodeList.add(cNode);

                    hNode.setHeight(currentTime);
                    hNode.getChildren().add(cNode);
                    cNode.getParents().add(hNode);
                }
            }
        }

        // reached the present, set tip heights to zero
        for (NetworkNode node: networkNodeList) {
            node.setHeight(0.0);
        }

        // TODO: deal with branch numbers and relationships


        return speciesNetwork;
    }




}
