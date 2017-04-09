package speciesnetwork.tools;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multiset;
import com.google.common.collect.Multisets;
import com.google.common.collect.Table;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Runnable;
import beast.evolution.tree.Tree;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;
import speciesnetwork.NetworkParser;

/**
 * @author Huw Ogilvie
 */

@Description("Summarize the posterior distribution of species networks.")
public class SummarizePosterior extends Runnable {
    public final Input<String> inputFileNameInput = new Input<>("inputFileName",
            "Name of the file that contains networks in extended newick format.", Validate.REQUIRED);
    public final Input<String> outputFileNameInput = new Input<>("outputFileName",
            "If provided, write to this file rather than to standard out.");
    public final Input<Integer> burninInput = new Input<>("burnin", "The absolute burn-in.", 0);

    private Map<Integer, Integer> rSubnetworks;
    private Table<Integer, Integer, Integer> sSubnetworks;

    private int nextSubnetworkNumber;

    @Override
    public void initAndValidate() {
    	rSubnetworks = new HashMap<>(); // subnetworks defined by a reticulation node
    	sSubnetworks = HashBasedTable.create(); // subnetworks defined by a speciation node
    }

    @Override
    public void run() throws IOException {
        final String inputFileName = inputFileNameInput.get();
        /* final String outputFileName = outputFileNameInput.get();

        PrintStream out;  // where to print
        if (outputFileName == null) {
            out = System.out;
        } else {
            String msg = "Writing";
            if (new File(outputFileName).exists())
                msg = "Warning: Overwriting";
            System.err.println(msg + " file " + outputFileName);
            out = new PrintStream(outputFileName);
        } */

        final int burnin = burninInput.get();
        NexusParser nexusParser = new NexusParser();
        nexusParser.parseFile(new File(inputFileName));
        List<Tree> parsedTrees = nexusParser.trees;

        nextSubnetworkNumber = nexusParser.taxa.size();
        final Multimap<Integer, Network> binnedNetworks = HashMultimap.create(); // binned by topology

        // bin networks by topology
        for (int i = 0; i < parsedTrees.size(); i++) {
        	if (i >= burnin) {
        		final Tree tree = parsedTrees.get(i);
        		final Network network = new NetworkParser(tree);
        		final NetworkNode root = network.getRoot();
        		final int networkNr = findSubnetworks(root);
    			binnedNetworks.put(networkNr, network);
        	}
        }

        // in descending order of frequency, calculate mean node heights by topology
        Set<Integer> topologyOrder = Multisets.copyHighestCountFirst(binnedNetworks.keys()).elementSet();
        for (Integer networkNr: topologyOrder) {
            final Multimap<Integer, Double> networkHeights = HashMultimap.create();
            final Table<Integer, Integer, List<Double>> networkGammas = HashBasedTable.create();
            for (Network network: binnedNetworks.get(networkNr)) {
            	NetworkNode root = network.getRoot();
            	collateParameters(root, null, null, networkHeights, networkGammas);
            }

            for (Network network: binnedNetworks.get(networkNr)) {
            	NetworkNode root = network.getRoot();
            	averageParameters(root, null, null, networkHeights, networkGammas);
            	System.out.println(network.toString());
            	break;
            }
        }
    }

    /* Figure out what subnetwork this node defines,
     * and whether it already has an assigned number.
     * If not, assign nextSubnetworkNumber.
     */
    private Integer findSubnetworks(NetworkNode node) {
    	Multiset<NetworkNode> children = node.getChildren();

    	Integer subnetworkNr = -1;
    	if (children.size() == 0) { // leaf node
    		subnetworkNr = node.getNr();
    	} else if (children.size() == 1) { // reticulation node
    		for (Integer branchNr: node.childBranchNumbers) { // only runs once
    			final NetworkNode child = node.getChildByBranch(branchNr);
    			final int childSubnetworkNr = findSubnetworks(child);
    			subnetworkNr = rSubnetworks.get(childSubnetworkNr);
	    		if (subnetworkNr == null) {
	    			subnetworkNr = nextSubnetworkNumber;
	    			rSubnetworks.put(childSubnetworkNr, subnetworkNr);
	    			nextSubnetworkNumber++;
	    		}
    		}
    	} else { // speciation node
    		int leftSubnetworkNr = -1;
    		int rightSubnetworkNr = -1;
    		for (Integer branchNr: node.childBranchNumbers) { // only runs twice
    			final NetworkNode child = node.getChildByBranch(branchNr);
    			final int childSubnetwork = findSubnetworks(child);

    			// "left" subnetwork is always the larger number
    			if (leftSubnetworkNr < childSubnetwork) {
    				rightSubnetworkNr = leftSubnetworkNr;
    				leftSubnetworkNr = childSubnetwork;
    			} else {
    				rightSubnetworkNr = childSubnetwork;
    			}
    		}

    		subnetworkNr = sSubnetworks.get(leftSubnetworkNr, rightSubnetworkNr);
    		if (subnetworkNr == null) {
    			subnetworkNr = nextSubnetworkNumber;
    			sSubnetworks.put(leftSubnetworkNr, rightSubnetworkNr, subnetworkNr);
    			nextSubnetworkNumber++;
    		}
    	}

    	node.subnetworkNr = subnetworkNr;
    	return subnetworkNr;
    }

    private void collateParameters(NetworkNode node, Integer parentSubnetworkNr, Integer parentBranchNr, Multimap<Integer, Double> heights, Table<Integer, Integer, List<Double>> gammas) {
    	Multiset<NetworkNode> children = node.getChildren();

    	final Integer subnetworkNr = node.subnetworkNr;
    	final Double nodeHeight = node.getHeight();
		heights.put(subnetworkNr, nodeHeight);

    	if (children.size() == 1) { // reticulation node
    		for (Integer branchNr: node.childBranchNumbers) { // only runs once
    			final NetworkNode child = node.getChildByBranch(branchNr);
    			collateParameters(child, subnetworkNr, branchNr, heights, gammas);

    			if (!gammas.contains(subnetworkNr, parentSubnetworkNr))
    				gammas.put(subnetworkNr, parentSubnetworkNr, new ArrayList<>());

    			final Double nodeGamma = node.getGammaProb();
    			if (parentBranchNr == node.gammaBranchNumber)
    				gammas.get(subnetworkNr, parentSubnetworkNr).add(nodeGamma);
    			else
    				gammas.get(subnetworkNr, parentSubnetworkNr).add(1.0 - nodeGamma);
    		}
    	} else if (children.size() == 2) { // speciation node
    		for (Integer branchNr: node.childBranchNumbers) { // only runs twice
    			final NetworkNode child = node.getChildByBranch(branchNr);
    			collateParameters(child, subnetworkNr, branchNr, heights, gammas);
    		}
    	}
    }

    /*
     * Sets the node heights and gammas to the mean average across all samples sharing the same network topology
     */
    private void averageParameters(NetworkNode node, Integer parentSubnetworkNr, Integer parentBranchNr, Multimap<Integer, Double> heights, Table<Integer, Integer, List<Double>> gammas) {
    	Multiset<NetworkNode> children = node.getChildren();

    	final Integer subnetworkNr = node.subnetworkNr;
    	final Double[] sampledHeights = (Double[]) heights.get(subnetworkNr).toArray();
    	final Double meanHeight = calculateMean(sampledHeights);
		node.setHeight(meanHeight);

    	if (children.size() == 1) { // reticulation node
    		for (Integer branchNr: node.childBranchNumbers) { // only runs once
    			final NetworkNode child = node.getChildByBranch(branchNr);
    			averageParameters(child, subnetworkNr, branchNr, heights, gammas);

    			if (parentBranchNr == node.gammaBranchNumber) {
    		    	final Double[] sampledGammas = (Double[]) gammas.get(subnetworkNr, parentSubnetworkNr).toArray();
    		    	final Double meanGamma = calculateMean(sampledGammas);
    				node.setGammaProb(meanGamma);
    			}
    		}
    	} else if (children.size() == 2) { // speciation node
    		for (Integer branchNr: node.childBranchNumbers) { // only runs twice
    			final NetworkNode child = node.getChildByBranch(branchNr);
    			averageParameters(child, subnetworkNr, branchNr, heights, gammas);
    		}
    	}
    }

    private double calculateMean(Double[] values) {
    	double sumValues = 0.0;
    	for (Double v: values) sumValues += v;

    	return sumValues / values.length;
    }
}
