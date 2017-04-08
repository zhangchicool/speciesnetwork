package speciesnetwork.tools;

import java.io.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;
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
 * @author Huw A. Ogilvie
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
    	rSubnetworks = new HashMap<>();
    	sSubnetworks = HashBasedTable.create();
    }

    @Override
    public void run() throws IOException {
        final String inputFileName = inputFileNameInput.get();
        final String outputFileName = outputFileNameInput.get();

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
        out.println("nHybrid  length  height  tHybrid  gamma");

        final int burnin = burninInput.get();
        NexusParser nexusParser = new NexusParser();
        nexusParser.parseFile(new File(inputFileName));
        List<Tree> parsedTrees = nexusParser.trees;

        nextSubnetworkNumber = nexusParser.taxa.size();
        Multiset<Integer> networkCounter = HashMultiset.create();
        for (int i = 0; i < parsedTrees.size(); i++) {
        	if (i >= burnin) {
        		Tree tree = parsedTrees.get(i);
        		Network network = new NetworkParser(tree);
        		NetworkNode root = network.getRoot();
        		final int networkNr = findSubnetworks(root);
        		networkCounter.add(networkNr);
        		System.out.println(String.format("Sample %d = network %d", i, networkNr));
        	}
        }

        for (Integer networkNr: networkCounter.elementSet()) {
        	System.out.println(String.format("Network %d = %d counts", networkNr, networkCounter.count(networkNr)));
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
    		for (NetworkNode child: children) { // only runs once
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
    		for (NetworkNode child: children) { // only runs twice
    			final int childSubnetwork = findSubnetworks(child);
    			if (leftSubnetworkNr == -1)
    				leftSubnetworkNr = childSubnetwork;
    			else
    				rightSubnetworkNr = childSubnetwork;
    		}

    		subnetworkNr = sSubnetworks.get(leftSubnetworkNr, rightSubnetworkNr);
    		if (subnetworkNr == null) {
    			subnetworkNr = nextSubnetworkNumber;
    			sSubnetworks.put(leftSubnetworkNr, rightSubnetworkNr, subnetworkNr);
    			sSubnetworks.put(rightSubnetworkNr, leftSubnetworkNr, subnetworkNr);
    			nextSubnetworkNumber++;
    		} 
    	}

    	return subnetworkNr;
    }
}
