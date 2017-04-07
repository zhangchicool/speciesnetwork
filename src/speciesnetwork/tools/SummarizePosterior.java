package speciesnetwork.tools;

import java.io.*;
import java.util.List;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;

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

    // need to do this more efficiently - uses like 800MB of RAM
    final private int initialArraySize = 10000;

    private Integer[] rSubnetworks;
    private Integer[][] sSubnetworks;

    private int nextSubnetworkNumber;

    @Override
    public void initAndValidate() {
    	rSubnetworks = new Integer[initialArraySize];
    	sSubnetworks = new Integer[initialArraySize][];
    	for (int i = 0; i < sSubnetworks.length; i++)
    		sSubnetworks[i] = new Integer[initialArraySize];
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
    			subnetworkNr = rSubnetworks[childSubnetworkNr];
	    		if (subnetworkNr == null) {
	    			subnetworkNr = nextSubnetworkNumber;
	    			rSubnetworks[childSubnetworkNr] = subnetworkNr;
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

    		subnetworkNr = sSubnetworks[leftSubnetworkNr][rightSubnetworkNr];
    		if (subnetworkNr == null) {
    			subnetworkNr = nextSubnetworkNumber;
    			sSubnetworks[leftSubnetworkNr][rightSubnetworkNr] = subnetworkNr;
    			sSubnetworks[rightSubnetworkNr][leftSubnetworkNr] = subnetworkNr;
    			nextSubnetworkNumber++;
    		} 
    	}

    	return subnetworkNr;
    }
}
