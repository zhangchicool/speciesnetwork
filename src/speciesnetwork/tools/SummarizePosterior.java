package speciesnetwork.tools;

import java.io.*;
import java.util.*;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multiset;
import com.google.common.collect.Multisets;
import com.google.common.collect.Table;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.util.Log;
import beast.core.Runnable;
import beast.evolution.tree.Tree;
import beast.util.HeapSort;
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
    public final Input<Double> burninInput = new Input<>("burnin",
            "Absolute burn-in if >=1, or relative burn-in if <1.", 0.);
    public final Input<Boolean> medianInput = new Input<>("useMedian",
            "Use median instead of mean for node heights and gamma probs in summary networks.", false);

    private Map<Integer, Integer> rSubnetworks;
    private Table<Integer, Integer, Integer> sSubnetworks;

    private int nextSubnetworkNumber;

    private static PrintStream progressStream = Log.info;

    @Override
    public void initAndValidate() {
    	rSubnetworks = new HashMap<>(); // subnetworks defined by a reticulation node
    	sSubnetworks = HashBasedTable.create(); // subnetworks defined by a speciation node
    }

    @Override
    public void run() throws IOException {
        final String inputFileName = inputFileNameInput.get();
        final String outputFileName = outputFileNameInput.get();

        PrintStream out;  // where to print
        if (outputFileName == null) {
            out = System.out;
        } else {
            String msg = "Writing to";
            if (new File(outputFileName).exists())
                msg = "Warning: Overwriting";
            progressStream.println(msg + " file " + outputFileName);
            out = new PrintStream(outputFileName);
        }

        progressStream.println("Reading in posterior network samples...");
        ExtNexusParser nexusParser = new ExtNexusParser();
        nexusParser.parseFile(new File(inputFileName));
        List<Tree> parsedTrees = nexusParser.trees;
        final int numNetworks = parsedTrees.size();

        // get the absolute burn-in
        final double burnin = burninInput.get();
        int absoluteBurnin = 0;
        if (burnin > 0.0 && burnin < 1.0)
            absoluteBurnin = (int)(burnin * numNetworks);
        else if (burnin >= 1.0)
            absoluteBurnin = (int) burnin;

        nextSubnetworkNumber = nexusParser.taxa.size();
        final Multimap<Integer, Network> binnedNetworks = HashMultimap.create(); // binned by topology

        // bin networks by topology
        for (int i = 0; i < numNetworks; i++) {
        	if (i >= absoluteBurnin) {
        		final Tree tree = parsedTrees.get(i);
        		final Network network = new NetworkParser(tree);
        		final NetworkNode origin = network.getOrigin();
        		final int networkNr = findSubnetworks(origin, true);
    			binnedNetworks.put(networkNr, network);
        	}
        }
        progressStream.println("Parsed " + numNetworks + " networks totally, " + absoluteBurnin + " discarded as burn-in.");

        // in descending order of frequency, calculate mean node heights by topology
        final Multiset<Integer> allNetworkNrs = binnedNetworks.keys();
        final ImmutableMultiset<Integer> orderedNetworkNrs = Multisets.copyHighestCountFirst(allNetworkNrs);
        final Set<Integer> uniqueNetworkNrs = new LinkedHashSet<>();
        // for (Integer networkNr: orderedNetworkNrs) uniqueNetworkNrs.add(networkNr);
		uniqueNetworkNrs.addAll(orderedNetworkNrs);

        progressStream.println("Writing summary networks with mean heights and gammas...");
        for (Integer networkNr: uniqueNetworkNrs) {
        	final double topologySupport = (double) allNetworkNrs.count(networkNr) / (double) allNetworkNrs.size();
            final ListMultimap<Integer, Double> networkHeights = ArrayListMultimap.create();
            final Table<Integer, Integer, List<Double>> networkGammas = HashBasedTable.create();

            for (Network network: binnedNetworks.get(networkNr)) {
            	NetworkNode origin = network.getOrigin();
            	collateParameters(origin, null, null, networkHeights, networkGammas);
            }

            for (Network network: binnedNetworks.get(networkNr)) {  // does not loop ??
            	NetworkNode origin = network.getOrigin();
            	origin.topologySupport = topologySupport;
            	averageParameters(origin, null, null, networkHeights, networkGammas);
            	out.println(String.format("%s;", network.toString()));
            	break;
            }
        }

        out.close();
    }

    /* Figure out what subnetwork this node defines,
     * and whether it already has an assigned number.
     * If not, assign nextSubnetworkNumber.
     */
    private Integer findSubnetworks(NetworkNode node, boolean isOrigin) {
    	Multiset<NetworkNode> children = node.getChildren();

    	Integer subnetworkNr = -1;
    	if (children.size() == 2 || isOrigin) { // speciation or origin node
    		int leftSubnetworkNr = -1;
    		int rightSubnetworkNr = -1;
    		for (Integer branchNr: node.childBranchNumbers) { // only runs twice
    			final NetworkNode child = node.getChildByBranch(branchNr);
    			final int childSubnetwork = findSubnetworks(child, false);

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
    	} else if (children.size() == 1) { // reticulation node
    		for (Integer branchNr: node.childBranchNumbers) { // only runs once
    			final NetworkNode child = node.getChildByBranch(branchNr);
    			final int childSubnetworkNr = findSubnetworks(child, false);
    			subnetworkNr = rSubnetworks.get(childSubnetworkNr);
	    		if (subnetworkNr == null) {
	    			subnetworkNr = nextSubnetworkNumber;
	    			rSubnetworks.put(childSubnetworkNr, subnetworkNr);
	    			nextSubnetworkNumber++;
	    		}
    		}
    	} else { // is a leaf node
    		subnetworkNr = node.getNr();
    	}

    	node.subnetworkNr = subnetworkNr;
    	return subnetworkNr;
    }

    /*
     * Collate all node heights and reticulation node gammas, for networks sharing a common topology
     */
    private void collateParameters(NetworkNode node, Integer parentSubnetworkNr, Integer parentBranchNr, Multimap<Integer, Double> heights, Table<Integer, Integer, List<Double>> gammas) {
    	Multiset<NetworkNode> children = node.getChildren();

    	final Integer subnetworkNr = node.subnetworkNr;
    	final Double nodeHeight = node.getHeight();
		heights.put(subnetworkNr, nodeHeight);

		if (children.size() == 2 || parentSubnetworkNr == null) { // speciation node or origin
    		for (Integer branchNr: node.childBranchNumbers) { // only runs twice
    			final NetworkNode child = node.getChildByBranch(branchNr);
    			collateParameters(child, subnetworkNr, branchNr, heights, gammas);
    		}
    	} else if (children.size() == 1) { // reticulation node
    		for (Integer branchNr: node.childBranchNumbers) { // only runs once
    			final NetworkNode child = node.getChildByBranch(branchNr);
    			collateParameters(child, subnetworkNr, branchNr, heights, gammas);

    			if (!gammas.contains(subnetworkNr, parentSubnetworkNr))
    				gammas.put(subnetworkNr, parentSubnetworkNr, new ArrayList<>());

    			final Double nodeGamma = node.getGammaProb();
    			if (parentBranchNr.equals(node.gammaBranchNumber))
    				gammas.get(subnetworkNr, parentSubnetworkNr).add(nodeGamma);
    			else
    				gammas.get(subnetworkNr, parentSubnetworkNr).add(1.0 - nodeGamma);
    		}
    	}
    }

    /*
     * Sets the node heights and gammas to the mean average across all samples sharing the same network topology
     */
    private void averageParameters(NetworkNode node, Integer parentSubnetworkNr, Integer parentBranchNr, ListMultimap<Integer, Double> heights, Table<Integer, Integer, List<Double>> gammas) {
    	Multiset<NetworkNode> children = node.getChildren();

    	final Integer subnetworkNr = node.subnetworkNr;
    	final List<Double> sampledHeights = heights.get(subnetworkNr);
    	final Double meanHeight = calculateMean(sampledHeights);
		node.setHeight(meanHeight);

		if (children.size() == 2 || parentSubnetworkNr == null) { // speciation or origin node
    		for (Integer branchNr: node.childBranchNumbers) { // only runs twice
    			final NetworkNode child = node.getChildByBranch(branchNr);
    			averageParameters(child, subnetworkNr, branchNr, heights, gammas);
    		}
    	} else if (children.size() == 1) { // reticulation node
    		for (Integer branchNr: node.childBranchNumbers) { // only runs once
    			final NetworkNode child = node.getChildByBranch(branchNr);
    			averageParameters(child, subnetworkNr, branchNr, heights, gammas);

    			if (parentBranchNr.equals(node.gammaBranchNumber)) {
    		    	final List<Double> sampledGammas = gammas.get(subnetworkNr, parentSubnetworkNr);
    		    	final Double meanGamma = calculateMean(sampledGammas);
    				node.setGammaProb(meanGamma);
    			}
    		}
    	}
    }

    private double calculateMean(List<Double> sample) {
        double avg = 0.0;
    	int n = 1;
    	for (Double v : sample){
    		avg += (v - avg) / n;
    		n++;
		}  // more numerically stable
		return avg;
    }

    private double calculateMedian(List<Double> sample) {
        final int length = sample.size();
        int[] indices = new int[length];
        HeapSort.sort(sample, indices);

        int pos = length / 2;
        if (length % 2 == 1) {
            return sample.get(indices[pos]);
        } else {
            return (sample.get(indices[pos - 1]) + sample.get(indices[pos])) / 2.0;
        }
    }

    private double[] calculateHPDInterval(double prop, List<Double> sample) {
        final int length = sample.size();
        int[] indices = new int[length];
        HeapSort.sort(sample, indices);

        double minRange = Double.MAX_VALUE;
        int hpdIndex = 0;

        int diff = (int) Math.round(prop * length);
        for (int i = 0; i <= (length - diff); i++) {
            double minValue = sample.get(indices[i]);
            double maxValue = sample.get(indices[i + diff - 1]);
            double range = Math.abs(maxValue - minValue);
            if (range < minRange) {
                minRange = range;
                hpdIndex = i;
            }
        }
        double lower = sample.get(indices[hpdIndex]);
        double upper = sample.get(indices[hpdIndex + diff - 1]);

        return new double[]{lower, upper};
    }
}
