package speciesnetwork.tools;

import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;
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
 * @author Chi Zhang
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
    public final Input<Integer> decimalPlacesInput = new Input<>("dp",
            "The number of decimal places to use (default -1 for full precision)", -1);

    private Map<Integer, Integer> rSubnetworks;
    private Table<Integer, Integer, Integer> sSubnetworks;
    private int nextSubnetworkNumber;

    private boolean useMedian;
    private DecimalFormat df;
    private static PrintStream progressStream = Log.info;

    @Override
    public void initAndValidate() {
    	rSubnetworks = new HashMap<>(); // subnetworks defined by a reticulation node
    	sSubnetworks = HashBasedTable.create(); // subnetworks defined by a speciation node
        useMedian = medianInput.get();

        int dp = decimalPlacesInput.get();
        if (dp < 0) {
            df = null;
        } else {
            // just new DecimalFormat("#.######") (with dp time '#' after the decimal)
            df = new DecimalFormat("#." + new String(new char[dp]).replace('\0', '#'));
            df.setRoundingMode(RoundingMode.HALF_UP);
        }
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
                final int networkNr = findSubnetworks(network.getOrigin());
                binnedNetworks.put(networkNr, network);
            }
        }
        progressStream.println("Parsed " + numNetworks + " networks totally, " + absoluteBurnin + " discarded as burn-in.");

        // in descending order of topology frequencies, calculate node summaries
        final Multiset<Integer> allNetworkNrs = binnedNetworks.keys();
        final ImmutableMultiset<Integer> orderedNetworkNrs = Multisets.copyHighestCountFirst(allNetworkNrs);
        final Set<Integer> uniqueNetworkNrs = new LinkedHashSet<>();
        uniqueNetworkNrs.addAll(orderedNetworkNrs);

        progressStream.println("Writing summary networks with heights and gamma probs...");
        for (Integer networkNr: uniqueNetworkNrs) {
            final ListMultimap<Integer, Double> networkHeights = ArrayListMultimap.create();
            final Table<Integer, Integer, List<Double>> networkGammas = HashBasedTable.create();

            for (Network network: binnedNetworks.get(networkNr)) {
                NetworkNode origin = network.getOrigin();
                network.resetAllVisited();
                collateParameters(origin, null, null, networkHeights, networkGammas);
            }

            for (Network network: binnedNetworks.get(networkNr)) {
                NetworkNode origin = network.getOrigin();
                origin.topologySupport = (double) allNetworkNrs.count(networkNr) / (double) allNetworkNrs.size();
                network.resetAllVisited();
                summarizeParameters(origin, null, null, networkHeights, networkGammas);

                out.println(network.toString(df));
                break;
            }
        }

        out.close();
    }

    /*
     * Figure out what subnetwork this node defines. If it doesn't have an assigned number, assign nextSubnetworkNumber.
     */
    private Integer findSubnetworks(NetworkNode node) {
    	Integer subnetworkNr = -1;
    	if (node.isSpeciation() || node.isOrigin()) {  // speciation or origin node
    		int leftSubnetworkNr = -1;
    		int rightSubnetworkNr = -1;
    		for (Integer branchNr: node.childBranchNumbers) {
    			final NetworkNode child = node.getChildByBranch(branchNr);
    			final int childSubnetworkNr = findSubnetworks(child);

    			// "left" subnetwork is always the larger number
    			if (leftSubnetworkNr < childSubnetworkNr) {
    				rightSubnetworkNr = leftSubnetworkNr;
    				leftSubnetworkNr = childSubnetworkNr;
    			} else {
    				rightSubnetworkNr = childSubnetworkNr;
    			}
    		}

    		subnetworkNr = sSubnetworks.get(leftSubnetworkNr, rightSubnetworkNr);
    		if (subnetworkNr == null) {
    			subnetworkNr = nextSubnetworkNumber;
    			sSubnetworks.put(leftSubnetworkNr, rightSubnetworkNr, subnetworkNr);
    			nextSubnetworkNumber++;
    		}
    	} else if (node.isReticulation()) {  // reticulation node
    		Integer branchNr = node.childBranchNumbers.get(0); // one child
            final NetworkNode child = node.getChildByBranch(branchNr);
            final int childSubnetworkNr = findSubnetworks(child);

            subnetworkNr = rSubnetworks.get(childSubnetworkNr);
            if (subnetworkNr == null) {
                subnetworkNr = nextSubnetworkNumber;
                rSubnetworks.put(childSubnetworkNr, subnetworkNr);
                nextSubnetworkNumber++;
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
    private void collateParameters(NetworkNode node, Integer parentSubnetworkNr, Integer parentBranchNr,
                                   ListMultimap<Integer, Double> heights, Table<Integer, Integer, List<Double>> gammas) {
        final Integer subnetworkNr = node.subnetworkNr;

        if (node.isReticulation()) {
            if (!gammas.contains(subnetworkNr, parentSubnetworkNr))
                gammas.put(subnetworkNr, parentSubnetworkNr, new ArrayList<>());
            final double nodeGamma = node.getGammaProb();
            if (node.gammaBranchNumber.equals(parentBranchNr))
                gammas.get(subnetworkNr, parentSubnetworkNr).add(nodeGamma);
            else
                gammas.get(subnetworkNr, parentSubnetworkNr).add(1.0 - nodeGamma);
        }

        if (node.isVisited())
            return;
        // mark visited to avoid duplicated recursion
        node.setVisited(true);

        final Double nodeHeight = node.getHeight();
        heights.put(subnetworkNr, nodeHeight);

        for (Integer branchNr : node.childBranchNumbers) {
            final NetworkNode child = node.getChildByBranch(branchNr);
            collateParameters(child, subnetworkNr, branchNr, heights, gammas);
        }
    }

    /*
     * Summarize the node heights and gammas across all samples sharing the same network topology
     */
    private void summarizeParameters(NetworkNode node, Integer parentSubnetworkNr, Integer parentBranchNr,
                                     ListMultimap<Integer, Double> heights, Table<Integer, Integer, List<Double>> gammas) {
        final Integer subnetworkNr = node.subnetworkNr;

        if (node.isReticulation() && node.gammaBranchNumber.equals(parentBranchNr)) {
            final List<Double> sampledGammas = gammas.get(subnetworkNr, parentSubnetworkNr);
            final double meanGamma = calculateMean(sampledGammas);
            final double medianGamma = calculateMedian(sampledGammas);
            final double[] hpdGamma = calculateHPDInterval(0.95, sampledGammas);
            final double[] rangeGamma = calculateRange(sampledGammas);
            node.setMetaData("gamma_mean", meanGamma);
            node.setMetaData("gamma_median", medianGamma);
            node.setMetaData("gamma_95%HPD", new Object[]{hpdGamma[0], hpdGamma[1]});
            node.setMetaData("gamma_range", new Object[]{rangeGamma[0], rangeGamma[1]});

            if (useMedian) {
                node.setGammaProb(medianGamma);
            } else {
                node.setGammaProb(meanGamma);
            }
        }

        if (node.isVisited())
            return;
        // mark visited to avoid duplicated recursion
        node.setVisited(true);

        final List<Double> sampledHeights = heights.get(subnetworkNr);
        final double meanHeight = calculateMean(sampledHeights);
        final double medianHeight = calculateMedian(sampledHeights);
        final double[] hpdHeight = calculateHPDInterval(0.95, sampledHeights);
        final double[] rangeHeight = calculateRange(sampledHeights);
        node.setMetaData("height_mean", meanHeight);
        node.setMetaData("height_median", medianHeight);
        node.setMetaData("height_95%HPD", new Object[]{hpdHeight[0],hpdHeight[1]});
        node.setMetaData("height_range", new Object[]{rangeHeight[0], rangeHeight[1]});

        if (useMedian) {
            node.setHeight(medianHeight);
        } else {
            node.setHeight(meanHeight);
        }

        for (Integer branchNr: node.childBranchNumbers) {
            final NetworkNode child = node.getChildByBranch(branchNr);
            summarizeParameters(child, subnetworkNr, branchNr, heights, gammas);
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

    private double[] calculateRange(List<Double> sample) {
        double min = sample.get(0);
        double max = min;
        for (Double v : sample) {
            if (v > max) max = v;
            if (v < min) min = v;
        }
        return new double[]{min, max};
    }

    // 1st and 3rd quartiles
    private double[] calculateQuartiles(List<Double> sample) {
        final int length = sample.size();
        int[] indices = new int[length];
        HeapSort.sort(sample, indices);

        int mid = length / 2;
        int pos = mid / 2;
        double first, third;
        if (mid % 2 == 1 || length == 1) {
            first = sample.get(indices[pos]);
            third = sample.get(indices[length - pos - 1]);
        } else {
            first = (sample.get(indices[pos - 1]) + sample.get(indices[pos])) / 2.0;
            third = (sample.get(indices[length - pos - 1]) + sample.get(indices[length - pos])) / 2.0;
        }

        return new double[]{first, third};
    }
}
