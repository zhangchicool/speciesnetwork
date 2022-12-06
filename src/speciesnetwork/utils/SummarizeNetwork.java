package speciesnetwork.utils;

import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.*;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.util.Log;
import beast.base.core.Runnable;
import beast.base.util.HeapSort;
import beast.base.util.TreeParser;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;
import speciesnetwork.NetworkParser;

@Deprecated
public class SummarizeNetwork extends Runnable {
    public final Input<String> inputFileNameInput = new Input<>("inputFileName",
            "Name of the file that contains networks in extended newick format.", Validate.REQUIRED);
    public final Input<String> outputFileNameInput = new Input<>("outputFileName",
            "If provided, write to this file rather than to standard out.");
    public final Input<Integer> burninInput = new Input<>("burnin", "The absolute burn-in.", 0);
    public final Input<Boolean> medianInput = new Input<>("useMedian",
            "Use median instead of mean for node heights and gamma probs in summary networks.", false);
    public final Input<Integer> decimalPlacesInput = new Input<>("dp",
            "The number of decimal places to use (default -1 for full precision)", -1);

    private boolean useMedian;
    private DecimalFormat df;
    private static PrintStream progressStream = Log.info;

    @Override
    public void initAndValidate() {
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

        // get the absolute burn-in
        final int burnin = burninInput.get();

        progressStream.print("Parsing network samples ");
        List<Network> networks = new ArrayList<>();
        int numNetworks = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(inputFileName))) {
            String line;
            while ((line = br.readLine()) != null) {
                // process the line (newick network string)
                //if (line.trim().toLowerCase().startsWith("tree ")) {
                final int i = line.indexOf('(');
                if (i >= 0) {
                    line = line.substring(i);
                    TreeParser tree = new TreeParser(line);
                    NetworkParser network = new NetworkParser(tree);
                    numNetworks++;

                    if (numNetworks > burnin) {
                        networks.add(network);
                    }
                    if (numNetworks % 10000 == 0) progressStream.print(".");
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        progressStream.println("\nParsed " + numNetworks + " networks totally, " + burnin + " discarded as burn-in.");

        printSummary_net(networks, out);

        out.close();
    }

    private void printSummary_net (List<Network> networks, PrintStream out){
        // print header
        out.println("length  height  nHybrid  tHybrid");

        for (Network network : networks) {
            double rootHeight = network.getRoot().getHeight();
            double originHeight = network.getOrigin().getHeight();
            double length = network.getNetworkLength() - (originHeight - rootHeight);

            int nHybrid = network.getReticulationNodeCount();

            double tHybrid = 0.0;  // youngest hybridization time
            for (NetworkNode hNode : network.getReticulationNodes()) {
                if (tHybrid == 0.0 || tHybrid > hNode.getHeight())
                    tHybrid = hNode.getHeight();
            }

            out.println(length + "\t" + rootHeight + "\t" + nHybrid + "\t" + tHybrid);
        }
    }

    private void printSummary_h1 (List<Network> networks, PrintStream out){
        List<Double> gamma1 = new ArrayList<>();
        List<Double> hgt_H1 = new ArrayList<>();
        List<Double> height = new ArrayList<>();
        List<Double> length = new ArrayList<>();

        int nTrees = 0;
        int nTruth = 0;
        for (Network network : networks) {
            // network height and length
            double rootHeight = network.getRoot().getHeight();
            double originHeight = network.getOrigin().getHeight();
            double netLength = network.getNetworkLength() - (originHeight - rootHeight);
            height.add(rootHeight);
            length.add(netLength);

            NetworkNode tipA = null, tipB = null, tipC = null;
            for (NetworkNode tip: network.getLeafNodes()) {
                switch (tip.getLabel()) {
                    case "A":  tipA = tip;  break;
                    case "B":  tipB = tip;  break;
                    case "C":  tipC = tip;  break;
                    default:
                }
            }
            if (tipA == null || tipB == null || tipC == null)
                throw new RuntimeException("Check tip label!");

            // tree
            if (network.getReticulationNodeCount() == 0)
                nTrees++;

            // true network
            NetworkNode parentA = tipA.getParentByBranch(tipA.gammaBranchNumber);
            NetworkNode parentB = tipB.getParentByBranch(tipB.gammaBranchNumber);
            NetworkNode parentC = tipC.getParentByBranch(tipC.gammaBranchNumber);
            if (network.getReticulationNodeCount() == 1 && parentB.isReticulation() &&
                    parentA.getChildren().contains(parentB) && parentC.getChildren().contains(parentB)) {
                nTruth++;

                if (parentB.getParentByBranch(parentB.gammaBranchNumber) == parentA)
                    gamma1.add(parentB.getGammaProb());
                else
                    gamma1.add(1.0 - parentB.getGammaProb());
                hgt_H1.add(parentB.getHeight());
            }
        }

        // percentage of trees
        final int nNetworks = networks.size();
        out.print((double)nTrees/nNetworks + "\t");
        // percentage of true networks
        out.print((double)nTruth/nNetworks + "\t");

        // network height
        double median = calculateMedian(height);
        double[] hpd = calculateHPDInterval(0.95, height);
        out.print(median + "\t" + hpd[0] + "\t" + hpd[1] + "\t");

        // network length
        median = calculateMedian(length);
        hpd = calculateHPDInterval(0.95, length);
        out.print(median + "\t" + hpd[0] + "\t" + hpd[1] + "\t");

        // gamma1 if true network
        if (gamma1.size() == 0)
            out.print(Double.NaN + "\t" + Double.NaN + "\t" + Double.NaN + "\t");
        else {
            median = calculateMedian(gamma1);
            hpd = calculateHPDInterval(0.95, gamma1);
            out.print(median + "\t" + hpd[0] + "\t" + hpd[1] + "\t");
        }

        // h1 height if true network
        if (hgt_H1.size() == 0)
            out.print(Double.NaN + "\t" + Double.NaN + "\t" + Double.NaN + "\n");
        else {
            median = calculateMedian(hgt_H1);
            hpd = calculateHPDInterval(0.95, hgt_H1);
            out.print(median + "\t" + hpd[0] + "\t" + hpd[1] + "\n");
        }
    }

    private void printSummary_h2 (List<Network> networks, PrintStream out){
        List<Double> gamma1 = new ArrayList<>();
        List<Double> hgt_H1 = new ArrayList<>();
        List<Double> gamma2 = new ArrayList<>();
        List<Double> hgt_H2 = new ArrayList<>();
        List<Double> height = new ArrayList<>();
        List<Double> length = new ArrayList<>();

        int nTrees = 0;
        int nTruth = 0;  // number of true networks
        int nOneHy = 0;  // number of networks with 1 hybridization
        int nTwoHy = 0;  // number of networks with 1 hybridization
        int nBCDHy = 0;  // number of networks with the BCDH structure
        for (Network network : networks) {
            // network height and length
            double rootHeight = network.getRoot().getHeight();
            double originHeight = network.getOrigin().getHeight();
            double netLength = network.getNetworkLength() - (originHeight - rootHeight);
            height.add(rootHeight);
            length.add(netLength);

            NetworkNode tipA = null, tipB = null, tipC = null, tipD = null;
            for (NetworkNode tip: network.getLeafNodes()) {
                switch (tip.getLabel()) {
                    case "A":  tipA = tip;  break;
                    case "B":  tipB = tip;  break;
                    case "C":  tipC = tip;  break;
                    case "D":  tipD = tip;  break;
                    default:
                }
            }
            if (tipA == null || tipB == null || tipC == null || tipD == null)
                throw new RuntimeException("Check tip label!");

            // tree
            if (network.getReticulationNodeCount() == 0)
                nTrees++;
            // network with 1 hybridization
            if (network.getReticulationNodeCount() == 1)
                nOneHy++;
            // network with 2 hybridizations
            if (network.getReticulationNodeCount() == 2)
                nTwoHy++;

            NetworkNode parentA = tipA.getParentByBranch(tipA.gammaBranchNumber);  // S1
            NetworkNode parentB = tipB.getParentByBranch(tipB.gammaBranchNumber);  // S3
            NetworkNode parentC = tipC.getParentByBranch(tipC.gammaBranchNumber);  // H2
            NetworkNode parentD = tipD.getParentByBranch(tipD.gammaBranchNumber);  // S4
            NetworkNode parentParentB = parentB.getParentByBranch(parentB.gammaBranchNumber);  // H1
            NetworkNode parentParentD = parentD.getParentByBranch(parentD.gammaBranchNumber);  // S2

            // networks with the BCDH structure
            if (parentC.isReticulation() && parentB.getChildren().contains(parentC) && parentD.getChildren().contains(parentC)) {
                nBCDHy++;

                if (network.getReticulationNodeCount() == 1) {
                    // networks with 1 hybridization
                    if (parentC.getParentByBranch(parentC.gammaBranchNumber) == parentB)
                        gamma2.add(parentC.getGammaProb());
                    else
                        gamma2.add(1.0 - parentC.getGammaProb());
                    hgt_H2.add(parentC.getHeight());
                }
                else if (network.getReticulationNodeCount() == 2 && parentParentB.isReticulation() &&
                        parentA.getChildren().contains(parentParentB) && parentParentD.getChildren().contains(parentParentB)) {
                    // true networks
                    nTruth++;

                    if (parentParentB.getParentByBranch(parentParentB.gammaBranchNumber) == parentA)
                        gamma1.add(parentParentB.getGammaProb());
                    else
                        gamma1.add(1.0 - parentParentB.getGammaProb());
                    hgt_H1.add(parentParentB.getHeight());

                    if (parentC.getParentByBranch(parentC.gammaBranchNumber) == parentB)
                        gamma2.add(parentC.getGammaProb());
                    else
                        gamma2.add(1.0 - parentC.getGammaProb());
                    hgt_H2.add(parentC.getHeight());
                }
            }
        }

        final int nNetworks = networks.size();
        // percentage of trees
        out.print((double)nTrees/nNetworks + "\t");
        // percentage of 1 hybridization
        out.print((double)nOneHy/nNetworks + "\t");
        // percentage of 2 hybridizations
        out.print((double)nTwoHy/nNetworks + "\t");
        // percentage of true networks
        out.print((double)nTruth/nNetworks + "\t");
        // percentage of BCDH structure
        out.print((double)nBCDHy/nNetworks + "\t");

        // network height
        double median = calculateMedian(height);
        double[] hpd = calculateHPDInterval(0.95, height);
        out.print(median + "\t" + hpd[0] + "\t" + hpd[1] + "\t");

        // network length
        median = calculateMedian(length);
        hpd = calculateHPDInterval(0.95, length);
        out.print(median + "\t" + hpd[0] + "\t" + hpd[1] + "\t");

        // gamma1 if true network
        if (gamma1.size() == 0)
            out.print(Double.NaN + "\t" + Double.NaN + "\t" + Double.NaN + "\t");
        else {
            median = calculateMedian(gamma1);
            hpd = calculateHPDInterval(0.95, gamma1);
            out.print(median + "\t" + hpd[0] + "\t" + hpd[1] + "\t");
        }

        // h1 height if true network
        if (hgt_H1.size() == 0)
            out.print(Double.NaN + "\t" + Double.NaN + "\t" + Double.NaN + "\t");
        else {
            median = calculateMedian(hgt_H1);
            hpd = calculateHPDInterval(0.95, hgt_H1);
            out.print(median + "\t" + hpd[0] + "\t" + hpd[1] + "\t");
        }

        // gamma2 if BCDH structure
        if (gamma2.size() == 0)
            out.print(Double.NaN + "\t" + Double.NaN + "\t" + Double.NaN + "\t");
        else {
            median = calculateMedian(gamma2);
            hpd = calculateHPDInterval(0.95, gamma2);
            out.print(median + "\t" + hpd[0] + "\t" + hpd[1] + "\t");
        }

        // h2 height if BCDH structure
        if (hgt_H2.size() == 0)
            out.print(Double.NaN + "\t" + Double.NaN + "\t" + Double.NaN + "\n");
        else {
            median = calculateMedian(hgt_H2);
            hpd = calculateHPDInterval(0.95, hgt_H2);
            out.print(median + "\t" + hpd[0] + "\t" + hpd[1] + "\n");
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
