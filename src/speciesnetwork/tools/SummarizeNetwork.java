package speciesnetwork.tools;

import java.io.*;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Runnable;
import beast.util.TreeParser;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;
import speciesnetwork.NetworkParser;
import speciesnetwork.NodeHeightComparator;

/**
 * @author Chi Zhang
 */

@Description("Summarize the posterior distribution of species networks.")
@Deprecated
public class SummarizeNetwork extends Runnable {
    public final Input<String> inputFileNameInput = new Input<>("inputFileName",
            "Name of the file that contains networks in extended newick format.", Validate.REQUIRED);
    public final Input<String> outputFileNameInput = new Input<>("outputFileName",
            "If provided, write to this file rather than to standard out.");
    public final Input<Integer> burninInput = new Input<>("burnin", "The absolute burn-in.", 0);

    // map the number of reticulations with the networks
    private Multimap<Integer, Network> nHybridInNetworkMap = HashMultimap.create();

    private static Comparator<NetworkNode> hc = new NodeHeightComparator();

    @Override
    public void initAndValidate() {
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
        int numNetworks = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(inputFileName))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.trim().toLowerCase().startsWith("tree ")) {
                    // process the line (newick network string)
                    final int i = line.indexOf('(');
                    if (i > 0) line = line.substring(i);
                    TreeParser tree = new TreeParser(line);
                    NetworkParser network = new NetworkParser(tree);

                    numNetworks++;
                    if (numNetworks > burnin)
                        printSummary(out, network);  // summarize
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.err.println(numNetworks + " networks totally, " + burnin + " discarded as burn-in.");
    }

    private double getGammaProb(Network network) {
        if (network.getReticulationNodeCount() != 1)
            return Double.NaN;

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

        NetworkNode parentA = tipA.getParentByBranch(tipA.gammaBranchNumber);
        NetworkNode parentB = tipB.getParentByBranch(tipB.gammaBranchNumber);
        NetworkNode parentC = tipC.getParentByBranch(tipC.gammaBranchNumber);
        if (parentB.isReticulation() && parentA.getChildren().contains(parentB) && parentC.getChildren().contains(parentB)) {
            // System.out.println(parentB.getHeight() + "\t" + parentA.getHeight() + "\t" + parentC.getHeight()  + "\t" +
            //                    network.getRoot().getHeight() + "\t" + network.getOrigin().getHeight());
            if (parentB.getParentByBranch(parentB.gammaBranchNumber) == parentA)
                return parentB.getGammaProb();
            else
                return 1.0 - parentB.getGammaProb();
        }
        else return Double.NaN;
    }

    /**
     * print some summary statics, including
     * 1. number of hybridization
     * 2. network length
     * 3. root height
     * 4. time of youngest hybridization
     * 5. gamma (NA if not the true network)
     */
    private void printSummary(PrintStream out, Network network) {
        final int nHybrid = network.getReticulationNodeCount();
        final double length = network.getNetworkLength();
        final double height = network.getRoot().getHeight();
        final double gamma = getGammaProb(network);

        double tHybrid = 0.0;
        for (NetworkNode hNode: network.getReticulationNodes()) {
            if (tHybrid == 0.0 || tHybrid > hNode.getHeight())
                tHybrid = hNode.getHeight();
        }

        DecimalFormat df = new DecimalFormat("#.########");
        out.println(nHybrid + "\t" + df.format(length) + "\t" + df.format(height) + "\t" + df.format(tHybrid) + "\t" + gamma);
    }
}
