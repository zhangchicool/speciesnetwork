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
public class SummarizeNetwork extends Runnable {
    public final Input<String> inputFileNameInput = new Input<>("inputFileName",
            "Name of the file that contains networks in extended newick format.", Validate.REQUIRED);
    public final Input<String> outputFileNameInput = new Input<>("outputFileName",
            "If provided, write to this file rather than to standard out.");

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
        out.println("nHybrid  length  height  tHybrid  tSpecia");

        int numNetworks = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(inputFileName))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.trim().toLowerCase().startsWith("tree ")) {
                    // process the line.
                    final int i = line.indexOf('(');
                    if (i > 0) line = line.substring(i);
                    TreeParser tree = new TreeParser(line);
                    NetworkParser network = new NetworkParser(tree);

                    printSummary(out, network);  // summarize

                    // final int nHybrid = network.getReticulationNodeCount();
                    // nHybridInNetworkMap.put(nHybrid, network);
                    numNetworks++;
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.err.println(numNetworks + " networks processed");
    }

    /**
     * print some summary statics, including
     * 1. number of hybridization
     * 2. network length
     * 3. root height
     * 4. time of youngest hybridization
     * 5. time of youngest speciation
     */
    private void printSummary(PrintStream out, Network network) {
        final int nHybrid = network.getReticulationNodeCount();
        final double length = network.getNetworkLength();
        final double height = network.getRoot().getHeight();

        double tHybrid = 0.0;
        for (NetworkNode hNode: network.getReticulationNodes()) {
            if (tHybrid == 0.0 || tHybrid > hNode.getHeight())
                tHybrid = hNode.getHeight();
        }
        double tSpecia = 0.0;
        for (NetworkNode sNode: network.getSpeciationNodes()) {
            if (tSpecia == 0.0 || tSpecia > sNode.getHeight())
                tSpecia = sNode.getHeight();
        }

        DecimalFormat df = new DecimalFormat("#.########");
        out.println(nHybrid + "\t" + df.format(length) + "\t" + df.format(height) + "\t" + df.format(tHybrid) + "\t" + df.format(tSpecia));
    }
}
