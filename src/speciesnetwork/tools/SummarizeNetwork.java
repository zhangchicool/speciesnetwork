package speciesnetwork.tools;

import java.io.*;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Runnable;
import beast.util.TreeParser;
import speciesnetwork.Network;
import speciesnetwork.NetworkParser;

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

    @Override
    public void initAndValidate() {

    }

    @Override
    public void run() throws IOException {
        final String inFileName = inputFileNameInput.get();
        final String outFileName = outputFileNameInput.get();

        int numNetworks = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(inFileName))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.trim().toLowerCase().startsWith("tree ")) {
                    // process the line.
                    final int i = line.indexOf('(');
                    if (i > 0) line = line.substring(i);
                    TreeParser tree = new TreeParser(line);
                    NetworkParser network = new NetworkParser(tree);

                    final int nHybrid = network.getReticulationNodeCount();
                    nHybridInNetworkMap.put(nHybrid, network);

                    numNetworks++;
                }
            }

            writeOutput(outFileName);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void writeOutput(String outputFileName) throws IOException {
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

        for (Integer i = 0; i < nHybridInNetworkMap.keySet().size(); i++) {
            out.println(i + "-->" + nHybridInNetworkMap.get(i).size());
        }
    }
}
