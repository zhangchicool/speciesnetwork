package speciesnetwork.tools;

import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.*;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.util.Log;
import beast.core.Runnable;
import beast.util.HeapSort;
import beast.util.TreeParser;
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
        // print header
        out.println("length  height  nHybrid  tHybrid");
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

                    if (numNetworks > burnin) {
                        printSummary (network, out);
                    }
                    if (numNetworks % 10000 == 0) progressStream.print(".");
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        progressStream.println("\nParsed " + numNetworks + " networks totally, " + burnin + " discarded as burn-in.");

        out.close();
    }

    private void printSummary (Network network, PrintStream out){
        double rootHeight = network.getRoot().getHeight();
        double originHeight = network.getOrigin().getHeight();
        double length = network.getNetworkLength() - (originHeight - rootHeight);

        int nHybrid = network.getReticulationNodeCount();

        double tHybrid = 0.0;
        for (NetworkNode hNode: network.getReticulationNodes()) {
            if (tHybrid == 0.0 || tHybrid > hNode.getHeight())
                tHybrid = hNode.getHeight();
        }

        out.println(length + "\t" + rootHeight + "\t" + nHybrid + "\t" + tHybrid);
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
