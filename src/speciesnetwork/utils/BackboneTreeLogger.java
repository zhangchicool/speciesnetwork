package speciesnetwork.utils;

import java.io.PrintStream;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.StateNode;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;
import speciesnetwork.NodeHeightComparator;

@Description("Logs backbone tree annotated with metadata")
public class BackboneTreeLogger extends BEASTObject implements Loggable {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network to be logged.", Validate.REQUIRED);
    public final Input<Integer> decimalPlacesInput = new Input<>("dp",
            "The number of decimal places to use (default -1 for full precision)", -1);

    private DecimalFormat df;
    private static Comparator<NetworkNode> hc = new NodeHeightComparator();

    @Override
    public void initAndValidate() {
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
    public void init(PrintStream out) {
        Network speciesNetwork = speciesNetworkInput.get();
        speciesNetwork.init(out);
    }

    @Override
    public void log(long sample, PrintStream out) {
        // make sure we get the current version of the inputs
        Network network = (Network) speciesNetworkInput.get().getCurrent();
        Network backbone = getBackboneTree(network);

        // write out the backbone tree with meta data
        out.print("tree STATE_" + sample + " = ");
        out.print(backbone.toString(df));
        out.print(";");
    }

    @Override
    public void close(PrintStream out) {
        speciesNetworkInput.get().close(out);
    }

    private Network getBackboneTree (Network network) {
        // make a copy of the current network to modify
        Network backbone = new Network();
        backbone.assignFrom(network);

        // get and sort the hybridization nodes
        List<NetworkNode> hNodes = Arrays.asList(backbone.getReticulationNodes());
        hNodes.sort(hc);

        // start deleting from the oldest to the youngest
        // this guarantees every hybridization branch can be safely deleted
        for (int i = hNodes.size() - 1; i >= 0; i--) {
            final NetworkNode node = hNodes.get(i);
            if (node.getGammaProb() < 0.5) {
                // delete the gamma branch
                backbone.deleteReticulationBranch(node.gammaBranchNumber);
            }
            else {
                // delete the alternative branch
                backbone.deleteReticulationBranch(node.gammaBranchNumber + 1);
            }

        }

        return backbone;
    }
}
