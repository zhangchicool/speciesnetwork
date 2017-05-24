package speciesnetwork;

import java.io.PrintStream;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.StateNode;

@Description("Based on the TreeWithMetaDataLogger class, but with support for network")
public class NetworkWithMetaDataLogger extends BEASTObject implements Loggable {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network to be logged.", Validate.REQUIRED);
    public final Input<PopulationSizeModel> populationModelInput =
            new Input<>("populationmodel", "Population sizes to be logged for branches.");
    public final Input<List<Function>> parameterInput =
            new Input<>("metadata", "meta data to be logged with the nodes.",new ArrayList<>());
    public final Input<Integer> decimalPlacesInput = new Input<>("dp",
            "The number of decimal places to use (default -1 for full precision)", -1);

    private DecimalFormat df;

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
    public void log(int nSample, PrintStream out) {
        // make sure we get the current version of the inputs
        Network network = (Network) speciesNetworkInput.get().getCurrent();
        List<Function> metadata = parameterInput.get();
        for (int i = 0; i < metadata.size(); i++) {
            if (metadata.get(i) instanceof StateNode) {
                metadata.set(i, ((StateNode) metadata.get(i)).getCurrent());
            }
        }
        // BranchRateModel branchRateModel = clockModelInput.get();
        // PopulationSizeModel populationModel = populationModelInput.get();
        // write out the log tree with meta data
        out.print("tree STATE_" + nSample + " = ");
        out.print(network.toString(df));
        out.print(";");
    }

    @Override
    public void close(PrintStream out) {
        Network speciesNetwork = speciesNetworkInput.get();
        speciesNetwork.close(out);
    }
}
