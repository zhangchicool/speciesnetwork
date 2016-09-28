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
import beast.evolution.branchratemodel.BranchRateModel;

@Description("Based on the TreeWithMetaDataLogger class, but with support for population sizes")
public class SpeciesNetworkLogger extends BEASTObject implements Loggable {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network to be logged.", Validate.REQUIRED);
    public final Input<PopulationSizeModel> populationModelInput =
            new Input<>("populationmodel", "Population sizes to be logged with branches of the tree");
    // TODO: make this input a list of valuables
    public final Input<List<Function>> parameterInput =
            new Input<>("metadata", "meta data to be logged with the tree nodes",new ArrayList<>());
    public final Input<BranchRateModel> clockModelInput =
            new Input<>("branchratemodel", "Rate to be logged with branches of the tree");
    public final Input<Boolean> substitutionsInput = new Input<>("substitutions",
            "Report branch lengths as substitutions (branch length times clock rate, default false)", false);
    public final Input<Integer> decimalPlacesInput = new Input<>("dp",
            "The number of decimal places to use writing branch lengths and rates (default -1 for full precision)", -1);

    boolean someMetaDataNeedsLogging;
    boolean substitutions = false;

    private DecimalFormat df;

    @Override
    public void initAndValidate() {
        if (parameterInput.get().size() == 0 && clockModelInput.get() == null && populationModelInput.get() == null) {
            someMetaDataNeedsLogging = false;
            return;
        }
        someMetaDataNeedsLogging = true;
        // without substitution model, reporting substitutions == reporting branch lengths 
        if (clockModelInput.get() != null) {
            substitutions = substitutionsInput.get();
        }

        int dp = decimalPlacesInput.get();

        if (dp < 0) {
            df = null;
        } else {
            // just new DecimalFormat("#.######") (with dp time '#' after the decimal)
            df = new DecimalFormat("#."+new String(new char[dp]).replace('\0', '#'));
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
        // ??? network.getRoot().sort();
        // out.print(toNewick(network.getRoot(), metadata, branchRateModel, populationModel));
        // out.print(tree.getRoot().toShortNewick(false));
        out.print(network.toString());
        out.print(";");
    }

    @Override
    public void close(PrintStream out) {
        Network speciesNetwork = speciesNetworkInput.get();
        speciesNetwork.close(out);
    }
}
