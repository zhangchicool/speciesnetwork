package speciesnetwork.utils;

import java.io.PrintStream;
import java.math.RoundingMode;
import java.text.DecimalFormat;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import speciesnetwork.Network;

@Description("Logs backbone tree annotated with metadata")
public class BackboneTreeLogger extends BEASTObject implements Loggable {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network to be logged.", Validate.REQUIRED);
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
    public void log(long sample, PrintStream out) {
        // make sure we get the current version of the inputs
        Network network = (Network) speciesNetworkInput.get().getCurrent();
        Network backbone = network.getBackboneTree();

        // write out the backbone tree with meta data
        out.print("tree STATE_" + sample + " = ");
        out.print(backbone.toString(df));
        out.print(";");
    }

    @Override
    public void close(PrintStream out) {
        speciesNetworkInput.get().close(out);
    }
}
