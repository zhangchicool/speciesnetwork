package speciesnetwork.utils;

import java.io.PrintStream;

import beast.base.inference.CalculationNode;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.core.Loggable;
import speciesnetwork.Network;

@Description("Logs statistics of a species network")
public class NetworkStatLogger extends CalculationNode implements Loggable, Function {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network to be logged.", Validate.REQUIRED);
    public final Input<Boolean> logHeightInput = new Input<>("logHeight", "If true, root height will be logged.", true);
    public final Input<Boolean> logLengthInput = new Input<>("logLength", "If true, network length will be logged.", true);

    @Override
    public void initAndValidate() {
    	if (!logHeightInput.get() && !logLengthInput.get()) {
    		Log.warning.println("NetworkStatLogger " + getID() + " logs nothing. Set logHeigth=true or logLength=true");
    	}
    }

    @Override
    public void init(PrintStream out) {
        final Network speciesNetwork = speciesNetworkInput.get();
        if (logHeightInput.get()) {
            out.print(speciesNetwork.getID() + ".height\t");
        }
        if (logLengthInput.get()) {
            out.print(speciesNetwork.getID() + ".length\t");
        }
    }

    @Override
    public void log(long sample, PrintStream out) {
        final Network speciesNetwork = speciesNetworkInput.get();
        if (logHeightInput.get()) {
            // network root height
        	out.print(speciesNetwork.getRoot().getHeight() + "\t");
        }
        if (logLengthInput.get()) {
            // network length (excluding the root branch length)
            final double rootBrl = speciesNetwork.getOrigin().getHeight() - speciesNetwork.getRoot().getHeight();
            out.print(speciesNetwork.getNetworkLength() - rootBrl + "\t");
        }
    }

	@Override
    public void close(PrintStream out) {
        // nothing to do
    }

    @Override
    public int getDimension() {
        return 2;
    }

    @Override
    public double getArrayValue() {
        return speciesNetworkInput.get().getRoot().getHeight();
    }

    @Override
    public double getArrayValue(int dim) {
        if (dim == 0)
            return speciesNetworkInput.get().getRoot().getHeight();
        else
            return speciesNetworkInput.get().getNetworkLength();
    }
}
