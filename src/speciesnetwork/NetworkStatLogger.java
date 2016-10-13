package speciesnetwork;

import java.io.PrintStream;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.util.Log;
import beast.core.Loggable;


@Description("Logger to report statistics of a network")
public class NetworkStatLogger extends CalculationNode implements Loggable, Function {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network to be logged.", Validate.REQUIRED);
    public final Input<Boolean> logHeigthInput = new Input<>("logHeigth", "If true, MRCA height will be logged.", true);
    public final Input<Boolean> logLengthInput = new Input<>("logLength", "If true, network length will be logged.", true);

    @Override
    public void initAndValidate() {
    	if (!logHeigthInput.get() && !logLengthInput.get()) {
    		Log.warning.println("NetworkStatLogger " + getID() + " logs nothing. Set logHeigth=true or logLength=true");
    	}
    }

    @Override
    public void init(PrintStream out) {
        final Network speciesNetwork = speciesNetworkInput.get();
        if (logHeigthInput.get()) {
            out.print(speciesNetwork.getID() + ".tmrca\t");
        }
        if (logLengthInput.get()) {
            out.print(speciesNetwork.getID() + ".length\t");
        }
    }

    @Override
    public void log(int sample, PrintStream out) {
        final Network speciesNetwork = speciesNetworkInput.get();
        if (logHeigthInput.get()) {
        	out.print(speciesNetwork.getRoot().getHeight() + "\t");
        }
        if (logLengthInput.get()) {
            out.print(speciesNetwork.getNetworkLength() + "\t");
        }
    }

	@Override
    public void close(PrintStream out) {
        // nothing to do
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getArrayValue() {
        return speciesNetworkInput.get().getRoot().getHeight();
    }

    @Override
    public double getArrayValue(int dim) {
        return speciesNetworkInput.get().getRoot().getHeight();
    }
}
