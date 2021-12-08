package speciesnetwork.utils;

import java.io.PrintStream;
import java.math.RoundingMode;
import java.text.DecimalFormat;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import speciesnetwork.Network;

@Description("Logs backbone tree annotated with metadata")
public class ReticulationLogger extends BEASTObject implements Loggable {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network to be logged.", Validate.REQUIRED);


    @Override
    public void initAndValidate() {
    	
    }

    @Override
    public void init(PrintStream out) {
        final Network speciesNetwork = speciesNetworkInput.get();
        out.print(speciesNetwork.getID() + ".reticulationNodes\t");
        
    }

    @Override
    public void log(long sample, PrintStream out) {
        final Network speciesNetwork = speciesNetworkInput.get();
    	out.print(speciesNetwork.getReticulationNodeCount() + "\t");
    }

    @Override
    public void close(PrintStream out) {
    }
}
