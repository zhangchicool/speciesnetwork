package speciesnetwork.simulator;

import java.io.*;
import java.util.*;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Runnable;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import speciesnetwork.Network;
import speciesnetwork.NetworkNode;

/**
 * @author Chi Zhang
 */

@Description("Simulate a species network under the pure birth and hybridization process.")
public class NetworkPureBirthHybrid extends Runnable {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "Species network to be simulated.", Validate.REQUIRED);
    public final Input<RealParameter> originInput =
            new Input<>("origin", "The time when the process started.", Validate.REQUIRED);

    @Override
    public void initAndValidate() {
    }

    @Override
    public void run() throws IOException {
        Network speciesNetwork = speciesNetworkInput.get();

    }


}
