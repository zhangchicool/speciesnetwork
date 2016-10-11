package speciesnetwork;

import java.io.*;
import java.util.*;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Runnable;

/**
 * @author Chi Zhang
 */

@Description("Simulate gene trees given a species network (multispecies coalescent).")
public class SummarizeNetwork extends Runnable {
    public final Input<String> fileNameInput = new Input<>("fileName",
            "Name of the file that contains networks in extended newick format.", Validate.REQUIRED);

    @Override
    public void initAndValidate() {

    }

    @Override
    public void run() throws IOException {

        try (BufferedReader br = new BufferedReader(new FileReader(fileNameInput.get()))) {
            String line;
            while ((line = br.readLine()) != null) {
                // process the line.
            }
        }
    }
}
