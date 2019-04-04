package speciesnetwork;

import java.io.PrintStream;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;

@Description("Logs embedding of a gene tree")
public class EmbeddingLogger extends CalculationNode implements Loggable {
    public final Input<GeneTreeInSpeciesNetwork> geneTreeInput =
            new Input<>("geneTree", "Gene tree embedded in the species network.", Validate.REQUIRED);

    @Override
    public void initAndValidate() {
    }

    @Override
    public void init(PrintStream out) {
        final Embedding embedding = geneTreeInput.get().embeddingInput.get();
        for (int r = 0; r < embedding.geneNodeCount; r++)
            for (int c = 0; c < embedding.traversalNodeCount; c++)
                out.print(geneTreeInput.get().getID() + "_" + r + "_" + c + "\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        final Embedding embedding = geneTreeInput.get().embeddingInput.get();
        out.print(embedding.toString());
    }

	@Override
    public void close(PrintStream out) {
        // nothing to do
    }
}
