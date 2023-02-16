package speciesnetwork.utils;

import java.io.PrintStream;

import beast.base.inference.CalculationNode;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import speciesnetwork.EmbeddedTree;
import speciesnetwork.Embedding;

@Description("Logs embedding of a gene tree")
public class EmbeddingLogger extends CalculationNode implements Loggable, Function {
    public final Input<EmbeddedTree> geneTreeInput =
            new Input<>("geneTree", "Gene tree embedded in the species network.", Validate.REQUIRED);

    @Override
    public void initAndValidate() {
    }

    @Override
    public void init(PrintStream out) {
        final Embedding embedding = geneTreeInput.get().embedding;
        for (int r = 0; r < embedding.geneNodeCount; r++)
            for (int c = 0; c < embedding.traversalNodeCount; c++)
                out.print(geneTreeInput.get().getID() + "_" + r + "_" + c + "\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        final Embedding embedding = geneTreeInput.get().embedding;
        out.print(embedding.toString());
    }

	@Override
    public void close(PrintStream out) {
        // nothing to do
    }

    @Override
    public int getDimension() {
        final Embedding embedding = geneTreeInput.get().embedding;
        return embedding.geneNodeCount * embedding.traversalNodeCount;
    }

    @Override
    public double getArrayValue() {
        return -1;
    }

    @Override
    public double getArrayValue(int i) {
        final Embedding embedding = geneTreeInput.get().embedding;
        final int[] e = embedding.getEmbedding();
        return e[i];
    }
}
