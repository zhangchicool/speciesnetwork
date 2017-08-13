package speciesnetwork;

import java.io.PrintStream;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.util.Log;
import beast.core.Loggable;

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
        final int row = embedding.geneNodeCount;
        final int col = embedding.traversalNodeCount;

        for (int r = 0; r < row; r++)
            for (int c = 0; c < col; c++)
                out.print(geneTreeInput.get().getID() + "_" + r + "_" + c + "\t");
    }

    @Override
    public void log(int sample, PrintStream out) {
        final Embedding embedding = geneTreeInput.get().embedding;
        final int row = embedding.geneNodeCount;
        final int col = embedding.traversalNodeCount;

        for (int r = 0; r < row; r++)
            for (int c = 0; c < col; c++)
                out.print(embedding.getDirection(r, c) + "\t");
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
        return embedding.embedding[i];
    }
}
