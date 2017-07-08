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
        final int[][] embedding = geneTreeInput.get().embedding;
        int row = embedding.length;
        int col = embedding[0].length;

        for (int r = 0; r < row; r++)
            for (int c = 0; c < col; c++)
                out.print(geneTreeInput.get().getID() + "_" + r + "_" + c + "\t");
    }

    @Override
    public void log(int sample, PrintStream out) {
        final int[][] embedding = geneTreeInput.get().embedding;
        int row = embedding.length;
        int col = embedding[0].length;

        for (int r = 0; r < row; r++)
            for (int c = 0; c < col; c++)
                out.print(embedding[r][c] + "\t");
    }

	@Override
    public void close(PrintStream out) {
        // nothing to do
    }

    @Override
    public int getDimension() {
        final int[][] embedding = geneTreeInput.get().embedding;
        return embedding.length * embedding[0].length;
    }

    @Override
    public double getArrayValue() {
        return -1;
    }

    @Override
    public double getArrayValue(int dim) {
        final int[][] embedding = geneTreeInput.get().embedding;
        int row = dim / embedding[0].length;
        int col = dim % embedding[0].length;
        return embedding[row][col];
    }
}
