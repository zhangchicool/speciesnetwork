package speciesnetwork;

import java.io.PrintStream;

import org.w3c.dom.Node;

import beast.core.StateNode;
import beast.util.TreeParser;

public class Embedding extends StateNode {
	protected int geneNodeCount;
	protected int traversalNodeCount;
	protected int[][] embedding;
	public double probability = 1.0;
	public double probabilitySum = 1.0;
	
	Embedding stored = null;
	
	public Embedding() {
		geneNodeCount = 0;
		traversalNodeCount = 1;
		embedding = new int[geneNodeCount][traversalNodeCount];		
	}

	public Embedding(int gnc) {
		geneNodeCount = gnc;
		traversalNodeCount = 1;
		embedding = new int[geneNodeCount][traversalNodeCount];
		for (int[] row : embedding) {
			java.util.Arrays.fill(row, -1);
		}
	}

	public Embedding(int gnc, int tnc) {
		geneNodeCount = gnc;
		traversalNodeCount = tnc;
		embedding = new int[geneNodeCount][traversalNodeCount];
		for (int[] row : embedding) {
			java.util.Arrays.fill(row, -1);
		}
	}

	public int getDirection(int geneNode, int traversalNode) {
		return embedding[geneNode][traversalNode];
	}

	public void setDirection(int geneNode, int traversalNode, int value) {
		embedding[geneNode][traversalNode] = value;
	}

	public void reset(int tnc) {
		if (traversalNodeCount != tnc) {
			traversalNodeCount = tnc;
			embedding = new int[geneNodeCount][traversalNodeCount];
		}
		for (int g = 0; g < geneNodeCount; g++) {
			java.util.Arrays.fill(embedding[g], -1);
		}
	}

	public void mergeWith(Embedding src) {
		assert src.geneNodeCount == geneNodeCount;
		assert src.traversalNodeCount == traversalNodeCount;

		probability *= src.probability;
		probabilitySum *= src.probabilitySum;
		for (int g = 0; g < geneNodeCount; g++) {
			for (int t = 0; t < traversalNodeCount; t++) {
				if (embedding[g][t] == -1 && src.embedding[g][t] != -1) {
					embedding[g][t] = src.embedding[g][t];
				}
			}
		}
	}

	public String rowToString(int row) {
		StringBuilder str = new StringBuilder();
		for (int t = 0; t < traversalNodeCount; t++) {
			if (t > 0) {
				str.append(' ');
			}
			str.append(embedding[row][t]);
		}
		return str.toString();
	}

	@Override
	public String toString() {
		StringBuilder str = new StringBuilder();
		for (int g = 0; g < geneNodeCount; g++) {
			if (g > 0) {
				str.append("//");
			}
			str.append(rowToString(g));
		}
		return str.toString();
	}

	@Override
	public void init(PrintStream out) {
		out.append(getID() + "\t");
	}

	@Override
	public void log(long sample, PrintStream out) {
		out.append(toString() + "\t");
	}

	@Override
	public void close(PrintStream out) {
	}

	@Override
	public int getDimension() {
		return geneNodeCount * traversalNodeCount;
	}

	@Override
	public double getArrayValue() {
		// This method makes no sense.
		throw new RuntimeException("Useless call to Embedding.getArrayValue");
	}

	@Override
	public double getArrayValue(int dim) {
		// This method would make slightly more sense as int method.
		return embedding[dim / traversalNodeCount][dim % traversalNodeCount];
	}

	@Override
	public void initAndValidate() {
	}

	@Override
	public Embedding copy() {
		Embedding e = new Embedding(geneNodeCount, traversalNodeCount);

		e.probability = probability;
		e.probabilitySum = probabilitySum;
		e.reset(traversalNodeCount);
		e.mergeWith(this);
		return e;
	}

	@Override
	public void assignTo(StateNode other) {
		throw new RuntimeException("Useless call to Embedding.assignTo");
	}

	@Override
	public void assignFrom(StateNode other) {
		if (other instanceof Embedding) {
			Embedding e = (Embedding) other;
			
			geneNodeCount = e.geneNodeCount;
			traversalNodeCount = e.traversalNodeCount;
			embedding = new int[geneNodeCount][traversalNodeCount];
			for (int g = 0; g < geneNodeCount; g++) {
				System.arraycopy(e.embedding[g], 0, embedding[g], 0, embedding[g].length);
			}
			probability = e.probability;
			probabilitySum = e.probabilitySum;
		} else {
			throw new RuntimeException("Useless call to Embedding.assignFrom");
		}
	}

	@Override
	public void assignFromFragile(StateNode other) {
		throw new RuntimeException("Useless call to Embedding.assignFromFragile");
	}

	@Override
	public void fromXML(Node node) {
		String[] tNodes = node.getTextContent().split("//");
		geneNodeCount = tNodes.length;
		traversalNodeCount = 0;
		embedding = new int[geneNodeCount][];
		for (int g = 0; g < geneNodeCount; ++g) {
			String row = tNodes[g];
			String[] values = row.trim().split("\\s");
			if (traversalNodeCount == 0) {
				traversalNodeCount = values.length;
			} else {
				if (traversalNodeCount != values.length) {
					throw new RuntimeException("Embedding is not rectangular");
				}
			}
			embedding[g] = new int[traversalNodeCount];
			for (int t = 0; t < traversalNodeCount; ++t) {
				embedding[g][t] = Integer.parseInt(values[t]);
			}
		}
	}

	@Override
	public int scale(double scale) {
		throw new RuntimeException("Useless call to Embedding.scale");
	}

	@Override
	protected void store() {
		stored = copy();
	}

	@Override
	public void restore() {
		assignFrom(stored);
	}

	@Override
	public void setEverythingDirty(boolean isDirty) {
		setSomethingIsDirty(isDirty);
	}
}
