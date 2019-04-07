package speciesnetwork;

import java.io.PrintStream;

import org.w3c.dom.Node;

import beast.core.StateNode;
import java.util.ArrayList;

public class Embedding {
	protected int geneNodeCount;
	protected int traversalNodeCount;
	protected int[] embedding;
	public double probability = 1.0;
	public double probabilitySum = 1.0;
		
	public Embedding() {
		geneNodeCount = 0;
		traversalNodeCount = 1;
		embedding = new int[geneNodeCount * traversalNodeCount];	
    java.util.Arrays.fill(embedding, -1);
	}

	public Embedding(int gnc) {
		geneNodeCount = gnc;
		traversalNodeCount = 1;
		embedding = new int[geneNodeCount * traversalNodeCount];
    java.util.Arrays.fill(embedding, -1);
	}

	public Embedding(int gnc, int tnc) {
		geneNodeCount = gnc;
		traversalNodeCount = tnc;
		embedding = new int[geneNodeCount * traversalNodeCount];
    java.util.Arrays.fill(embedding, -1);
	}

	public int getDirection(int geneNode, int traversalNode) {
		return embedding[geneNode * traversalNodeCount + traversalNode];
	}

	public void setDirection(int geneNode, int traversalNode, int value) {
		embedding[geneNode * traversalNodeCount + traversalNode] = value;
	}

	public void reset(int tnc) {
		if (traversalNodeCount != tnc) {
			traversalNodeCount = tnc;
			embedding = new int[geneNodeCount * traversalNodeCount];
		}
    java.util.Arrays.fill(embedding, -1);
	}

	public void mergeWith(Embedding src) {
		assert src.geneNodeCount == geneNodeCount;
		assert src.traversalNodeCount == traversalNodeCount;

		probability *= src.probability;
		probabilitySum *= src.probabilitySum;
		for (int g = 0; g < geneNodeCount; g++) {
			for (int t = 0; t < traversalNodeCount; t++) {
				if (embedding[g * traversalNodeCount + t] == -1 && src.embedding[g * traversalNodeCount + t] != -1) {
					embedding[g * traversalNodeCount + t] = src.embedding[g * traversalNodeCount + t];
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
			str.append(embedding[row * traversalNodeCount + t]);
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

	public void log(long sample, PrintStream out) {
		out.append(toString() + "\t");
	}

	public void close(PrintStream out) {
	}

	public int getDimension() {
		return geneNodeCount * traversalNodeCount;
	}

	public double getArrayValue() {
		// This method makes no sense.
		throw new RuntimeException("Useless call to Embedding.getArrayValue");
	}

	public double getArrayValue(int dim) {
		// This method would make slightly more sense as int method.
		return embedding[dim];
	}

	public Embedding copy() {
		Embedding e = new Embedding(geneNodeCount, traversalNodeCount);

		e.probability = probability;
		e.probabilitySum = probabilitySum;
		e.reset(traversalNodeCount);
		e.mergeWith(this);
		return e;
	}

	public void assignFrom(Object other) {
		if (other instanceof Embedding) {
			Embedding e = (Embedding) other;
			
			geneNodeCount = e.geneNodeCount;
			traversalNodeCount = e.traversalNodeCount;
			embedding = new int[geneNodeCount * traversalNodeCount];
			System.arraycopy(e.embedding, 0, embedding, 0, geneNodeCount * traversalNodeCount);
			probability = e.probability;
			probabilitySum = e.probabilitySum;
		} else {
			throw new RuntimeException("Useless call to Embedding.assignFrom");
		}
	}

	public void fromXML(Node node) {
		String[] tNodes = node.getTextContent().split("//");
		geneNodeCount = tNodes.length;
    ArrayList<Integer> em = new ArrayList<>();
    int i = 0;
		for (int g = 0; g < geneNodeCount; ++g) {
			String row = tNodes[g];
			String[] values = row.trim().split("\\s");
      if (g == 0) {
        traversalNodeCount = values.length;
        embedding = new int[geneNodeCount * traversalNodeCount];
      }
			for (String v: values) {
				embedding[i] = Integer.parseInt(v);
        ++i;
			}
		}
	}
}
