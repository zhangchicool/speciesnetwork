package speciesnetwork;

public class Embedding {
	protected int geneNodeCount;
	protected int traversalNodeCount;
	protected int[][] embedding;
	public double probability = 1.0;
	public double probabilitySum = 1.0;

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

	public Embedding(Embedding src) {
		geneNodeCount = src.geneNodeCount;
		traversalNodeCount = src.traversalNodeCount;
		embedding = new int[src.embedding.length][src.embedding[0].length];
		System.arraycopy(src.embedding, 0, embedding, 0, embedding.length);
		probability = src.probability;
		probabilitySum = src.probabilitySum;
	}

	public int getGenes() {
		return geneNodeCount;
	}

	public int getTravs() {
		return traversalNodeCount;
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
		java.util.Arrays.fill(embedding, -1);
	}

	// assumes geneNodeCount is unchanged
	public void copyFrom(Embedding src) {
		geneNodeCount = src.geneNodeCount;
		traversalNodeCount = src.traversalNodeCount;
		embedding = new int[geneNodeCount][traversalNodeCount];
		for (int g = 0; g < geneNodeCount; g++) {
			System.arraycopy(src.embedding[g], 0, embedding[g], 0, embedding[g].length);
		}
		probability = src.probability;
		probabilitySum = src.probabilitySum;
	}

	public void mergeWith(Embedding src) {
		assert src.geneNodeCount == geneNodeCount;
		assert src.traversalNodeCount == traversalNodeCount;

		probability *= src.probability;
		probabilitySum *= src.probabilitySum;
		for (int g = 0; g < geneNodeCount; g++) {
			for (int t = 0; t < traversalNodeCount; t++) {
				if (embedding[g][t] == -1) {
					embedding[g][t] = src.embedding[g][t];
				}
			}
		}
	}

	public String rowToString(int row) {
		StringBuilder str = new StringBuilder();
		for (int t = 0; t < traversalNodeCount; t++) {
			str.append(' ');
			str.append(embedding[row][t]);
		}
		return str.toString();
	}

	@Override
	public String toString() {
		StringBuilder str = new StringBuilder();
		for (int g = 0; g < geneNodeCount; g++) {
			for (int t = 0; t < traversalNodeCount; t++) {
				str.append(' ');
				str.append(embedding[g][t]);
			}
		}
		return str.toString();
	}
}
