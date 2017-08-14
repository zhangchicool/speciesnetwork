package speciesnetwork;

public class Embedding {
	protected int geneNodeCount;
	protected int traversalNodeCount;
	protected int[] embedding;
	protected double probability;

	public Embedding(int gnc) {
		geneNodeCount = gnc;
		traversalNodeCount = 1;
		embedding = new int[gnc];
	}

	public Embedding(int gnc, int tnc) {
		geneNodeCount = gnc;
		traversalNodeCount = tnc;
		embedding = new int[geneNodeCount * traversalNodeCount];
	}

	public Embedding(Embedding src) {
		geneNodeCount = src.geneNodeCount;
		traversalNodeCount = src.traversalNodeCount;
		if (src.embedding != null) {
			embedding = new int[src.embedding.length];
			System.arraycopy(src.embedding, 0, embedding, 0, embedding.length);
		}
		probability = src.probability;
	}

	public int getDirection(int geneNode, int traversalNode) {
		final int i = (traversalNodeCount * geneNode) + traversalNode;
		return embedding[i];
	}

	public void setDirection(int geneNode, int traversalNode, int value) {
		final int i = (traversalNodeCount * geneNode) + traversalNode;
		embedding[i] = value;
	}

	public void reset(int tnc, int value) {
		if (traversalNodeCount != tnc) {
			traversalNodeCount = tnc;
			embedding = new int[geneNodeCount * traversalNodeCount];
		}

		java.util.Arrays.fill(embedding, value);
	}

	// assumes geneNodeCount is unchanged
	public void copyFrom(Embedding src) {
		if (src.traversalNodeCount != traversalNodeCount) {
			traversalNodeCount = src.traversalNodeCount;
			embedding = new int[src.embedding.length];
		}

		System.arraycopy(src.embedding, 0, embedding, 0, embedding.length);
	}

	public void appendRowToString(StringBuilder embeddingBuf, int row) {
		final int offset = row * traversalNodeCount;
		embeddingBuf.append(offset);
		for (int i = 1; i < geneNodeCount; i++) {
			embeddingBuf.append(' ');
			embeddingBuf.append(offset + i);
		}
	}
}
