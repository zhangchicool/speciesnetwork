package speciesnetwork;

public class Embedding {
	protected int geneNodeCount;
	protected int traversalNodeCount;
	protected int[] embedding;
	public double probability;

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
		embedding = new int[src.embedding.length];
		System.arraycopy(src.embedding, 0, embedding, 0, embedding.length);
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
		if (traversalNodeCount != src.traversalNodeCount) {
			traversalNodeCount = src.traversalNodeCount;
			embedding = new int[src.embedding.length];
		}

		System.arraycopy(src.embedding, 0, embedding, 0, embedding.length);
		probability = src.probability;
	}

	public String rowToString(int row) {
        StringBuilder str = new StringBuilder();
        final int offset = row * traversalNodeCount;
		str.append(embedding[offset]);
		for (int i = 1; i < traversalNodeCount; i++) {
			str.append(' ');
			str.append(embedding[offset + i]);
		}
		return str.toString();
	}

    @Override
	public String toString() {
        StringBuilder str = new StringBuilder();
        str.append(embedding[0]);
	    for (int i = 1; i < geneNodeCount * traversalNodeCount; i++) {
            str.append(' ');
            str.append(embedding[i]);
        }
        return str.toString();
    }
}
