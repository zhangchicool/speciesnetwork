package speciesnetwork;

import beast.core.StateNode;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

public class EmbeddedTree extends Tree implements EmbeddableTree {
	int[][] embedding;
	int[][] storedEmbedding;

	public EmbeddedTree() {
	}

	// based on Tree(final Node rootNode)
    public EmbeddedTree(final Node rootNode) {
        setRoot(rootNode);
        initArrays();
		embedding = new int[nodeCount][];
    }

	@Override
	public void initAndValidate() {
		super.initAndValidate();
		embedding = new int[nodeCount][];
	}

	@Override
	public void resetEmbedding(int nCol, int value) {
		final boolean reinitialize = embedding[0] == null || embedding[0].length != nCol;
		for (int r = 0; r < nodeCount; r++) {
			if (reinitialize) embedding[r] = new int[nCol];
			java.util.Arrays.fill(embedding[r], value);
		}
	}

	@Override
	public void setEmbedding(int row, int col, int value) {
		embedding[row][col] = value;
	}

	@Override
	public int getEmbedding(int row, int col) {
		return embedding[row][col];
	}

	// based on Tree.copy()
    @Override
    public EmbeddedTree copy() {
        EmbeddedTree etree = new EmbeddedTree();
        etree.setID(getID());
        etree.index = index;
        etree.root = root.copy();
        etree.nodeCount = nodeCount;
        etree.internalNodeCount = internalNodeCount;
        etree.leafNodeCount = leafNodeCount;
    	copyEmbedding(embedding, etree.embedding);

    	return etree;
    }

    @Override
    public void assignTo(final StateNode other) {
    	super.assignTo(other);
        final EmbeddedTree etree = (EmbeddedTree) other;
        copyEmbedding(embedding, etree.embedding);
    }

    @Override
    public void assignFrom(final StateNode other) {
    	super.assignFrom(other);
        final EmbeddedTree etree = (EmbeddedTree) other;
        copyEmbedding(etree.embedding, embedding);
    }

    @Override
    public void assignFromFragile(final StateNode other) {
    	super.assignFromFragile(other);
        final EmbeddedTree etree = (EmbeddedTree) other;
        copyEmbedding(etree.embedding, embedding);
    }

    @Override
    public void store() {
    	super.store();
    	copyEmbedding(embedding, storedEmbedding);
    }

    @Override
    public void restore() {
    	super.restore();
    	final int[][] tmpEmbedding = embedding;
    	embedding = storedEmbedding;
    	storedEmbedding = tmpEmbedding;
    }

    static void copyEmbedding(int[][] src, int[][] dst) {
		if (dst == null) dst = new int[src.length][];

		if (src[0] != null) {
			final int traversalNodeCount = src[0].length;
			final boolean reinitialize = dst[0] == null || dst[0].length != src[0].length;
			for (int i = 0; i < src.length; i++) {
				if (reinitialize) dst[i] = new int[traversalNodeCount];
				System.arraycopy(src[i], 0, dst[i], 0, traversalNodeCount);
			}
		}
    }
}
