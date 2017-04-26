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
		storedEmbedding = new int[nodeCount][];
    }

	@Override
	public void initAndValidate() {
		super.initAndValidate();
		embedding = new int[nodeCount][];
		storedEmbedding = new int[nodeCount][];
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
		//printEmbedding();
		return embedding[row][col];
	}

	/* private void printEmbedding() {
		for (int r = 0; r < embedding.length; r++) {
			final StringBuffer rowStr = new StringBuffer();
			final int rowLen = embedding[r].length;
			for (int c = 0; c < rowLen - 1; c++) {
				rowStr.append(embedding[r][c]);
				rowStr.append(" ");
			}
			rowStr.append(embedding[r][rowLen - 1]);
			System.out.println(rowStr);
		}
	} */

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
			final int geneNodeCount = src.length;
			final int traversalNodeCount = src[0].length;
			final boolean reinitialize = dst[0] == null || dst[0].length != traversalNodeCount;
			for (int i = 0; i < geneNodeCount; i++) {
				if (reinitialize) dst[i] = new int[traversalNodeCount];
				System.arraycopy(src[i], 0, dst[i], 0, traversalNodeCount);
			}
		}
    }

    @Override
    public String toString() {
    	for (int i = 0; i < nodeCount; i++) {
    		final int lastTraversalNode = embedding[i].length - 1;
    		final StringBuffer embeddingBuf = new StringBuffer();
    		embeddingBuf.append("embedding=\"");
    		for (int j = 0; j < lastTraversalNode; j++) {
    			embeddingBuf.append(embedding[i][j]);
    			embeddingBuf.append(' ');
    		}
			embeddingBuf.append(embedding[i][lastTraversalNode]);
			embeddingBuf.append('"');
    		m_nodes[i].metaDataString = embeddingBuf.toString();
    	}
    	return super.toString();
    }

    @Override
    public void fromXML(final org.w3c.dom.Node node) {
    	super.fromXML(node);

    	embedding = new int[nodeCount][];

    	for (int i = 0; i < nodeCount; i++) {
    		final String embedStr = (String) m_nodes[i].getMetaData("embedding");
    		// adapted from stackoverflow.com/a/7781050
    		String[] parts = embedStr.split(" ");
    		final int[] n1 = new int[parts.length];
    		for(int n = 0; n < parts.length; n++) n1[n] = Integer.parseInt(parts[n]);
    		embedding[i] = n1;
    	}
    }
}
