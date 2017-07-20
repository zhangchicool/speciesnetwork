package speciesnetwork;

import beast.core.StateNode;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

public class EmbeddedTree extends Tree {
	protected int[][] embedding;
	private int[][] storedEmbedding;

	// for the RebuildEmbedding operator
	public int embeddingCount;  // the number of possible embeddings
    protected int storedEmbeddingCount;
    public int choicesCount;  // the number of alternative traversing choices
    protected int storedChoicesCount;

    public EmbeddedTree() {
    }

	// based on Tree(final Node rootNode)
    public EmbeddedTree(final Node rootNode) {
        setRoot(rootNode);
        initArrays();
        embedding = new int[nodeCount][];
        storedEmbedding = new int[nodeCount][];
        embeddingCount = -1;
        storedEmbeddingCount = -1;
        choicesCount = -1;
        storedChoicesCount = -1;
    }

	@Override
	public void initAndValidate() {
        super.initAndValidate();
        embedding = new int[nodeCount][];
        storedEmbedding = new int[nodeCount][];
        embeddingCount = -1;
        storedEmbeddingCount = -1;
        choicesCount = -1;
        storedChoicesCount = -1;
    }

	public void resetEmbedding(int nCol, int value) {
		final boolean reinitialize = embedding[0] == null || embedding[0].length != nCol;
		for (int r = 0; r < nodeCount; r++) {
			if (reinitialize) embedding[r] = new int[nCol];
			java.util.Arrays.fill(embedding[r], value);
		}
	}

	public void setEmbedding(int row, int col, int value) {
		embedding[row][col] = value;
	}

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
        copyEmbedding(storedEmbedding, etree.storedEmbedding);
        etree.embeddingCount = embeddingCount;
        etree.storedEmbeddingCount = storedEmbeddingCount;
        etree.choicesCount = choicesCount;
        etree.storedChoicesCount = storedChoicesCount;
        return etree;
    }

    @Override
    public void assignTo(final StateNode other) {
    	super.assignTo(other);
        final EmbeddedTree etree = (EmbeddedTree) other;
        copyEmbedding(embedding, etree.embedding);
        copyEmbedding(storedEmbedding, etree.storedEmbedding);
        etree.embeddingCount = embeddingCount;
        etree.storedEmbeddingCount = storedEmbeddingCount;
        etree.choicesCount = choicesCount;
        etree.storedChoicesCount = storedChoicesCount;
    }

    @Override
    public void assignFrom(final StateNode other) {
    	super.assignFrom(other);
        final EmbeddedTree etree = (EmbeddedTree) other;
        copyEmbedding(etree.embedding, embedding);
        copyEmbedding(etree.storedEmbedding, storedEmbedding);
        embeddingCount = etree.embeddingCount;
        storedEmbeddingCount = etree.storedEmbeddingCount;
        choicesCount = etree.choicesCount;
        storedChoicesCount = etree.storedChoicesCount;
    }

    @Override
    public void assignFromFragile(final StateNode other) {
    	super.assignFromFragile(other);
        final EmbeddedTree etree = (EmbeddedTree) other;
        copyEmbedding(etree.embedding, embedding);
        copyEmbedding(etree.storedEmbedding, storedEmbedding);
        embeddingCount = etree.embeddingCount;
        storedEmbeddingCount = etree.storedEmbeddingCount;
        choicesCount = etree.choicesCount;
        storedChoicesCount = etree.storedChoicesCount;
    }

	public void assignFromTree(final StateNode other) {
		super.assignFrom(other);
	}

    public void assignToTree(final StateNode other) {
        super.assignTo(other);
    }

    @Override
    public void store() {
    	super.store();
    	copyEmbedding(embedding, storedEmbedding);
        storedEmbeddingCount = embeddingCount;
        storedChoicesCount = choicesCount;
    }

    @Override
    public void restore() {
    	super.restore();
    	final int[][] tmpEmbedding = embedding;
    	embedding = storedEmbedding;
    	storedEmbedding = tmpEmbedding;
    	int tmpCount = embeddingCount;
        embeddingCount = storedEmbeddingCount;
        storedEmbeddingCount = tmpCount;
        tmpCount = choicesCount;
        choicesCount = storedChoicesCount;
        storedChoicesCount = tmpCount;
    }

    private static void copyEmbedding(int[][] src, int[][] dst) {
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

	public void printEmbedding() {
		System.out.println(String.format("EMBEDDING: %s", getID()));
		for (int r = 0; r < embedding.length; r++) {
			final StringBuilder rowStr = new StringBuilder();
			final int rowLen = embedding[r].length;
			for (int c = 0; c < rowLen - 1; c++) {
				rowStr.append(embedding[r][c]);
				rowStr.append(" ");
			}
			rowStr.append(embedding[r][rowLen - 1]);
			System.out.println(rowStr);
		}
	}

    @Override
    public String toString() {
    	//printEmbedding();
    	for (int i = 0; i < nodeCount; i++) {
    		final int lastTraversalNode = embedding[i].length - 1;
    		final StringBuilder embeddingBuf = new StringBuilder();
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
    	System.out.println(getID() + "\n" + toString());
    	//printEmbedding();
    }
}
