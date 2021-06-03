package speciesnetwork;

import beast.core.StateNode;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

public class EmbeddedTree extends Tree {
    public Embedding embedding;
    private Embedding storedEmbedding;

    public EmbeddedTree() {
    }

    // based on Tree(final Node rootNode)
    public EmbeddedTree(final Node rootNode) {
        setRoot(rootNode);
        initArrays();
        embedding = new Embedding(nodeCount);
        storedEmbedding = new Embedding(nodeCount);
    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        embedding = new Embedding(nodeCount);
        storedEmbedding = new Embedding(nodeCount);
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
        etree.embedding = new Embedding(embedding);
        etree.storedEmbedding = new Embedding(storedEmbedding);
        return etree;
    }

    /**
     * copy of all values into existing tree *
     */
    @Override
    public void assignTo(final StateNode other) {
        super.assignTo(other);
        final EmbeddedTree etree = (EmbeddedTree) other;
        etree.embedding = new Embedding(embedding);
        etree.storedEmbedding = new Embedding(storedEmbedding);
    }

    /**
     * copy of all values from existing tree *
     */
    @Override
    public void assignFrom(final StateNode other) {
        super.assignFrom(other);
        final EmbeddedTree etree = (EmbeddedTree) other;
        embedding = new Embedding(etree.embedding);
        storedEmbedding = new Embedding(etree.storedEmbedding);
    }

    /**
     * as assignFrom, but only copy tree structure *
     */
    @Override
    public void assignFromFragile(final StateNode other) {
        super.assignFromFragile(other);
        final EmbeddedTree etree = (EmbeddedTree) other;
        embedding = new Embedding(etree.embedding);
        storedEmbedding = new Embedding(etree.storedEmbedding);
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
        storedEmbedding.copyFrom(embedding);
    }

    @Override
    public void restore() {
        super.restore();
        final Embedding tmpEmbedding = embedding;
        embedding = storedEmbedding;
        storedEmbedding = tmpEmbedding;
    }

    @Override
    public String toString() {
        //printEmbedding();
        for (int i = 0; i < nodeCount; i++) {
            m_nodes[i].metaDataString = "embedding=\"" + embedding.rowToString(i) + "\"";
        }
        return super.toString();
    }

    @Override
    public void fromXML(final org.w3c.dom.Node node) {
        super.fromXML(node);

        // TODO: not working as node numbers might change when restored and embedding becomes invalid
        for (int i = 0; i < nodeCount; i++) {
            final String embedStr = (String) m_nodes[i].getMetaData("embedding");
            String[] parts = embedStr.split(" ");
            if (i == 0) embedding = new Embedding(nodeCount, parts.length);
            for(int j = 0; j < parts.length; j++) {
                final int value = Integer.parseInt(parts[j]);
                embedding.setDirection(i, j, value);
            }
        }
        // System.out.println(getID() + "\n" + toString());
    }
}
