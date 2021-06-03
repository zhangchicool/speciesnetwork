package speciesnetwork;

public class Embedding {
    public int geneNodeCount;
    public int traversalNodeCount;      // number of traversable species network nodes
    protected int[] embedding;          // matrix with nrow=geneNodeCount and ncol=traversalNodeCount
    public double probability = 1.0;    // probability of this embedding
    public double probabilitySum = 1.0; // sum of probabilities of all alternative embeddings

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

    public Embedding(Embedding src) {
        geneNodeCount = src.geneNodeCount;
        traversalNodeCount = src.traversalNodeCount;
        embedding = new int[src.embedding.length];
        System.arraycopy(src.embedding, 0, embedding, 0, embedding.length);
        probability = src.probability;
        probabilitySum = src.probabilitySum;
    }

    public int[] getEmbedding() {
        return embedding;
    }

    public int getDirection(int geneNode, int traversalNode) {
        final int i = (traversalNodeCount * geneNode) + traversalNode;
        return embedding[i];
    }

    public void setDirection(int geneNode, int traversalNode, int value) {
        final int i = (traversalNodeCount * geneNode) + traversalNode;
        embedding[i] = value;
    }

    public void reset(int tnc) {
        // assume that geneNodeCount is not changed and only check traversalNodeCount
        if (traversalNodeCount != tnc) {
            traversalNodeCount = tnc;
            embedding = new int[geneNodeCount * traversalNodeCount];
        }
        java.util.Arrays.fill(embedding, -1);
    }

    public void copyFrom(Embedding src) {
        // assume that geneNodeCount is not changed and only check traversalNodeCount
        if (traversalNodeCount != src.traversalNodeCount) {
            traversalNodeCount = src.traversalNodeCount;
            embedding = new int[src.embedding.length];
        }
        System.arraycopy(src.embedding, 0, embedding, 0, embedding.length);
        probability = src.probability;
        probabilitySum = src.probabilitySum;
    }

    public void mergeWith(Embedding src) {
        assert src.geneNodeCount == geneNodeCount;
        assert src.traversalNodeCount == traversalNodeCount;

        probability *= src.probability;
        probabilitySum *= src.probabilitySum;
        for (int i = 0; i < embedding.length; i++) {
            if (embedding[i] == -1)
                embedding[i] = src.embedding[i];
        }
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
        for (int i = 1; i < embedding.length; i++) {
            str.append(' ');
            str.append(embedding[i]);
        }
        return str.toString();
    }
}
