package speciesnetwork;

import beast.evolution.tree.TreeInterface;

/**
 * each row represents the corresponding gene tree node
 * each column represents the corresponding species network traversal node
 */

public interface EmbeddableTree extends TreeInterface {
	public void resetEmbedding(int nCol, int value);

	public void setEmbedding(int row, int col, int value);

	public int getEmbedding(int row, int col);
}
