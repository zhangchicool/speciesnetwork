package speciesnetwork;

import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;

public interface EmbeddedTreeInterface extends TreeInterface {
	Embedding getEmbedding();

	void makeCaterpillar(double rootHeight, double d, boolean b);

	void assignFromTree(Tree tempTree);

	void setEmbedding(Embedding newEmbedding);
}
