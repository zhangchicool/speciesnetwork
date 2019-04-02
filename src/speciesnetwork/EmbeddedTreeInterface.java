package speciesnetwork;

import beast.evolution.tree.TreeInterface;

public interface EmbeddedTreeInterface extends TreeInterface {
	Embedding getEmbedding();
}
