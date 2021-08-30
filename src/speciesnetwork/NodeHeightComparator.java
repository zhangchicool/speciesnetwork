package speciesnetwork;

import java.util.Comparator;

import speciesnetwork.NetworkNode;

// unambiguous, so only identical nodes are considered equal
// if height is equal, consider distance from root (depth)
// if depth is equal, consider assigned node number
public final class NodeHeightComparator implements Comparator<NetworkNode> {
    @Override
    public int compare(NetworkNode nodeA, NetworkNode nodeB) {
        final double heightA = nodeA.getHeight();
        final double heightB = nodeB.getHeight();
        if (heightA == heightB) {
            final int nodeNumberA = nodeA.getNr();
            final int nodeNumberB = nodeB.getNr();
            if (nodeNumberA == nodeNumberB) {
                return 0;
            }
            return nodeNumberA > nodeNumberB ? 1 : -1;
        }
        return heightA > heightB ? 1 : -1;
    }
}
