package speciesnetwork;

import java.io.PrintStream;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.StateNode;
import beast.evolution.alignment.TaxonSet;
import beast.util.TreeParser;

/**
 * Network class to replace Tree class
 * It includes tip (in-degree 1/out-degree 0), bifurcation (1/2), and reticulation (2/1) nodes.
 * @author Chi Zhang
 * @author Huw Ogilvie
 */

@Description("Network representing reticulate evolution of species")
public class Network extends StateNode {
    final public Input<TaxonSet> taxonSetInput =
            new Input<>("taxonset", "Set of taxa at the leafs of the network.");

    /**
     * state of dirtiness of a node in the network
     * DIRTY means a property on the node has changed, but not the topology. "property" includes the node height
     *       and that branch length to its parent.
     * FILTHY means the node's parent or child has changed.
     */
    public static final int IS_CLEAN = 0, IS_DIRTY = 1, IS_FILTHY = 2;

    // number of nodes
    protected int nodeCount = -1;
    private int storedNodeCount = -1;
    protected int speciationNodeCount = -1;
    private int storedSpeciationNodeCount = -1;
    protected int leafNodeCount = -1;
    private int storedLeafNodeCount = -1;
    protected int reticulationNodeCount = -1;
    private int storedReticulationNodeCount = -1;

    /**
     * array of all nodes in the network
     */
    protected NetworkNode[] nodes = null;
    private NetworkNode[] storedNodes = null;

    @Override
    public void initAndValidate() {
        if (nodeCount < 0) {
            if (taxonSetInput.get() != null) {
                makeCaterpillar(0, 1);
                updateRelationships();
            } else {
                // make dummy network with a single root node
                nodes = new NetworkNode[1];
                nodes[1] = new NetworkNode(this);
                nodeCount = leafNodeCount = 1;
                speciationNodeCount = reticulationNodeCount = 0;
            }
        }
    }

    public void updateRelationships() {
        for (NetworkNode n: nodes) {
            n.updateRelationships();
        }
    }

    @Override
    public int scale(final double scale) {
        for (NetworkNode n: nodes) {
            n.height = n.height * scale;
        }
        return speciationNodeCount + reticulationNodeCount;
    }

    private void makeCaterpillar(final double minInternalHeight, final double step) {
        // make a caterpillar
        final List<String> taxa = taxonSetInput.get().asStringList();
        leafNodeCount = taxa.size();
        speciationNodeCount = leafNodeCount - 1;
        reticulationNodeCount = 0;
        nodeCount = leafNodeCount * 2;
        nodes = new NetworkNode[nodeCount];

        int leftNr = 0;
        nodes[leftNr] = new NetworkNode(this);
        NetworkNode left = nodes[leftNr];
        left.height = 0.0;
        left.label = taxa.get(leftNr);
        for (int rightNr = 1; rightNr < leafNodeCount; rightNr++) {
            nodes[rightNr] = new NetworkNode(this);
            final NetworkNode right = nodes[rightNr];
            right.height = 0.0;
            right.label = taxa.get(rightNr);
            final int parentNr = leafNodeCount + (rightNr - 1);
            nodes[parentNr] = new NetworkNode(this);
            final NetworkNode parent = nodes[parentNr];
            parent.height = minInternalHeight + rightNr * step;
            parent.childBranchNumbers.add(rightNr);
            parent.childBranchNumbers.add(leftNr);
            // left = parent;
            leftNr = parentNr;
        }

        // node of origin
        nodes[nodeCount - 1] = new NetworkNode(this);
        final NetworkNode origin = nodes[nodeCount - 1];
        origin.childBranchNumbers.add(leftNr);
    }

    public NetworkNode getOrigin() {
        return nodes[nodeCount - 1];
    }

    public NetworkNode getRoot() {
        NetworkNode origin = getOrigin();
        int childBrNr = origin.childBranchNumbers.get(0);
        return origin.getChildByBranch(childBrNr);
    }

    public void swapNodes(final int nodeI, final int nodeJ) {
        final NetworkNode tmp = nodes[nodeI];
        nodes[nodeI] = nodes[nodeJ];
        nodes[nodeJ] = tmp;
    }

    /**
     * @return the number of nodes
     */
    public int getNodeCount() {
        return nodeCount;
    }

    public int getSpeciationNodeCount() {
        return speciationNodeCount;
    }

    public int getLeafNodeCount() {
        return leafNodeCount;
    }

    public int getReticulationNodeCount() {
        return reticulationNodeCount;
    }

    /**
     * @return the index of the first reticulation node
     */
    public int getReticulationOffset() {
        return leafNodeCount + speciationNodeCount;
    }

    /**
     * @return get the total number of branches in the tree
     */
    public int getBranchCount() {
        return nodeCount + reticulationNodeCount - 1;
    }

    /**
     * @return the number of branches at the given time
     */
    public int getBranchCount(double time) {
        int nB = 0;
        for (NetworkNode node: nodes) {
            for (NetworkNode child: node.getChildren()) {
                if (node.getHeight() > time && child.getHeight() <= time)
                    nB++;
            }
        }
        return nB;
    }

    public int getTraversalNodeCount() {
        return speciationNodeCount + reticulationNodeCount;
    }

    /**
     * @return an array of all the nodes in this network
     */
    public NetworkNode[] getAllNodes() {
        // Q2HO: why not return nodes directly?
        final NetworkNode[] nodesCopy = new NetworkNode[nodeCount];
        System.arraycopy(nodes, 0, nodesCopy, 0, nodeCount);
        return nodesCopy;
    }

    public NetworkNode[] getAllNodesExceptOrigin() {
        final NetworkNode[] nodesCopy = new NetworkNode[nodeCount - 1];
        System.arraycopy(nodes, 0, nodesCopy, 0, nodeCount - 1);
        return nodesCopy;
    }

    /**
     * @return an array of leaf nodes in this network
     */
    public NetworkNode[] getLeafNodes() {
        final NetworkNode[] leafNodes = new NetworkNode[leafNodeCount];
        System.arraycopy(nodes, 0, leafNodes, 0, leafNodeCount);
        return leafNodes;
    }

    /**
     * speciation and reticulation nodes, do not include origin node
     * @return an array of internal nodes in this network
     */
    public NetworkNode[] getInternalNodes() {
        final int internalNodeCount = speciationNodeCount + reticulationNodeCount;
        final NetworkNode[] internalNodes = new NetworkNode[internalNodeCount];
        System.arraycopy(nodes, leafNodeCount, internalNodes, 0, internalNodeCount);
        return internalNodes;
    }

    public NetworkNode[] getInternalNodesWithOrigin() {
        final int internalNodeCount = speciationNodeCount + reticulationNodeCount;
        final NetworkNode[] internalNodes = new NetworkNode[internalNodeCount + 1];
        System.arraycopy(nodes, leafNodeCount, internalNodes, 0, internalNodeCount + 1);
        return internalNodes;
    }

    /**
     * @return an array of reticulation nodes in this network
     */
    public NetworkNode[] getReticulationNodes() {
        final int reticulationOffset = getReticulationOffset();
        final NetworkNode[] reticulationNodes = new NetworkNode[reticulationNodeCount];
        System.arraycopy(nodes, reticulationOffset, reticulationNodes, 0, reticulationNodeCount);
        return reticulationNodes;
    }

    public NetworkNode getNode(final int nodeI) {
        return nodes[nodeI];
    }

    public int getNodeNumber(final String query) {
        for (int i = 0; i < nodeCount; i++) {
            final NetworkNode n = nodes[i];
            if (n != null && n.label != null) {
                if (n.label.equals(query)) return i;
            }
        }
        return -1; // no match
    }

    public double getNetworkLength() {
        double netLength = 0;
        for (NetworkNode n: nodes) {
            for (NetworkNode p: n.parents) {
                netLength += p.height - n.height;
            }
        }
        return netLength;
    }

    @Override
    public String toString() {
        return getOrigin().toString();
    }

    @Override
    public void setEverythingDirty(final boolean isDirty) {
        setSomethingIsDirty(isDirty);
        if (!isDirty) {
            for(NetworkNode node : nodes) {
                node.isDirty = IS_CLEAN;
            }
        } else {
            for(NetworkNode node : nodes) {
                node.isDirty = IS_FILTHY;
            }
        }
    }

    /**
     * @return a deep copy of this network
     */
    @Override
    public Network copy() {
        Network copy = new Network();
        copy.index = index;
        copyNetwork(this, copy);
        return copy;
    }
    
    /**
     * copy of all values into existing network
     */
    @Override
    public void assignTo(final StateNode other) {
        final Network dst = (Network) other;
        copyNetwork(this, dst);
    }

    /**
     * copy of all values from existing network
     */
    @Override
    public void assignFrom(final StateNode other) {
        final Network src = (Network) other;
        copyNetwork(src, this);
    }

    protected static void copyNetwork(Network src, Network dst) {
        final int copyNodeCount = src.nodeCount;

        dst.setID(src.getID());
        dst.nodeCount = copyNodeCount;
        dst.speciationNodeCount = src.speciationNodeCount;
        dst.leafNodeCount = src.leafNodeCount;
        dst.reticulationNodeCount = src.reticulationNodeCount;

        dst.nodes = new NetworkNode[copyNodeCount];
        for (int i = 0; i < copyNodeCount; i++) {
            dst.nodes[i] = new NetworkNode(dst);
            dst.nodes[i].copyFrom(src.nodes[i]);
        }
        dst.updateRelationships();
    }

    /**
     * as assignFrom, but assumes this network has been initialized
     * with the same dimensions as the source network
     */
    @Override
    public void assignFromFragile(final StateNode other) {
        final Network src = (Network) other;

        for (int i = 0; i < nodeCount; i++) {
            nodes[i].copyFrom(src.nodes[i]);
        }
        updateRelationships();
    }

    /**
     * reconstruct tree from XML fragment in the form of a DOM node
     */
    @Override
    public void fromXML(final org.w3c.dom.Node node) {
        final String newick = node.getTextContent();
        final TreeParser parser = new TreeParser();
        try {
            parser.thresholdInput.setValue(1e-10, parser);
        } catch (Exception e1) {
            e1.printStackTrace();
        }
        try {
            parser.offsetInput.setValue(0, parser);
            parser.parseNewick(newick); // TODO covert to network
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Override
    protected void store() {
        // System.out.println(String.format("Storing network state... %d = %d + %d + %d", nodeCount, leafNodeCount, speciationNodeCount, reticulationNodeCount));
        storedNodeCount = nodeCount;
        storedSpeciationNodeCount = speciationNodeCount;
        storedLeafNodeCount = leafNodeCount;
        storedReticulationNodeCount = reticulationNodeCount;
        storedNodes = new NetworkNode[nodeCount];

        for (int i = 0; i < nodeCount; i++) {
            storedNodes[i] = new NetworkNode(this);
            storedNodes[i].copyFrom(nodes[i]);
        }
    }

    @Override
    public void restore() {
        // System.out.println("Restoring network state...");
        int tmpNodeCount = nodeCount;
        nodeCount = storedNodeCount;
        storedNodeCount = tmpNodeCount;

        int tmpSpeciationNodeCount = speciationNodeCount;
        speciationNodeCount = storedSpeciationNodeCount;
        storedSpeciationNodeCount = tmpSpeciationNodeCount;

        int tmpLeafNodeCount = leafNodeCount;
        leafNodeCount = storedLeafNodeCount;
        storedLeafNodeCount = tmpLeafNodeCount;

        int tmpReticulationNodeCount = reticulationNodeCount;
        reticulationNodeCount = storedReticulationNodeCount;
        storedReticulationNodeCount = tmpReticulationNodeCount;

        NetworkNode[] tmpNodes = nodes;
        nodes = storedNodes;
        storedNodes = tmpNodes;

        hasStartedEditing = false;

        for(NetworkNode n: nodes) {
            n.updateRelationships();
            n.isDirty = IS_CLEAN;
        }
    }

    /** Loggable interface implementation follows **/

    @Override
    public void init(PrintStream out) {
        out.println("#NEXUS\n");
        out.println("Begin taxa;");
        out.println("\tDimensions ntax=" + leafNodeCount + ";");
        out.println("\t\tTaxlabels");
        for (int i = 0; i < leafNodeCount; i++) {
            out.println("\t\t\t" + nodes[i].label);
        }
        out.println("\t\t\t;");
        out.println("End;");
        out.println("Begin trees;");
    }

    @Override
    public void log(int sample, PrintStream out) {
        Network network = (Network) getCurrent();
        out.print("tree STATE_" + sample + " = ");
        final String newick = network.toString();
        out.print(newick);
        out.print(";");
    }

    /**
     * @see beast.core.Loggable *
     */
    @Override
    public void close(PrintStream out) {
        out.print("End;");
    }

    @Override
    public int getDimension() {
        return nodeCount;
    }

    @Override
    public double getArrayValue() {
        return getRoot().height;
    }

    @Override
    public double getArrayValue(final int nodeI) {
        return nodes[nodeI].height;
    }

    public boolean isDirty() {
        for (NetworkNode n: nodes) {
            if (n.isDirty != IS_CLEAN) return true;
        }
        return false;
    }

    public void resetAllVisited() {
        for (NetworkNode n: nodes) {
            n.visited = false;
        }
    }

    /**
     * @return (gamma) branch number that corresponds to a node number
     */
    public int getBranchNumber(final int nodeNumber) {
        final int reticulationOffset = getReticulationOffset();
        if (nodeNumber < reticulationOffset) {
            return nodeNumber;
        } else {
            return (nodeNumber * 2) - reticulationOffset;
        }
    }

    /**
     * @return node number that corresponds to a branch number
     */
    public int getNodeNumber(final int branchNumber) {
        final int reticulationOffset = getReticulationOffset();
        if (branchNumber < reticulationOffset) {
            return branchNumber;
        } else {
            return (branchNumber - reticulationOffset) / 2 + reticulationOffset;
        }
    }

    /**
     * add a reticulation branch to the species network
     */
    public void addReticulationBranch(NetworkNode reticulationNode, NetworkNode bifurcationNode,
                                      int retAttachBranchNr, int bifAttachBranchNr) {

        final int pickedNodeNr1 = getNodeNumber(retAttachBranchNr);
        NetworkNode pickedNode1 = getNode(pickedNodeNr1);
        final int pickedNodeNr2 = getNodeNumber(bifAttachBranchNr);
        NetworkNode pickedNode2 = getNode(pickedNodeNr2);
        NetworkNode parentNode1 = pickedNode1.getParentByBranch(retAttachBranchNr);
        NetworkNode parentNode2 = pickedNode2.getParentByBranch(bifAttachBranchNr);

        // add the two nodes to the network node array
        NetworkNode[] tempNodes = new NetworkNode[nodeCount+2];
        System.arraycopy(nodes, 0, tempNodes, 0, leafNodeCount+speciationNodeCount);
        System.arraycopy(nodes, leafNodeCount+speciationNodeCount, tempNodes, leafNodeCount+speciationNodeCount+2, reticulationNodeCount+1);
        tempNodes[leafNodeCount+speciationNodeCount] = bifurcationNode;
        tempNodes[leafNodeCount+speciationNodeCount+1] = reticulationNode;
        nodes = tempNodes;
        nodeCount += 2;
        speciationNodeCount += 1;
        reticulationNodeCount += 1;

        // update child branch numbers
        for (NetworkNode node: nodes) {
            for (Integer nr: node.childBranchNumbers) {
                if (nr >= leafNodeCount+speciationNodeCount-1) {
                    node.childBranchNumbers.remove(nr);
                    node.childBranchNumbers.add(nr+3);
                    if (retAttachBranchNr == nr) retAttachBranchNr += 3;
                    if (bifAttachBranchNr == nr) bifAttachBranchNr += 3;
                }
            }
        }
        if (pickedNode1 == pickedNode2 && parentNode1 == parentNode2) {
            // the two nodes are on the same branch
            bifurcationNode.childBranchNumbers.add(leafNodeCount+speciationNodeCount);
            bifurcationNode.childBranchNumbers.add(leafNodeCount+speciationNodeCount+1);
            reticulationNode.childBranchNumbers.add(retAttachBranchNr);
            parentNode2.childBranchNumbers.remove(bifAttachBranchNr);
            parentNode2.childBranchNumbers.add(leafNodeCount+speciationNodeCount-1);
        } else {
            reticulationNode.childBranchNumbers.add(retAttachBranchNr);
            parentNode1.childBranchNumbers.remove(retAttachBranchNr);
            parentNode1.childBranchNumbers.add(leafNodeCount+speciationNodeCount);
            bifurcationNode.childBranchNumbers.add(leafNodeCount+speciationNodeCount+1);
            bifurcationNode.childBranchNumbers.add(bifAttachBranchNr);
            parentNode2.childBranchNumbers.remove(bifAttachBranchNr);
            parentNode2.childBranchNumbers.add(leafNodeCount+speciationNodeCount-1);
        }

        // update relationships
        for (NetworkNode node: nodes) {
            node.updateRelationships();
        }
    }

    /**
     * delete a reticulation branch from the species network
     * @param hybridBranchNr reticulation branch number
     */
    public void deleteReticulationBranch(int hybridBranchNr) {
        final int hybridNodeNr = getNodeNumber(hybridBranchNr);
        NetworkNode hybridNode = getNode(hybridNodeNr);

        // branch hybridBranchNr is connecting hybridNode and pickedParent
        NetworkNode pickedParent = hybridNode.getParentByBranch(hybridBranchNr);
        final int pickedParentNr = pickedParent.getNr();

        // get the child node and another parent node of hybridNode
        final int childBranchNr1 = hybridNode.childBranchNumbers.get(0);
        NetworkNode childNode1 = hybridNode.getChildByBranch(childBranchNr1);
        final int parentBranchNr1;
        if (hybridBranchNr == hybridNode.gammaBranchNumber) {
            parentBranchNr1 = hybridBranchNr + 1;
        } else {
            parentBranchNr1 = hybridBranchNr;
        }
        NetworkNode parentNode1 = hybridNode.getParentByBranch(parentBranchNr1);

        // get the parent node and another child node of pickedParent
        final int childBranchNr2;
        if (hybridBranchNr == pickedParent.childBranchNumbers.get(0)) {
            childBranchNr2 = pickedParent.childBranchNumbers.get(1);
        } else {
            childBranchNr2 = pickedParent.childBranchNumbers.get(0);
        }
        NetworkNode childNode2 = pickedParent.getChildByBranch(childBranchNr2);
        final int parentBranchNr2 = pickedParent.gammaBranchNumber;
        NetworkNode parentNode2 = pickedParent.getParentByBranch(parentBranchNr2);

        // update children parents relationship
        NetworkNode[] tempNodes = new NetworkNode[nodeCount-2];

        if (hybridNode == childNode2 && pickedParent == parentNode1) {
            // the two nodes are on the same branch
            parentNode2.childBranchNumbers.remove(parentBranchNr2);
            parentNode2.childBranchNumbers.add(childBranchNr1);
        } else {
            parentNode1.childBranchNumbers.remove(parentBranchNr1);
            parentNode1.childBranchNumbers.add(childBranchNr1);
            parentNode2.childBranchNumbers.remove(parentBranchNr2);
            parentNode2.childBranchNumbers.add(childBranchNr2);
        }

        // remove the two nodes from the network
        System.arraycopy(nodes, 0, tempNodes, 0, pickedParentNr);
        System.arraycopy(nodes, pickedParentNr+1, tempNodes, pickedParentNr, hybridNodeNr-pickedParentNr-1);
        System.arraycopy(nodes, hybridNodeNr+1, tempNodes, hybridNodeNr-1, nodeCount-hybridNodeNr-1);
        nodes = tempNodes;
        nodeCount -= 2;
        speciationNodeCount -= 1;
        reticulationNodeCount -= 1;

        // update child branch numbers
        for (NetworkNode node: nodes) {
            for (Integer nr: node.childBranchNumbers) {
                if (nr > pickedParentNr && nr < hybridNodeNr) {
                    node.childBranchNumbers.remove(nr);
                    node.childBranchNumbers.add(nr-1);
                } else if (nr > hybridNodeNr) {
                    node.childBranchNumbers.remove(nr);
                    node.childBranchNumbers.add(nr-3);
                }
            }
        }

        // update relationships
        for (NetworkNode node: nodes) {
            node.updateRelationships();
        }
    }
}
