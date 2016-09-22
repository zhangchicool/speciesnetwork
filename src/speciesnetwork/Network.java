package speciesnetwork;

import java.io.PrintStream;
import java.util.ArrayList;
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
        for (NetworkNode node: nodes) {
            node.updateRelationships();
        }
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

    public int getTraversalNodeCount() {
        return speciationNodeCount + reticulationNodeCount;
    }

    /**
     * @return get the total number of branches in the tree
     */
    public int getBranchCount() {
        return reticulationNodeCount + nodeCount - 1;
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

    public NetworkNode getNode(final int nodeI) {
        return nodes[nodeI];
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

    public void resetAllVisited() {
        for (NetworkNode node: nodes) {
            node.visited = false;
        }
    }

    /**
     * @return (gamma) branch number that corresponds to a node number
     */
    public Integer getBranchNumber(final int nodeNumber) {
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
    public int getNodeNumber(final Integer branchNumber) {
        final int reticulationOffset = getReticulationOffset();
        if (branchNumber < reticulationOffset) {
            return branchNumber;
        } else {
            return (branchNumber - reticulationOffset) / 2 + reticulationOffset;
        }
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
        for (NetworkNode node: nodes) {
            for (NetworkNode parent: node.parents) {
                netLength += parent.height - node.height;
            }
        }
        return netLength;
    }

    public boolean isDirty() {
        for (NetworkNode node: nodes) {
            if (node.isDirty != IS_CLEAN) return true;
        }
        return false;
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

    @Override
    public int scale(final double scale) {
        for (NetworkNode node: nodes) {
            node.height *= scale;
        }
        return speciationNodeCount + reticulationNodeCount;
    }

    @Override
    public String toString() {
        return getOrigin().toString();
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

    private static void copyNetwork(Network src, Network dst) {
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

        for(NetworkNode node: nodes) {
            node.updateRelationships();
            node.isDirty = IS_CLEAN;
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

    /**
     * add a reticulation branch to the species network
     */
    public void addReticulationBranch(NetworkNode reticulationNode, NetworkNode bifurcationNode,
                                      Integer retAttachBranchNr, Integer bifAttachBranchNr) {
        NetworkNode pickedNode1 = getNode(getNodeNumber(retAttachBranchNr));
        NetworkNode pickedNode2 = getNode(getNodeNumber(bifAttachBranchNr));
        NetworkNode parentNode1 = pickedNode1.getParentByBranch(retAttachBranchNr);
        NetworkNode parentNode2 = pickedNode2.getParentByBranch(bifAttachBranchNr);

        // increase the existing reticulation branch numbers by 3
        Integer oldRetBranchNr = retAttachBranchNr;
        Integer oldBifBranchNr = bifAttachBranchNr;
        for (NetworkNode node: getInternalNodes()) {
            List<Integer> newBranchNrs = new ArrayList<>();
            for (Integer bNr: node.childBranchNumbers) {
                if (bNr >= leafNodeCount+speciationNodeCount) {
                    newBranchNrs.add(bNr + 3);
                    if (bNr.equals(oldRetBranchNr)) retAttachBranchNr += 3;
                    if (bNr.equals(oldBifBranchNr)) bifAttachBranchNr += 3;
                } else {
                    newBranchNrs.add(bNr);
                }
            }
            node.childBranchNumbers = newBranchNrs;
        }

        // add the two nodes to the network node array, between the speciation and reticulation nodes
        NetworkNode[] tempNodes = new NetworkNode[nodeCount+2];
        System.arraycopy(nodes, 0, tempNodes, 0, leafNodeCount+speciationNodeCount);
        System.arraycopy(nodes, leafNodeCount+speciationNodeCount, tempNodes, leafNodeCount+speciationNodeCount+2, reticulationNodeCount+1);
        tempNodes[leafNodeCount+speciationNodeCount] = bifurcationNode;
        tempNodes[leafNodeCount+speciationNodeCount+1] = reticulationNode;
        nodes = tempNodes;
        // update the node counts
        nodeCount += 2;
        speciationNodeCount += 1;
        reticulationNodeCount += 1;
        reticulationNode.setLabel("#H" + reticulationNodeCount);
        // bifurcationNode.setLabel("S" + speciationNodeCount);

        if (retAttachBranchNr.equals(bifAttachBranchNr)) {
            // the two nodes are on the same branch
            reticulationNode.childBranchNumbers.add(retAttachBranchNr);
            bifurcationNode.childBranchNumbers.add(leafNodeCount+speciationNodeCount);
            bifurcationNode.childBranchNumbers.add(leafNodeCount+speciationNodeCount+1);
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
        updateRelationships();
    }

    /**
     * delete a reticulation branch from the species network
     * @param hybridBranchNr reticulation branch number
     */
    public void deleteReticulationBranch(Integer hybridBranchNr) {
        // branch hybridBranchNr is connecting hybridNode and pickedParent
        final int hybridNodeNr = getNodeNumber(hybridBranchNr);
        NetworkNode hybridNode = getNode(hybridNodeNr);
        NetworkNode bifurcNode = hybridNode.getParentByBranch(hybridBranchNr);
        final int bifurcNodeNr = bifurcNode.getNr();
        assert (bifurcNode.isSpeciation() && hybridNode.isReticulation());

        // get the child node and another parent node of hybridNode
        final Integer childBranchNr1 = hybridNode.childBranchNumbers.get(0);
        // NetworkNode childNode1 = hybridNode.getChildByBranch(childBranchNr1);
        final Integer parentBranchNr1;
        if (hybridNode.gammaBranchNumber.equals(hybridBranchNr)) {
            parentBranchNr1 = hybridNode.gammaBranchNumber + 1;
        } else {
            parentBranchNr1 = hybridNode.gammaBranchNumber;
        }
        NetworkNode parentNode1 = hybridNode.getParentByBranch(parentBranchNr1);

        // get the parent node and another child node of bifurcNode
        final Integer childBranchNr2;
        if (bifurcNode.childBranchNumbers.get(0).equals(hybridBranchNr)) {
            childBranchNr2 = bifurcNode.childBranchNumbers.get(1);
        } else {
            childBranchNr2 = bifurcNode.childBranchNumbers.get(0);
        }
        NetworkNode childNode2 = bifurcNode.getChildByBranch(childBranchNr2);
        final Integer parentBranchNr2 = bifurcNode.gammaBranchNumber;  // should be equal to bifurcNodeNr
        NetworkNode parentNode2 = bifurcNode.getParentByBranch(parentBranchNr2);

        // update child branch numbers
        if (hybridNode == childNode2 && bifurcNode == parentNode1) {
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
        NetworkNode[] tempNodes = new NetworkNode[nodeCount-2];
        System.arraycopy(nodes, 0, tempNodes, 0, bifurcNodeNr);
        System.arraycopy(nodes, bifurcNodeNr+1, tempNodes, bifurcNodeNr, hybridNodeNr-bifurcNodeNr-1);
        System.arraycopy(nodes, hybridNodeNr+1, tempNodes, hybridNodeNr-1, nodeCount-hybridNodeNr-1);
        nodes = tempNodes;
        // update the node counts
        nodeCount -= 2;
        speciationNodeCount -= 1;
        reticulationNodeCount -= 1;

        // decreasing the child branch numbers between bifurcNodeNr and hybridBranchNr by 1, larger than hybridBranchNr by 3
        for (NetworkNode node: getInternalNodesWithOrigin()) {
            List<Integer> newBranchNrs = new ArrayList<>();
            for (Integer bNr: node.childBranchNumbers) {
                if (bNr > bifurcNodeNr && bNr < hybridBranchNr) {
                    newBranchNrs.add(bNr - 1);
                } else if (bNr > hybridBranchNr) {
                    newBranchNrs.add(bNr - 3);
                } else {
                    newBranchNrs.add(bNr);
                }
            }
            node.childBranchNumbers = newBranchNrs;
        }

        // update the hybrid node labels
        for (int i = 0; i < reticulationNodeCount; i++) {
            nodes[leafNodeCount + speciationNodeCount + i].setLabel("#H" + (i+1));
        }

        // update relationships
        updateRelationships();
    }
}
