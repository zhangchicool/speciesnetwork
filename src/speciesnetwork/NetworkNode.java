package speciesnetwork;

import java.text.DecimalFormat;
import java.util.*;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;

import beast.base.core.Description;
import beast.base.evolution.tree.Node;

/**
 * NetworkNode: tip (one parent no child), bifurcation (one parent two children), reticulation (two parents one child)
 * @author Chi Zhang
 * @author Huw Ogilvie
 */

@Description("Network node for binary rooted network")
public class NetworkNode {
    // the taxonomic name of this node
    protected String label;

    protected int nodeNumber;

    // height of this node
    protected double height;

    // inheritance probability associated with the gamma branch
    protected double inheritProb;

    // children and parents of this node
    public List<Integer> childBranchNumbers;
    public Integer gammaBranchNumber;
    protected Multiset<NetworkNode> children;
    protected Multiset<NetworkNode> parents;

    /**
     * status of this node after an operation is performed on the state
     * A Node IS_DIRTY if its value (like height) has changed.
     * A Node IS_FILTHY if its parent or child has changed.
     * Otherwise the node IS_CLEAN.
     */
    protected int isDirty;

    // used when summarizing posterior distribution
    public Integer subnetworkNr;
    public Double topologySupport;

    // meta-data contained in square brackets in Newick
    protected String metaDataString;

    // arbitrarily labeled metadata on this node
    protected Map<String, Object> metaData = new TreeMap<>();

    // the network that this node is a part of
    protected Network network;

    public NetworkNode() {
        initProps();
    }

    public NetworkNode(Network sNetwork) {
        initProps();
        network = sNetwork;
    }

    private void initProps() {
        label = null;
        height = 0.0;
        inheritProb = 0.5;
        childBranchNumbers = new ArrayList<>();
        children = HashMultiset.create();
        parents = HashMultiset.create();
        nodeNumber = -1;
        isDirty = Network.IS_DIRTY;
    }

    /* instantiate a new network node with the same height, labels and metadata as a tree node
       this does not copy the parents or children */
    public NetworkNode(Node treeNode) {
        label = treeNode.getID();
        height = treeNode.getHeight();
        metaDataString = treeNode.metaDataString;
        for (String key : treeNode.getMetaDataNames()) {
            Object value = treeNode.getMetaData(key);
            metaData.put(key, value);
        }
    }

    protected void copyTo(NetworkNode dst) {
        copyNode(this, dst);
    }

    protected void copyFrom(NetworkNode src) {
        copyNode(src, this);
    }

    protected static void copyNode(NetworkNode src, NetworkNode dst) {
        dst.label = src.label;
        dst.height = src.height;
        dst.inheritProb = src.inheritProb;
        dst.childBranchNumbers.clear();
        dst.childBranchNumbers.addAll(src.childBranchNumbers);
        dst.gammaBranchNumber = src.gammaBranchNumber;
        // children and parents will be sorted out using updateRelationships()
        dst.children.clear();
        dst.parents.clear();
        dst.nodeNumber = src.nodeNumber;
        // also copy meta data?
        dst.metaDataString = src.metaDataString;
        dst.metaData.clear();
        for (String key : src.getMetaDataNames()) {
            Object value = src.getMetaData(key);
            dst.metaData.put(key,value);
        }
        dst.isDirty = src.isDirty;
    }

    public Network getNetwork() {
        return network;
    }

    public int getNr() {
        return nodeNumber;
    }

    public void setNr(int nr) {
        nodeNumber = nr;
    }

    public double getHeight() {
        return height;
    }

    public void setHeight(final double height) {
        this.height = height;
        isDirty |= Network.IS_DIRTY;
        for (NetworkNode child: children) {
            child.isDirty |= Network.IS_DIRTY;
        }
    }

    public void setMetaData(final String pattern, final Object value) {
        metaData.put(pattern, value);
    }

    public Object getMetaData(final String pattern) {
        return metaData.get(pattern);
    }

    public Set<String> getMetaDataNames() {
        return metaData.keySet();
    }

    public Multiset<NetworkNode> getParents() {
        return parents;
    }

    public Multiset<NetworkNode> getChildren() {
        return children;
    }

    protected Multiset<NetworkNode> updateParents() {
        final Multiset<NetworkNode> parents = HashMultiset.create();

        for (NetworkNode node: network.nodes) {
            for (Integer i: node.childBranchNumbers) {
                final int childNodeNumber = network.getNodeNumber(i);
                final NetworkNode childNode = network.nodes[childNodeNumber];
                if (childNode == this) parents.add(node);
            }
        }

        return parents;
    }

    protected Multiset<NetworkNode> updateChildren() {
        final Multiset<NetworkNode> children = HashMultiset.create();

        for (Integer i: childBranchNumbers) {
            final int childNodeNumber = network.getNodeNumber(i);
            final NetworkNode childNode = network.nodes[childNodeNumber];
            children.add(childNode);
        }

        return children;
    }

    public void updateRelationships() {
        nodeNumber = -1;
        for (int i = 0; i < network.nodes.length; i++) {
            if (network.nodes[i] == this) {
                nodeNumber = i;
                break;
            }
        }
        if (nodeNumber < 0) {
            throw new RuntimeException("Node is not attached to the network!");
        }

        gammaBranchNumber = network.getBranchNumber(nodeNumber);
        parents = updateParents();
        children = updateChildren();

        isDirty |= Network.IS_DIRTY;
    }

    public NetworkNode getParentByBranch(Integer branchNr) {
        for (NetworkNode parent: parents) {
            if (parent.childBranchNumbers.contains(branchNr))
                return parent;
        }
        return null;
    }

    public NetworkNode getChildByBranch(Integer branchNr) {
        if (childBranchNumbers.contains(branchNr)) {
            final int childNodeNumber = network.getNodeNumber(branchNr);
            return network.nodes[childNodeNumber];
        }
        return null;
    }

    /**
     * @return true if current node is origin node
     */
    public boolean isOrigin() {
        return parents.size() == 0;
    }

    /**
     * @return true if current node is root node
     */
    public boolean isRoot() {
        return getParentByBranch(gammaBranchNumber).isOrigin();
    }

    /**
     * @return true if current node is leaf node
     */
    public boolean isLeaf() {
        return children.size() == 0;
    }

    /**
     * @return true if current node is reticulation node
     */
    public boolean isReticulation() {
        return parents.size() == 2;
    }

    /**
     * @return true if current node is speciation node
     */
    public boolean isSpeciation() {
        return children.size() == 2;
    }

    /**
     * whether the node has been visited, say by a recursive method
     */
    protected boolean visited = false; // this should be used outside of this class

    /* get and (re)set the visited indicator */
    public boolean isVisited() {
        return visited;
    }

    public void setVisited(boolean v) {
        visited = v;
    }

    private boolean touched = false; // this is only used inside this class

    private void resetAllTouched() {
        for (NetworkNode node: network.nodes) {
            node.touched = false;
        }
    }

    public String toString(DecimalFormat df, boolean inXML) {
        resetAllTouched();
        NetworkNode parent = getParentByBranch(gammaBranchNumber);
        final double parentHeight;
        if (parent == null) {
            parentHeight = Double.POSITIVE_INFINITY;
        } else {
            parentHeight = parent.getHeight();
        }
        return buildNewick(parentHeight, gammaBranchNumber, df, inXML);
    }

    public String toString() {
        return toString(null, false);
    }

    private String buildNewick(double parentHeight, Integer branchNumber, DecimalFormat df, boolean inXML) {
        final StringBuilder subStr = new StringBuilder();
        // only add children to a reticulation node once
        if (children.size() > 0 && !touched) {
            touched = true;
            subStr.append("(");
            int i = 0;
            for (Integer childBranchNr: childBranchNumbers) {
                if (i > 0) subStr.append(",");
                NetworkNode childNode = getChildByBranch(childBranchNr);
                subStr.append(childNode.buildNewick(height, childBranchNr, df, inXML));
                i++;
            }
            subStr.append(")");
        }

        if (label != null)
            subStr.append(label);
        // else
        //  subtreeString.append(nodeNumber);

        if (isOrigin() && topologySupport != null)
            setMetaData("topologySupport", topologySupport);

        if (isReticulation() && gammaBranchNumber.equals(branchNumber)) {
            setMetaData("gamma", inheritProb);
            processMetaData(true);  // write gamma prob associated with the branch
        } else {
            processMetaData(false); // do not write gamma prob at the other branch
        }
        subStr.append(getNewickMetaData(inXML));

        if (parentHeight < Double.POSITIVE_INFINITY) {
            final double branchLength = parentHeight - height;
            subStr.append(":");
            if (df == null) subStr.append(branchLength);
            else subStr.append(df.format(branchLength));
        }

        return subStr.toString();
    }

    /* put meta data in metaDataString TODO: decimal format for meta data */
    private void processMetaData(boolean withGamma) {
        StringBuilder metaStr = new StringBuilder();
        for (String name : getMetaDataNames()) {
            if (!name.contains("gamma") || withGamma) {
                Object value = getMetaData(name);
                metaStr.append(name).append("=");
                if (value instanceof Object[]) {
                    Object[] values = (Object[]) value;
                    metaStr.append("{");
                    for (int i = 0; i < values.length; i++) {
                        if (i > 0) metaStr.append(",");
                        metaStr.append(values[i].toString());
                    }
                    metaStr.append("}");
                } else {
                    metaStr.append(value.toString());
                }
                metaStr.append(",");
            }
        }
        if (metaStr.length() > 0)
            metaDataString = metaStr.substring(0, metaStr.length() - 1);
        else
            metaDataString = "";
    }

    private String getNewickMetaData(boolean inXML) {
        if (metaDataString != null && metaDataString.length() > 0) {
            if (inXML)
                return "[&amp;" + metaDataString + ']';
            else
                return "[&" + metaDataString + ']';
        } else {
            return "";
        }
    }

    public double getGammaProb() {
        return inheritProb;
    }

    public void setGammaProb(final double newGamma) {
        inheritProb = newGamma;
        isDirty |= Network.IS_DIRTY;
    }

    public String getLabel() {
        return label;
    }

    public void setLabel(String newLabel) {
         label = newLabel;
    }

    public int getTraversalNumber() {
        return nodeNumber - network.leafNodeCount;
    }

    /**
     * @return total node count in sub-network defined by this node (including this node itself)
     */
    public int getNodeCount() {
        resetAllTouched();
        return recurseNodeCount();
    }
    private int recurseNodeCount() {
        if (touched) return 0;

        int nodeCount = 1;
        for (NetworkNode child: children) {
            nodeCount += child.recurseNodeCount();
        }

        touched = true;
        return nodeCount;
    }

    public int getLeafNodeCount() {
        resetAllTouched();
        return recurseLeafNodeCount();
    }
    private int recurseLeafNodeCount() {
        if (touched)
            return 0;
        else if (children.size() == 0)
            return 1;

        int nodeCount = 0;
        for (NetworkNode child: children) {
            nodeCount += child.recurseLeafNodeCount();
        }

        touched = true;
        return nodeCount;
    }

    public int getSpeciationNodeCount() {
        resetAllTouched();
        return recurseSpeciationNodeCount();
    }
    private int recurseSpeciationNodeCount() {
        if (touched) return 0;

        // only count speciation nodes
        int nodeCount = (children.size() == 2) ? 1 : 0;
        for (NetworkNode child: children) {
            nodeCount += child.recurseSpeciationNodeCount();
        }

        touched = true;
        return nodeCount;
    }

    public int getReticulationNodeCount() {
        resetAllTouched();
        return recurseReticulationNodeCount();
    }
    private int recurseReticulationNodeCount() {
        if (touched) return 0;

        // only count reticulation nodes
        int nodeCount = (parents.size() == 2) ? 1 : 0;
        for (NetworkNode child: children) {
            nodeCount += child.recurseReticulationNodeCount();
        }

        touched = true;
        return nodeCount;
    }

    public void printDeets() {
        System.out.printf("%s: %s %n", "label", label);
        System.out.printf("%s: %f %n", "inheritProb", inheritProb);
        System.out.printf("%s: %f %n", "height", height);
        System.out.printf("%s: %d %n", "nodeNr", nodeNumber);
        System.out.printf("%s: %d %n", "branchNr", gammaBranchNumber);
        System.out.printf("%s: %d %n", "nParents", parents.size());
        System.out.printf("%s: %d %n", "nChildren", children.size());
        for (Integer i: childBranchNumbers) {
            System.out.printf("%s: %d %n", "childBranchNumber", i);
        }
        System.out.println();
    }
}
