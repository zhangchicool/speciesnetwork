package speciesnetwork;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.TreeParser;

/**
 * Network class modelling the phylogeny of species, which includes
 * tip (in-degree 1/out-degree 0), bifurcation (1/2), and reticulation (2/1) nodes.
 * @author Chi Zhang
 * @author Huw Ogilvie
 */

@Citation("Zhang C., Ogilvie H.A., Drummond A.J., Stadler T. (2018).\n" +
        "  Bayesian inference of species networks from multilocus sequence data.\n" +
        "  Molecular Biology and Evolution 35(2):504–517.")

@Description("Network representing reticulate evolution of species")
public class Network extends StateNode {
    public final Input<TaxonSet> taxonSetInput =
            new Input<>("taxonset", "Set of taxa at the leafs of the network.");
    public final Input<List<TraitSet>> traitListInput = new Input<>("trait",
            "trait information for initializing traits (like node dates) in the network", new ArrayList<>());

    /**
     * state of dirtiness of a node in the network
     * DIRTY means a property on the node has changed, but not the topology. "property" includes the node height
     *       and that branch length to its parent.
     * FILTHY means the node's parent or child has changed.
     */
    public static final int IS_CLEAN = 0, IS_DIRTY = 1, IS_FILTHY = 2;

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
     * the order must obey: leaf nodes | speciation nodes | reticulation nodes | origin node
     */
    protected NetworkNode[] nodes = null;
    private NetworkNode[] storedNodes = null;
    private static Comparator<NetworkNode> hc = new NodeHeightComparator();

    // trait set which specifies leaf node times
    protected TraitSet timeTraitSet = null;

    @Override
    public void initAndValidate() {
        if (nodeCount < 0) {
            makeCaterpillar(0, 1);
        }

        // process all available time traits
        processTraits(traitListInput.get());

        // ensure network is compatible with time traits
        if (timeTraitSet != null) {
            adjustNodeHeights(getOrigin());
        }
    }

    public void makeDummy() {
        // make dummy network with a single node
        nodes = new NetworkNode[1];
        nodes[0] = new NetworkNode(this);
        nodeCount = 1;
        leafNodeCount = speciationNodeCount = reticulationNodeCount = 0;
    }

    private void makeCaterpillar(final double minInternalHeight, final double step) {
        // make a caterpillar species tree
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
        nodes[nodeCount - 1].height = minInternalHeight + leafNodeCount * step;
        nodes[nodeCount - 1].childBranchNumbers.add(leftNr);

        // set internal node labels
        resetInternalNodeLabels();

        updateRelationships();
    }

    public void updateRelationships() {
        for (NetworkNode node: nodes) {
            node.updateRelationships();
        }
    }

    /**
     * @return date trait set if available, null otherwise.
     */
    public TraitSet getDateTrait() {
        return timeTraitSet;
    }

    /**
     * @return true if network has a date/time trait set associated
     */
    public boolean hasDateTrait() {
        return timeTraitSet != null;
    }

    /* set the date trait set for this network.
       A null value simply removes the existing trait set. */
    public void setDateTrait(TraitSet traitSet) {
        if (hasDateTrait()) {
            traitListInput.get().remove(timeTraitSet);
        }

        if (traitSet != null)
            traitListInput.get().add(traitSet);

        timeTraitSet = traitSet;
    }

    protected void processTraits(List<TraitSet> traitList) {
        for (TraitSet traitSet : traitList) {
            for (NetworkNode node : getLeafNodes()) {
                final String speciesName = node.getLabel();
                if (speciesName != null) {
                    final double value = traitSet.getValue(speciesName);
                    // add trait to node meta data
                    node.setMetaData(traitSet.getTraitName(), value);

                    // also update tip node height
                    if (traitSet.isDateTrait()) {
                       node.setHeight(value);
                    }
                }
            }
            if (traitSet.isDateTrait())
                timeTraitSet = traitSet;
        }
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

    public int getInternalNodeCount() {
        // internal nodes include bifurcation and reticulation nodes, do not include origin node
        return speciationNodeCount + reticulationNodeCount;
    }

    /**
     * @return the index of the first reticulation node
     */
    public int getReticulationOffset() {
        return leafNodeCount + speciationNodeCount;
    }

    /**
     * @return the total number of branches in the network
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

    public NetworkNode getNode(final int idx) {
        return nodes[idx];
    }

    /**
     * @return an array of all the nodes in this network
     */
    public NetworkNode[] getAllNodes() {
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
     * speciation and reticulation nodes (do not include origin node)
     * @return an array of internal nodes in this network
     */
    public NetworkNode[] getInternalNodes() {
        final int internalNodeCount = speciationNodeCount + reticulationNodeCount;
        final NetworkNode[] internalNodes = new NetworkNode[internalNodeCount];
        System.arraycopy(nodes, leafNodeCount, internalNodes, 0, internalNodeCount);
        return internalNodes;
    }

    /**
     * speciation and reticulation nodes plus origin node
     */
    public NetworkNode[] getInternalNodesWithOrigin() {
        final int internalNodeCount = speciationNodeCount + reticulationNodeCount;
        final NetworkNode[] internalNodes = new NetworkNode[internalNodeCount + 1];
        System.arraycopy(nodes, leafNodeCount, internalNodes, 0, internalNodeCount + 1);
        return internalNodes;
    }

    /**
     * @return an array of speciation nodes in this network
     */
    public NetworkNode[] getSpeciationNodes() {
        final NetworkNode[] speciationNodes = new NetworkNode[speciationNodeCount];
        System.arraycopy(nodes, leafNodeCount, speciationNodes, 0, speciationNodeCount);
        return speciationNodes;
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

    /* set visited indicator to false for all nodes
       This is typically called before a recursive function to avoid duplicated traversal in the network */
    public void resetAllVisited() {
        for (NetworkNode node: nodes) {
            node.setVisited(false);
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

    /* adjust node heights to no negative branch lengths exist
       This can occur if leaf heights given as traits are incompatible with the existing network. */
    public void adjustNodeHeights(NetworkNode node) {
        if (!node.isLeaf()) {
            for (final NetworkNode child : node.getChildren()) {
                adjustNodeHeights(child);
            }
            for (final NetworkNode child : node.getChildren()) {
                // bump up a small amount
                final double minHeight = child.height + 0.001;
                if (node.height < minHeight)
                    node.height = minHeight;
            }
        }
    }

    @Override
    public int scale(final double scale) {
        int dof = 0;
        for (NetworkNode node : getInternalNodesWithOrigin()) {
            node.height *= scale;
            dof++;

            // check for negative branch length
            for (NetworkNode child : node.getChildren()) {
                if (child.isLeaf() && child.height > node.height)
                    throw new RuntimeException("Scale gives negative branch length!");
            }
        }
        return dof;
    }

    @Override
    public String toString() {
        return getOrigin().toString();
    }

    public String toString(DecimalFormat df) {
        return getOrigin().toString(df, false);
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
        dst.setID(src.getID());
        dst.nodeCount = src.nodeCount;
        dst.speciationNodeCount = src.speciationNodeCount;
        dst.leafNodeCount = src.leafNodeCount;
        dst.reticulationNodeCount = src.reticulationNodeCount;

        dst.nodes = new NetworkNode[src.nodeCount];
        dst.storedNodes = new NetworkNode[src.nodeCount];
        for (int i = 0; i < src.nodeCount; i++) {
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

        // assumes src.nodeCount == nodeCount
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
        final TreeParser treeParser = new TreeParser(newick);
        final NetworkParser networkParser = new NetworkParser(treeParser);
        assignFrom(networkParser);
        // networkParser.setID("dummy");
    }

    @Override
    protected void store() {
        // System.out.println("Storing network state...");
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
    public void log(long sample, PrintStream out) {
        Network network = (Network) getCurrent();
        out.print("tree STATE_" + sample + " = ");
        final String newick = network.toString();
        out.print(newick);
        out.print(";");
    }

    /**
     * @see beast.base.core.Loggable *
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
    public double getArrayValue(final int idx) {
        return nodes[idx].height;
    }

    public void resetInternalNodeLabels() {
        // reset the speciation and reticulation node labels
        for (int i = 0; i < speciationNodeCount; i++) {
            nodes[leafNodeCount + i].setLabel("S" + (i+1));
        }
        for (int i = 0; i < reticulationNodeCount; i++) {
            nodes[leafNodeCount + speciationNodeCount + i].setLabel("#H" + (i+1));
        }
    }

    /**
     * @return the backbone tree of this network
     */
    public Network getBackboneTree () {
        // make a copy of the current network to modify
        Network backbone = new Network();
        backbone.assignFrom(this);

        // get and sort the hybridization nodes
        List<NetworkNode> hNodes = Arrays.asList(backbone.getReticulationNodes());
        hNodes.sort(hc);

        // start deleting from the oldest to the youngest
        // this guarantees every hybridization branch can be safely deleted
        for (int i = hNodes.size() - 1; i >= 0; i--) {
            final NetworkNode node = hNodes.get(i);
            if (node.getGammaProb() < 0.5) {
                // delete the gamma branch
                backbone.deleteReticulationBranch(node.gammaBranchNumber);
            }
            else {
                // delete the alternative branch
                backbone.deleteReticulationBranch(node.gammaBranchNumber + 1);
            }
        }

        return backbone;
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
        bifurcationNode.setLabel("S" + speciationNodeCount);

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
     * @param reticuBranchNr reticulation branch number
     */
    public void deleteReticulationBranch(Integer reticuBranchNr) {
        // branch with reticuBranchNr is connecting hybridNode and bifurcNode
        final int hybridNodeNr = getNodeNumber(reticuBranchNr);
        NetworkNode hybridNode = getNode(hybridNodeNr);
        NetworkNode bifurcNode = hybridNode.getParentByBranch(reticuBranchNr);
        final int bifurcNodeNr = bifurcNode.getNr();
        assert (bifurcNode.isSpeciation() && hybridNode.isReticulation());

        // get the parent node and another child node of bifurcNode
        final Integer pNParentBranchNr = bifurcNode.gammaBranchNumber;  // should be equal to bifurcNodeNr
        NetworkNode pNParentNode = bifurcNode.getParentByBranch(pNParentBranchNr);
        final Integer pNChildBranchNr;
        if (bifurcNode.childBranchNumbers.get(0).equals(reticuBranchNr))
            pNChildBranchNr = bifurcNode.childBranchNumbers.get(1);
        else
            pNChildBranchNr = bifurcNode.childBranchNumbers.get(0);
        NetworkNode pNChildNode = bifurcNode.getChildByBranch(pNChildBranchNr);

        // get the child node and another parent node of hybridNode
        final Integer hNChildBranchNr = hybridNode.childBranchNumbers.get(0);
        // NetworkNode hNChildNode = hybridNode.getChildByBranch(hNChildBranchNr);
        final Integer hNParentBranchNr;
        if (hybridNode.gammaBranchNumber.equals(reticuBranchNr))
            hNParentBranchNr = hybridNode.gammaBranchNumber + 1;
        else
            hNParentBranchNr = hybridNode.gammaBranchNumber;
        NetworkNode hNParentNode = hybridNode.getParentByBranch(hNParentBranchNr);

        // update child branch numbers
        if (bifurcNode == hNParentNode && hybridNode == pNChildNode) {
            // the two nodes are on the same branch
            pNParentNode.childBranchNumbers.remove(pNParentBranchNr);
            pNParentNode.childBranchNumbers.add(hNChildBranchNr);
        } else {
            pNParentNode.childBranchNumbers.remove(pNParentBranchNr);
            pNParentNode.childBranchNumbers.add(pNChildBranchNr);
            hNParentNode.childBranchNumbers.remove(hNParentBranchNr);
            hNParentNode.childBranchNumbers.add(hNChildBranchNr);
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

        for (NetworkNode node: getInternalNodesWithOrigin()) {
            List<Integer> newBranchNrs = new ArrayList<>();
            for (Integer bNr: node.childBranchNumbers) {
                if (bNr > bifurcNodeNr && bNr < reticuBranchNr)
                    newBranchNrs.add(bNr - 1);  // decreasing the child branch numbers by 1
                else if (bNr > reticuBranchNr)
                    newBranchNrs.add(bNr - 3);  // decreasing the child branch numbers by 3
                else
                    newBranchNrs.add(bNr);
            }
            node.childBranchNumbers = newBranchNrs;
        }

        // update the speciation and reticulation node labels
        resetInternalNodeLabels();
        // update relationships
        updateRelationships();
    }

    /* add a speciation node to the nodes array */
    public void addSpeciationNode(NetworkNode sNode) {
        NetworkNode[] tmpNodes = new NetworkNode[nodeCount + 1];
        System.arraycopy(nodes, 0, tmpNodes, 0, leafNodeCount+speciationNodeCount);
        System.arraycopy(nodes, leafNodeCount+speciationNodeCount, tmpNodes, leafNodeCount+speciationNodeCount+1, reticulationNodeCount+1);
        // add the node to the end of the speciation nodes and before the reticulation nodes
        tmpNodes[leafNodeCount+speciationNodeCount] = sNode;
        nodes = tmpNodes;
        nodeCount++;
        speciationNodeCount++;
    }

    /* add a reticulation node to the nodes array */
    public void addReticulationNode(NetworkNode rNode) {
        NetworkNode[] tmpNodes = new NetworkNode[nodeCount + 1];
        System.arraycopy(nodes, 0, tmpNodes, 0, leafNodeCount+speciationNodeCount);
        System.arraycopy(nodes, leafNodeCount+speciationNodeCount, tmpNodes, leafNodeCount+speciationNodeCount+1, reticulationNodeCount+1);
        // add the node to the end of the speciation nodes and before the reticulation nodes
        tmpNodes[leafNodeCount+speciationNodeCount] = rNode;
        nodes = tmpNodes;
        nodeCount++;
        reticulationNodeCount++;
    }

    /* add a leaf node to the nodes array */
    public void addLeafNode(NetworkNode lNode) {
        NetworkNode[] tmpNodes = new NetworkNode[nodeCount + 1];
        System.arraycopy(nodes, 0, tmpNodes, 0, leafNodeCount);
        System.arraycopy(nodes, leafNodeCount, tmpNodes, leafNodeCount+1, speciationNodeCount+reticulationNodeCount+1);
        // add the node to the end of the leaf nodes and before the speciation nodes
        tmpNodes[leafNodeCount] = lNode;
        nodes = tmpNodes;
        nodeCount++;
        leafNodeCount++;
    }

    /* delete a node from the nodes array */
    public void deleteNode(NetworkNode node) {
        int index = -1;
        for (int i = 0; i < nodes.length; i++) {
            if (nodes[i] == node) {
                index = i;
                break;
            }
        }
        if (index < 0)  // node is not in nodes[]
            return;

        NetworkNode[] tmpNodes = new NetworkNode[nodeCount - 1];
        System.arraycopy(nodes, 0, tmpNodes, 0, index);
        System.arraycopy(nodes, index+1, tmpNodes, index, nodeCount-index-1);
        nodes = tmpNodes;

        // decrease the node count accordingly
        nodeCount--;
        if (index < leafNodeCount)
            leafNodeCount--;
        else if (index < leafNodeCount+speciationNodeCount)
            speciationNodeCount--;
        else if (index < nodeCount)
            reticulationNodeCount--;
    }

    /**
     * @return true if the network has a bubble
     */
    public boolean hasBubble() {
        for (NetworkNode hybridNode: getReticulationNodes()) {
            final int gammaBranchNr = hybridNode.gammaBranchNumber;
            if (hybridNode.getParentByBranch(gammaBranchNr) == hybridNode.getParentByBranch(gammaBranchNr + 1))
                return true;
        }
        return false;  // no parallel edges (bubble)
    }
}
