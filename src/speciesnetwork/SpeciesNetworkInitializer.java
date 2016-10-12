package speciesnetwork;

import java.util.*;
import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.*;
import beast.evolution.alignment.distance.Distance;
import beast.evolution.alignment.distance.JukesCantorDistance;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.ClusterTree;
import speciesnetwork.operators.RebuildEmbedding;

/**
 * @author Joseph Heled
 * @author Huw Ogilvie
 * @author Chi Zhang
 */

@Description("Set a starting point for a species-network analysis.")
public class SpeciesNetworkInitializer extends Tree implements StateNodeInitialiser {

    private enum Method {
        POINT("point"),   // point estimate from the sequences
        RANDOM("random"), // random starting trees (caterpillar)
        USER("user");     // user defined starting state

        Method(final String name) {
            this.ename = name;
        }

        @Override
		public String toString() {
            return ename;
        }

        private final String ename;
    }
    public final Input<Method> initMethod = new Input<>("method", "Initialise either with a totally random state" +
            "or a point estimate based on alignments data (default point)", Method.POINT, Method.values());
    public final Input<Network> speciesNetworkInput
            = new Input<>("speciesNetwork", "Species network to initialize.", Validate.REQUIRED);
    public final Input<RealParameter> originInput =
            new Input<>("origin", "The time when the process started.", Validate.REQUIRED);
    public final Input<List<Tree>> geneTreesInput =
            new Input<>("geneTree", "Gene tree to initialize.", new ArrayList<>());
    public final Input<List<RebuildEmbedding>> rebuildEmbeddingInput = new Input<>("rebuildEmbedding",
            "Operator which rebuilds embedding of gene trees within species network.", new ArrayList<>());
    public final Input<YuleHybridModel> hybridYuleInput = new Input<>("hybridYule",
            "The species network (with hybridization) to initialize.", Validate.XOR, speciesNetworkInput);
    public final Input<RealParameter> birthRateInput =
            new Input<>("birthRate", "Network birth rate to initialize.");
    public final Input<RealParameter> hybridRateInput =
            new Input<>("hybridRate", "Network hybridization rate to initialize.");
    public final Input<Function> clockRateInput =
            new Input<>("baseRate", "Main clock rate used to scale trees (default 1).");

    @Override
    public void initAndValidate() {
        // what does this do and is it dangerous to call it or not to call it at the start or at the end??????
        super.initAndValidate();
    }

    @Override
    public void initStateNodes() {
        final Method method = initMethod.get();
        switch( method ) {
            case POINT:
                pointInit();
                break;
            case RANDOM:
                randomInit();
                break;
            case USER:
                userStates();
                break;
        }

        final double tOrigin = originInput.get().getValue();
        final double tMRCA = speciesNetworkInput.get().getRoot().getHeight();
        if (tOrigin < tMRCA)
            throw new IllegalArgumentException("Time of origin (" + tOrigin + ") < time of MRCA (" + tMRCA + ")!");

        // initialize embedding for all gene trees
        for (RebuildEmbedding operator: rebuildEmbeddingInput.get()) {
            if (operator.initializeEmbedding() < 0)
                throw new RuntimeException("Failed to build gene tree embedding!");
        }
    }

    private void userStates() {
        // initialize parameters using user defined starting state
        // nothing to do here at this moment
    }

    private void randomInit() {
        // initialize caterpillar species tree, if no user defined network
        final Network sNetwork = speciesNetworkInput.get();
        // scale the species network according to the time of origin
        final double tOrigin = originInput.get().getValue();
        sNetwork.scale(tOrigin/sNetwork.getOrigin().getHeight());

        final double rootHeight = sNetwork.getRoot().getHeight();

        // initialize caterpillar gene trees
        final List<Tree> geneTrees = geneTreesInput.get();
        for (final Tree gtree : geneTrees) {
            gtree.makeCaterpillar(rootHeight, rootHeight/gtree.getInternalNodeCount(), true);
        }
    }

    /**
     * Build gene trees and species tree from alignments (copied from starbeast)
     */
    private void pointInit() {
        final double clockRate = (clockRateInput.get() != null) ? clockRateInput.get().getArrayValue() : 1;

        final Network sNetwork = speciesNetworkInput.get();
        final TaxonSet species = sNetwork.taxonSetInput.get();
        final List<String> speciesNames = species.asStringList();
        final int speciesCount = speciesNames.size();

        final List<Tree> geneTrees = geneTreesInput.get();
        double maxNsites = 0;
        for (final Tree gtree : geneTrees) {
            final Alignment alignment = gtree.m_taxonset.get().alignmentInput.get();
            final ClusterTree ctree = new ClusterTree();
            ctree.initByName("initial", gtree, "clusterType", "upgma", "taxa", alignment);
            gtree.scale(1 / clockRate);
            maxNsites = max(maxNsites, alignment.getSiteCount());
        }

        final Map<String, Integer> geneTips2Species = new LinkedHashMap<>();
        final List<Taxon> taxonSets = species.taxonsetInput.get();
        for(int k = 0; k < speciesNames.size(); ++k) {
            final Taxon nx = taxonSets.get(k);
            final List<Taxon> taxa = ((TaxonSet) nx).taxonsetInput.get();
            for(final Taxon n : taxa) {
                geneTips2Species.put(n.getID(), k);
            }
        }

        final double[] dg = new double[speciesCount*(speciesCount-1)/2];
        final double[][] genesDmins = new double[geneTrees.size()][];
        for(int ng = 0; ng < geneTrees.size(); ++ng) {
            final Tree g = geneTrees.get(ng);
            final double[] dmin = firstMeetings(g, geneTips2Species, speciesCount);
            genesDmins[ng] = dmin;

            for(int i = 0; i < dmin.length; ++i) {
                dg[i] += dmin[i];
                if (dmin[i] == Double.MAX_VALUE) {
                    // this happens when a gene tree has no taxa for some species-tree taxon.
                    // TODO: ensure that if this happens, there will always be an "infinite"
                    // distance between species-taxon 0 and the species-taxon with missing lineages,
                    // so i < speciesCount - 1.
                    // What if lineages for species-taxon 0 are missing? Then all entries will be 'infinite'.
                    String id = (i < speciesCount - 1? sNetwork.nodes[i+1].label : "unknown taxon");
                    if (i == 0) {
                        // test that all entries are 'infinite', which implies taxon 0 has lineages missing
                        boolean b = true;
                        for (int k = 1; b && k < speciesCount - 1; ++k) {
                            b = (dmin[k] == Double.MAX_VALUE);
                        }
                        if (b) {
                            // if all entries have 'infinite' distances, it is probably the first taxon that is at fault
                            id = sNetwork.nodes[0].label;
                        }
                    }
                    throw new RuntimeException("Gene tree " + g.getID() + " has no lineages for species taxon " + id + " ");
                }
            }
        }

        for(int i = 0; i < dg.length; ++i) {
            double d = dg[i] / geneTrees.size();
            if( d == 0 ) {
                d = (0.5/maxNsites) * (1/clockRate);
            } else {
                // heights to distances
                d *= 2;
            }
            dg[i] = d;
        }

        final ClusterTree stree = new ClusterTree();
        final Distance distance = new Distance() {
            @Override
            public double pairwiseDistance(final int s1, final int s2) {
                final int i = getDMindex(speciesCount, s1,s2);
                return dg[i];
            }
        };
        Tree speciesTree = new Tree();
        stree.initByName("initial", speciesTree, "taxonset", species,"clusterType", "upgma", "distance", distance);

        final Map<String, Integer> sptips2SpeciesIndex = new LinkedHashMap<>();
        for(int i = 0; i < speciesNames.size(); ++i) {
            sptips2SpeciesIndex.put(speciesNames.get(i), i);
        }
        final double[] spmin = firstMeetings(speciesTree, sptips2SpeciesIndex, speciesCount);

        for(int ng = 0; ng < geneTrees.size(); ++ng) {
            final double[] dmin = genesDmins[ng];
            boolean compatible = true;
            for(int i = 0; i < spmin.length; ++i) {
                if( dmin[i] <= spmin[i] ) {
                    compatible = false;
                    break;
                }
            }
            if( ! compatible ) {
                final Tree gtree = geneTrees.get(ng);
                final TaxonSet gtreeTaxa = gtree.m_taxonset.get();
                final Alignment alignment = gtreeTaxa.alignmentInput.get();
                final List<String> taxaNames = alignment.getTaxaNames();
                final int taxonCount =  taxaNames.size();
                // speedup
                final Map<Integer,Integer> g2s = new LinkedHashMap<>();
                for(int i = 0; i < taxonCount; ++i) {
                    g2s.put(i, geneTips2Species.get(taxaNames.get(i)));
                }

                final JukesCantorDistance jc = new JukesCantorDistance();
                jc.setPatterns(alignment);
                final Distance gdistance = new Distance() {
                    @Override
                    public double pairwiseDistance(final int t1, final int t2) {
                        final int s1 = g2s.get(t1);
                        final int s2 = g2s.get(t2);
                        double d = jc.pairwiseDistance(t1,t2) / clockRate;
                        if( s1 != s2 ) {
                            final int i = getDMindex(speciesCount, s1,s2);
                            final double minDist = 2 * spmin[i];
                            if( d <= minDist ) {
                                d = minDist * 1.001;
                            }
                        }
                        return d;
                    }
                };
                final ClusterTree ctree = new ClusterTree();
                ctree.initByName("initial", gtree, "taxonset", gtreeTaxa, "clusterType", "upgma", "distance", gdistance);
            }
        }

        // finally, convert the species tree to species network
        final Node root = speciesTree.getRoot();
        final Node origin = newNode();       // create origin
        origin.setHeight(1.2 * root.getHeight());
        origin.setLeft(root);                // set root as child of origin
        root.setParent(origin);              // set origin as parent of root
        speciesTree.getInternalNodeCount();  // make sure node counts are correct
        speciesTree.addNode(origin);         // add the origin
        speciesTree.setRoot(origin);         // set origin as new root
        NetworkParser networkParser = new NetworkParser(speciesTree);
        sNetwork.assignFrom(networkParser);
        sNetwork.resetInternalNodeLabels();

        // set the time of origin
        final double tOrigin = sNetwork.getOrigin().getHeight();
        final RealParameter originTime = originInput.get();
        originTime.setValue(tOrigin);
    }

    private double[] firstMeetings(final Tree gtree, final Map<String, Integer> tipName2Species, final int speciesCount) {
        final Node[] nodes = gtree.listNodesPostOrder(null, null);
        @SuppressWarnings("unchecked")
        final Set<Integer>[] tipsSpecies = new Set[nodes.length];
        for(int k = 0; k < tipsSpecies.length; ++k) {
            tipsSpecies[k] = new LinkedHashSet<>();
        }
        // d[i,j] = minimum height of node which has tips belonging to species i and j
        // d is is upper triangular
        final double[] dmin = new double[(speciesCount*(speciesCount-1))/2];
        Arrays.fill(dmin, Double.MAX_VALUE);

        for (final Node n : nodes) {
            if (n.isLeaf()) {
                tipsSpecies[n.getNr()].add(tipName2Species.get(n.getID()));
            } else {
                assert n.getChildCount() == 2;
                @SuppressWarnings("unchecked")
                final Set<Integer>[] sps = new Set[2];
                sps[0] = tipsSpecies[n.getChild(0).getNr()];
                sps[1] = tipsSpecies[n.getChild(1).getNr()];
                final Set<Integer> u = new LinkedHashSet<>(sps[0]);
                u.retainAll(sps[1]);
                sps[0].removeAll(u);
                sps[1].removeAll(u);

                for (final Integer s1 : sps[0]) {
                    for (final Integer s2 : sps[1]) {
                        final int i = getDMindex(speciesCount, s1, s2);
                        dmin[i] = min(dmin[i], n.getHeight());
                    }
                }
                u.addAll(sps[0]);
                u.addAll(sps[1]);
                tipsSpecies[n.getNr()] = u;
            }
        }
        return dmin;
    }

    private int getDMindex(final int speciesCount, final int s1, final int s2) {
        final int mij = min(s1,s2);
        return mij*(2*speciesCount-1 - mij)/2 + (abs(s1-s2)-1);
    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(speciesNetworkInput.get());
        for(final Tree gtree : geneTreesInput.get()) {
            stateNodes.add(gtree);
        }

        final RealParameter brate = birthRateInput.get();
        if(brate != null) stateNodes.add(brate);
        final RealParameter hrate = hybridRateInput.get();
        if(hrate != null) stateNodes.add(hrate);
    }
}
