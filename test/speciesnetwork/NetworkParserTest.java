package speciesnetwork;

import java.text.DecimalFormat;
import org.junit.Test;
import static org.junit.Assert.assertEquals;

import beast.util.TreeParser;
import junit.framework.TestCase;
import speciesnetwork.NetworkParser;

public class NetworkParserTest extends TestCase {
    TreeParser treeParser;
    NetworkParser networkParser;

    // each node should have a label to be parsed correctly
    final String testNetwork = "(((A:0.2,(B:0.1)#H1[&gamma=0.9]:0.1)S1:0.3,(#H1:0.2,C:0.3)S2:0.2)R:0.1)";
    final String testNetwork2= "((((A:0.02,(B:0.01)#H1[&gamma=0.5]:0.01)S3:0.01,(#H1:0.01)#H2[&gamma=0.5]:0.01)S2:0.02,(#H2:0.02,C:0.04)S1:0.01)R:0.01)";
    final String testNetwork3 = "(((((A:0.1)#H1[&gamma=0.9]:0.1)#H2[&gamma=0.8]:0.3,((#H2:0.1,(#H1:0.1)#H3[&gamma=0.7]:0.1)S1:0.1)#H4[&gamma=0.6]:0.1)S2:0.1,(#H4:0.1,#H3:0.3)S3:0.1)R:0.1)";

    public NetworkParserTest() {
        treeParser = new TreeParser();
        networkParser = new NetworkParser();

        try {
            treeParser.initByName("newick", testNetwork, "IsLabelledNewick", true, "adjustTipHeights", false);
            networkParser.initByName("tree", treeParser);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Test
    public void testParser() {
        DecimalFormat df = new DecimalFormat("0.####");;
        assertEquals(testNetwork, networkParser.toString(df));
    }
}
