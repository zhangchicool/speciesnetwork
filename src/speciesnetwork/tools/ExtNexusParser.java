package speciesnetwork.tools;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import beast.util.NexusParserListener;
import beast.util.TreeParser;
import beast.util.NexusParser;

/**
 * parses nexus file and grabs alignment and calibration from the file
 * Modified parseTreesBlock from NexusParser built into BEAST2 to support networks
 */
@Deprecated
public class ExtNexusParser extends NexusParser {

    @Override
    protected void parseTreesBlock(final BufferedReader fin) throws IOException {
        trees = new ArrayList<>();
        // read to first non-empty line within trees block
        String str = readLine(fin).trim();
        while (str.equals("")) {
            str = readLine(fin).trim();
        }

        int origin;
        // if first non-empty line is "translate" then parse translate block
        if (str.toLowerCase().contains("translate")) {
            translationMap = parseTranslateBlock(fin);
            origin = getIndexedTranslationMapOrigin(translationMap);
            if (origin != -1) {
                taxa = getIndexedTranslationMap(translationMap, origin);
            }
        }

        // read trees
        while (str != null) {
            if (str.toLowerCase().startsWith("tree ")) {
                final int i = str.indexOf('(');
                if (i > 0) {
                    str = str.substring(i);
                }
                final TreeParser treeParser = new TreeParser(str);  // modified here

                if (translationMap != null) treeParser.translateLeafIds(translationMap);

                // this needs to go after translation map or listeners have an incomplete tree!
                for (final NexusParserListener listener : listeners) {
                    listener.treeParsed(trees.size(), treeParser);
                }

                // this must come after listener or trees.size() gives the wrong index to treeParsed
                trees.add(treeParser);
            }
            str = fin.readLine();
            if (str != null) str = str.trim();
        }
    }
}
