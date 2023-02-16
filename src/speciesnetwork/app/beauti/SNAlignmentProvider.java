package speciesnetwork.app.beauti;

import java.io.File;
import java.util.List;

import beast.base.core.BEASTInterface;
import beastfx.app.inputeditor.BeautiAlignmentProvider;
import beastfx.app.inputeditor.BeautiDoc;

public class SNAlignmentProvider extends BeautiAlignmentProvider {

	@Override
	public List<BEASTInterface> getAlignments(BeautiDoc doc, File[] files) {
		doc.autoSetClockRate = false;
		doc.beauti.autoSetClockRate(false);
		return super.getAlignments(doc, files);
	}

}
