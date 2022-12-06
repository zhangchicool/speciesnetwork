package speciesnetwork.app.beauti;

import beastfx.app.inputeditor.BEASTObjectInputEditor;
import beastfx.app.inputeditor.BeautiDoc;

import java.text.DateFormat;

public class SNTipDatesInputEditor extends BEASTObjectInputEditor {

    public SNTipDatesInputEditor(BeautiDoc doc) {
        super(doc);
    }
    
    private static final long serialVersionUID = 1L;

    DateFormat dateFormat = DateFormat.getDateInstance();


}
