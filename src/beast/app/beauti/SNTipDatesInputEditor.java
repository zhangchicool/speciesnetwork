package beast.app.beauti;

import java.awt.*;
import java.text.DateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.EventObject;
import java.util.GregorianCalendar;
import java.util.List;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.event.CellEditorListener;
import javax.swing.table.TableCellEditor;
import javax.swing.table.TableCellRenderer;

import beast.app.draw.BEASTObjectInputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.core.util.Log;
import beast.evolution.alignment.Taxon;
import beast.evolution.tree.TraitSet;

public class SNTipDatesInputEditor extends BEASTObjectInputEditor {

    public SNTipDatesInputEditor(BeautiDoc doc) {
        super(doc);
    }
    
    private static final long serialVersionUID = 1L;

    DateFormat dateFormat = DateFormat.getDateInstance();


}
