/*
 * PmTableModel.java
 *
 * Created on 19. November 2006, 13:26
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package ms_newpars;

import javax.swing.table.*;




/**
 *
 * @author berreth
 */

class PmTableModel extends AbstractTableModel {

    
    public    String headers[] = { "Potential Model", "Mol. Fract." };
    public    Object[][] data = 
        {
            {"PM_1 test", new Double(0.6)},
            {"PM_2 ", new Double(0.4)},

        };
        
    
    public int getColumnCount() {
        return headers.length;
    }

    public int getRowCount() {
        return data.length;
    }

    public String getColumnName(int col) {
        return headers[col];
    }

    public Object getValueAt(int row, int col) {
        return data[row][col];
    }

    public Class getColumnClass(int c) {
        return getValueAt(0, c).getClass();
    }

    /*
     * Don't need to implement this method unless your table's
     * editable.
     */
    public boolean isCellEditable(int row, int col) {
        //Note that the data/cell address is constant,
        //no matter where the cell appears onscreen.
        if (col < 2) {
            return false;
        } else {
            return true;
        }
    }

    /*
     * Don't need to implement this method unless your table's
     * data can change.
     */
    public void setValueAt(Object value, int row, int col) {
        data[row][col] = value;
        fireTableCellUpdated(row, col);
    }
    
    public void addRow(Object newRow[]) {
        Object newData[][] = new Object[getRowCount() + 1][getColumnCount()];
        for(int i = 0; i < getRowCount(); i++)
        newData = data;
        newData[getRowCount()] = newRow;
        data = newData;
    }
}
