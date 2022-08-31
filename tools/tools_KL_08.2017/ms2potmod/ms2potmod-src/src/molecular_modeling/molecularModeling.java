/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package molecular_modeling;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Container;
import java.awt.Desktop;
import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.BufferedReader;
import java.io.BufferedWriter; 
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.StringTokenizer;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.DefaultCellEditor;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.KeyStroke;
import javax.swing.table.TableColumn;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.plaf.basic.BasicInternalFrameUI;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableColumnModel;
import javax.vecmath.Point3f;
import javax.vecmath.Vector3f;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

/**
 *
 * @author Syed ahsan ali
 * 
 * This is the main class for the overall functionality and interface provided
 * by ms2potmod.
 */
public class molecularModeling extends javax.swing.JFrame {

    /**
     * Creates new form molecularModeling
     */
    String UserDirectory = null;
    File OpenSelectedDir = null;
    File SaveSelectedDir = null;
    String[] DBvalues;
    JTable  ElementsTable = new JTable();
    EachRowEditor rowEditor,TrowEditor;
    List<JComboBox> Jbox = new ArrayList<JComboBox>();
    List<JComboBox> TJbox = new ArrayList<JComboBox>();
    HashMap<String,String> helpMap = new HashMap<String,String>();
    LJSites ljSite = new LJSites();
    
    int Bond3count = 0, Bond2count = 0, Bond4count = 0;
    int[][] bondlist, anglelist,dihedlist;
    double[] bondlength, bondconst,angle,angleconst;
    double[][] dihed,dihedconst;
    float max = 0.0f;
   // boolean ISGeom = false;
    
    /**
     * Constructor. Initializes and settings of some parameters
     * necessery for program.
     */
    public molecularModeling() {
        initComponents(); // Initializes UI
        SetFrame(); // setting frame size according to screen
        setTitle("MS2PotMod"); // set the title
        CreateElementsTable(); // create main table
        Toolkit.getDefaultToolkit().sync();
        addTablelistener(); // add listner to table
        UserDirectory = System.getProperty("user.dir"); // save user directory
        
        // add editor for dropdown
        try {
            rowEditor = new EachRowEditor(ElementsTable);
            TrowEditor = new EachRowEditor(ElementsTable);
            ElementsTable.getColumn("Database").setCellEditor(rowEditor);
            ElementsTable.getColumn("Type").setCellEditor(TrowEditor);
            DBvalues = getDBsNames("");
            addComboTotable();
        } catch (Exception ex) {
            System.out.println(ex.getStackTrace());
        }
        
        StartSceneGraphWindow(); // Start canvas to draw structure
        readHelpxlsxfile(); // read help file to add tootips
        settooltipforlabels(); // set tooltips for all controls
        addShortCutsForButtonActions(); // assign shortcut keys
    }
    
    /*
        assign shortcut keys for open, save and help buttons
    */
private void addShortCutsForButtonActions() {
    //Short cut for save button
    Action savebuttonAction = new AbstractAction("Save") {
        @Override
        public void actionPerformed(ActionEvent evt) {
            SavebtnActionPerformed(evt);
        }
    };
 
    this.Savebtn.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(
        KeyStroke.getKeyStroke(KeyEvent.VK_S, KeyEvent.CTRL_DOWN_MASK), "Save");
    this.Savebtn.getActionMap().put("Save", savebuttonAction);

    //Shortcut for load button
    Action loadbuttonAction = new AbstractAction("Load") {
        @Override
        public void actionPerformed(ActionEvent evt) {
           OpenBtnActionPerformed(evt);
        }
    };
 
    this.OpenBtn.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(
        KeyStroke.getKeyStroke(KeyEvent.VK_O, KeyEvent.CTRL_DOWN_MASK), "Load");
    this.OpenBtn.getActionMap().put("Load", loadbuttonAction);
    
   //shortcut for help button
    Action helpbuttonAction = new AbstractAction("Help") {

        @Override
        public void actionPerformed(ActionEvent evt) {
            helpbtnActionPerformed(evt);
        }
    };

    this.helpbtn.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(
            KeyStroke.getKeyStroke(KeyEvent.VK_F1, 0), "Help");
    this.helpbtn.getActionMap().put("Help", helpbuttonAction);

}
    /**
     * Get all the labels and buttons to assign tootip text 
     * from the hashmap (helpmap).
     * 
    */
    private void settooltipforlabels()
    {
        for(Component cm: mainmenu.getComponents())
        {
            if(cm instanceof JLabel)
            {
                JLabel lbl = (JLabel) cm;
                String tooltip = addLinebreaks(helpMap.get(lbl.getText().toLowerCase().trim()),30);
                lbl.setToolTipText("<html><body><p>"+tooltip+"</p></body></html>");
            }
            else
            if(cm instanceof JButton)
            {
                JButton btn = (JButton) cm;
                String tooltip = addLinebreaks(helpMap.get(btn.getText().toLowerCase().trim()),30);
                btn.setToolTipText("<html><body><p>"+tooltip+"</p></body></html>");
            }
        }
        
        for(Component cm: canvasmenu.getComponents())
        {
            if(cm instanceof JLabel)
            {
                JLabel lbl = (JLabel) cm;
                String tooltip = addLinebreaks(helpMap.get(lbl.getText().toLowerCase().trim()),30);
                lbl.setToolTipText("<html><body><p>"+tooltip+"</p></body></html>");
            }
            else
            if(cm instanceof JButton)
            {
                JButton btn = (JButton) cm;
                String tooltip = addLinebreaks(helpMap.get(btn.getText().toLowerCase().trim()),30);
                btn.setToolTipText("<html><body><p>"+tooltip+"</p></body></html>");
            }
        }
    }
    /*
    
    */

    /**
     * Add line breaks for the tooltip text to make it multiline easy to read.
     * 
     * @param input Help text
     * @param maxLineLength Length of a line
     * 
     * @return Formatted string
     */
    
     public String addLinebreaks(String input, int maxLineLength) {
    if(input != null)
    {
    StringTokenizer tok = new StringTokenizer(input, " ");
    StringBuilder output = new StringBuilder(input.length());
    int lineLen = 0;
    while (tok.hasMoreTokens()) {
        String word = tok.nextToken();

        if (lineLen + word.length() > maxLineLength) {
            output.append("<br>");
            lineLen = 0;
        }
        output.append(word+" ");
        lineLen += word.length();
    }
    return output.toString();
    }
    
    return "";
}
     /**
      * 
     * read the .xlsx help file
     * find the name and meaning coloumn
     * populate the hashmap (helpmap)
     * key: name coloumn from the help file
     * value: meaning coloumn from the help file.
     * 
     */
    private void readHelpxlsxfile()
    {
        try 
        {
        // Get the workbook instance for XLS file
        String filePath = UserDirectory+""+File.separator+"help"+File.separator+"ms2potmod_help_format.xlsx";
        XSSFWorkbook workbook = new XSSFWorkbook(new FileInputStream(filePath));
        // Get first sheet from the workbook
        XSSFSheet sheet = workbook.getSheetAt(0);
        Cell cell;
        Row row;
        int namePOS = -1, meanPOS = -1;
        row = sheet.getRow(0);
       for(int i=0; i<row.getLastCellNum();i++)
       {
           cell = row.getCell(i);
           if(cell.getStringCellValue().toLowerCase().equals("name"))
               namePOS = i;
           if(cell.getStringCellValue().toLowerCase().equals("meaning"))
               meanPOS = i;
       }
       
       // Iterate through each rows from first sheet
        Iterator<Row> rowIterator = sheet.iterator();
        while (rowIterator.hasNext()) 
        {
                row = rowIterator.next();
                cell = row.getCell(namePOS);
                Cell mcell = row.getCell(meanPOS); 
         
            if(cell != null && mcell != null)
               if(cell.getCellType() != Cell.CELL_TYPE_BLANK && mcell.getCellType() != Cell.CELL_TYPE_BLANK)
                   helpMap.put(cell.getStringCellValue().toLowerCase().trim(), mcell.getStringCellValue());
                  // System.out.println("Key : "+cell.getStringCellValue()+" | Value : "+mcell.getStringCellValue());
         }
        
    } 
        catch (FileNotFoundException e) 
        {
                System.err.println("Exception" + e.getMessage());
        }
        catch (IOException e) 
        {
                System.err.println("Exception" + e.getMessage());
        }
        catch(Exception e)
        {
            System.err.println("Exception" + e.getMessage());
        }
    }
    
    /**
     * 
     * Create the main elements table
     * set its properties accordingly.
     * 
    */
    private void CreateElementsTable()
    {
      
      JScrollPane JS = new JScrollPane();
        ElementsTable.setModel(new DefaultTableModel(
            new Object [][] {
                {"0", "0.0", "0.0", "0.0", "", "", null, null, null, null}
            },
            new String [] {
                "ID", "X", "Y", "Z", "Link", "Type", "Database", "Sigma", "Epsilon", "Mass"
            }
        ) {
            boolean[] canEdit = new boolean [] {
                false, true, true, true, false, true, true, false, false, false
            };

            public boolean isCellEditable(int rowIndex, int columnIndex) {
                return canEdit [columnIndex];
            }
            
        });
        ElementsTable.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);

        JS.setViewportView(ElementsTable);
        if (ElementsTable.getColumnModel().getColumnCount() > 0) {
            ElementsTable.getColumnModel().getColumn(0).setResizable(false);
            ElementsTable.getColumnModel().getColumn(0).setPreferredWidth(40);
            ElementsTable.getColumnModel().getColumn(1).setPreferredWidth(100);
            ElementsTable.getColumnModel().getColumn(2).setPreferredWidth(100);
            ElementsTable.getColumnModel().getColumn(3).setPreferredWidth(100);
            ElementsTable.getColumnModel().getColumn(4).setPreferredWidth(100);
            ElementsTable.getColumnModel().getColumn(5).setPreferredWidth(100);
            ElementsTable.getColumnModel().getColumn(6).setPreferredWidth(250);
            ElementsTable.getColumnModel().getColumn(7).setPreferredWidth(100);
            ElementsTable.getColumnModel().getColumn(8).setPreferredWidth(100);
            ElementsTable.getColumnModel().getColumn(9).setPreferredWidth(100);
        }
        
        this.leftPanelM.add(JS,BorderLayout.CENTER);
    }
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        ElementTablePopupmenu = new javax.swing.JPopupMenu();
        jPopupMenu1 = new javax.swing.JPopupMenu();
        jMenu3 = new javax.swing.JMenu();
        leftPanelM = new javax.swing.JPanel();
        mainmenu = new javax.swing.JToolBar();
        OpenBtn = new javax.swing.JButton();
        Savebtn = new javax.swing.JButton();
        addDBentry = new javax.swing.JButton();
        RefreshDBbtn = new javax.swing.JButton();
        ljsitebtn = new javax.swing.JButton();
        Restartbtn = new javax.swing.JButton();
        helpbtn = new javax.swing.JButton();
        canvasmenu = new javax.swing.JToolBar();
        undobtn = new javax.swing.JButton();
        resetbtn = new javax.swing.JButton();
        centerofmassbtn = new javax.swing.JButton();
        CalculateGeometry = new javax.swing.JButton();
        error = new javax.swing.JLabel();
        errorlbl = new javax.swing.JLabel();
        CanvasInternalFrame = new javax.swing.JInternalFrame();
        logpanel = new javax.swing.JPanel();
        jScrollPane2 = new javax.swing.JScrollPane();
        jScrollPane1 = new javax.swing.JScrollPane();
        log = new javax.swing.JTextArea();

        jMenu3.setText("jMenu3");

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        addComponentListener(new java.awt.event.ComponentAdapter() {
            public void componentResized(java.awt.event.ComponentEvent evt) {
                formComponentResized(evt);
            }
        });

        leftPanelM.setLayout(new java.awt.BorderLayout());

        mainmenu.setBorder(javax.swing.BorderFactory.createEtchedBorder());
        mainmenu.setFloatable(false);

        OpenBtn.setText("  Open  ");
        OpenBtn.setToolTipText("Open pm file");
        OpenBtn.setFocusable(false);
        OpenBtn.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        OpenBtn.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        OpenBtn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                OpenBtnActionPerformed(evt);
            }
        });
        mainmenu.add(OpenBtn);

        Savebtn.setText("  Save  ");
        Savebtn.setToolTipText("Save parameters to pm file");
        Savebtn.setFocusable(false);
        Savebtn.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        Savebtn.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        Savebtn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                SavebtnActionPerformed(evt);
            }
        });
        mainmenu.add(Savebtn);

        addDBentry.setText("Add Database entry");
        addDBentry.setToolTipText("Add new entries to the forcefield database");
        addDBentry.setFocusable(false);
        addDBentry.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        addDBentry.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        addDBentry.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                addDBentryActionPerformed(evt);
            }
        });
        mainmenu.add(addDBentry);

        RefreshDBbtn.setText("  Refresh Database");
        RefreshDBbtn.setToolTipText("Refresh the database column in the table");
        RefreshDBbtn.setFocusable(false);
        RefreshDBbtn.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        RefreshDBbtn.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        RefreshDBbtn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                RefreshDBbtnActionPerformed(evt);
            }
        });
        mainmenu.add(RefreshDBbtn);

        ljsitebtn.setText("  Change unit/site");
        ljsitebtn.setFocusable(false);
        ljsitebtn.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        ljsitebtn.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        ljsitebtn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ljsitebtnActionPerformed(evt);
            }
        });
        mainmenu.add(ljsitebtn);

        Restartbtn.setText("  Restart  ");
        Restartbtn.setToolTipText("Restart the program all over again");
        Restartbtn.setFocusable(false);
        Restartbtn.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        Restartbtn.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        Restartbtn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                RestartbtnActionPerformed(evt);
            }
        });
        mainmenu.add(Restartbtn);

        helpbtn.setText("  Help  ");
        helpbtn.setToolTipText("See help file");
        helpbtn.setBorder(null);
        helpbtn.setFocusable(false);
        helpbtn.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        helpbtn.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        helpbtn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                helpbtnActionPerformed(evt);
            }
        });
        mainmenu.add(helpbtn);

        canvasmenu.setBorder(javax.swing.BorderFactory.createEtchedBorder());
        canvasmenu.setFloatable(false);
        canvasmenu.setRollover(true);

        undobtn.setText("  Undo  ");
        undobtn.setToolTipText("delete only currently added atom");
        undobtn.setBorder(null);
        undobtn.setEnabled(false);
        undobtn.setFocusable(false);
        undobtn.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        undobtn.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        undobtn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                undobtnActionPerformed(evt);
            }
        });
        canvasmenu.add(undobtn);

        resetbtn.setText("Reset view");
        resetbtn.setToolTipText("Reset the structure on its orignal position");
        resetbtn.setBorder(null);
        resetbtn.setFocusable(false);
        resetbtn.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        resetbtn.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        resetbtn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                resetbtnActionPerformed(evt);
            }
        });
        canvasmenu.add(resetbtn);

        centerofmassbtn.setText("  Set Center of Mass");
        centerofmassbtn.setToolTipText("Set angle of rotation according to the center of mass");
        centerofmassbtn.setBorder(null);
        centerofmassbtn.setFocusable(false);
        centerofmassbtn.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        centerofmassbtn.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        centerofmassbtn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                centerofmassbtnActionPerformed(evt);
            }
        });
        canvasmenu.add(centerofmassbtn);

        CalculateGeometry.setText("Calculate Geometry");
        CalculateGeometry.setToolTipText("Calculate x,y,z coordinates according to the structure");
        CalculateGeometry.setBorder(null);
        CalculateGeometry.setFocusable(false);
        CalculateGeometry.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        CalculateGeometry.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        CalculateGeometry.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                CalculateGeometryActionPerformed(evt);
            }
        });
        canvasmenu.add(CalculateGeometry);

        error.setHorizontalAlignment(javax.swing.SwingConstants.RIGHT);
        error.setText("  Error");
        error.setToolTipText("");
        canvasmenu.add(error);

        errorlbl.setText("::[ BE: , AE: , DE: ]");
        canvasmenu.add(errorlbl);

        CanvasInternalFrame.setCursor(new java.awt.Cursor(java.awt.Cursor.HAND_CURSOR));
        CanvasInternalFrame.setVisible(true);

        javax.swing.GroupLayout CanvasInternalFrameLayout = new javax.swing.GroupLayout(CanvasInternalFrame.getContentPane());
        CanvasInternalFrame.getContentPane().setLayout(CanvasInternalFrameLayout);
        CanvasInternalFrameLayout.setHorizontalGroup(
            CanvasInternalFrameLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 1613, Short.MAX_VALUE)
        );
        CanvasInternalFrameLayout.setVerticalGroup(
            CanvasInternalFrameLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 1369, Short.MAX_VALUE)
        );

        log.setEditable(false);
        log.setBackground(new java.awt.Color(0, 0, 0));
        log.setColumns(20);
        log.setFont(log.getFont().deriveFont(log.getFont().getStyle() | java.awt.Font.BOLD, log.getFont().getSize()+6));
        log.setForeground(new java.awt.Color(255, 255, 255));
        log.setRows(5);
        log.setText("Log:");
        jScrollPane1.setViewportView(log);

        jScrollPane2.setViewportView(jScrollPane1);

        javax.swing.GroupLayout logpanelLayout = new javax.swing.GroupLayout(logpanel);
        logpanel.setLayout(logpanelLayout);
        logpanelLayout.setHorizontalGroup(
            logpanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jScrollPane2)
        );
        logpanelLayout.setVerticalGroup(
            logpanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, logpanelLayout.createSequentialGroup()
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(jScrollPane2, javax.swing.GroupLayout.PREFERRED_SIZE, 206, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(10, 10, 10))
        );

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(logpanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(leftPanelM, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(mainmenu, javax.swing.GroupLayout.DEFAULT_SIZE, 963, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(canvasmenu, javax.swing.GroupLayout.DEFAULT_SIZE, 1629, Short.MAX_VALUE)
                            .addComponent(CanvasInternalFrame))))
                .addGap(12, 12, 12))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addGap(8, 8, 8)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(mainmenu, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(canvasmenu, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(CanvasInternalFrame)
                    .addComponent(leftPanelM, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(logpanel, javax.swing.GroupLayout.PREFERRED_SIZE, 231, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(4, 4, 4))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents
    /**
     * 
     * To reset the SceneGraph to its original position
     * 
    */
    private void resetbtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_resetbtnActionPerformed

        CreateGraphics.ResetScene();
        
    }//GEN-LAST:event_resetbtnActionPerformed
    
    /**
     * 
     * undo button event handler to delete the most recent added atom or shape:
     * Calls undoscene function to delete the shapes from scenegraph 
     * delete the row from elements table 
     * And delete the its lisnks from the ID where it attached to.
     * 
    */
    private void undobtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_undobtnActionPerformed
       try{
        CreateGraphics.undoScene();
        DefaultTableModel tm = (DefaultTableModel) ElementsTable.getModel();
        int rowCount = tm.getRowCount();
        int LinkID = Integer.parseInt(tm.getValueAt(rowCount - 1, 4).toString());
        if(rowCount > 2)
            tm.setValueAt(tm.getValueAt(LinkID, 4).toString().substring(0,tm.getValueAt(LinkID, 4).toString().lastIndexOf(",")), LinkID, 4);
        else
            tm.setValueAt("",LinkID,4);
        
        tm.removeRow(rowCount-1);
        undobtn.setEnabled(false);
       }
       catch(Exception ex)
       {
           ex.printStackTrace();
           log.append("!!! Caught Exception while deleting latest atom\n");
       }
    }//GEN-LAST:event_undobtnActionPerformed
    
    /**
     * 
     * Save button event handler.
     * 
    */
    
    private void SavebtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_SavebtnActionPerformed
       
        JFileChooser pmfilechooser = new JFileChooser();
        pmfilechooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
        pmfilechooser.setAcceptAllFileFilterUsed(false);
        pmfilechooser.addChoosableFileFilter(new FileNameExtensionFilter("PotMoD Files", "pm"));
        pmfilechooser.setSelectedFile(new File("Default.pm"));
        if(SaveSelectedDir != null)
             pmfilechooser.setCurrentDirectory(SaveSelectedDir);
        if(pmfilechooser.showSaveDialog(this) == JFileChooser.APPROVE_OPTION)
        {
            SaveSelectedDir = pmfilechooser.getSelectedFile();
            log.append(">>> Saving directory selected...\n");
            if(!IsDatabaseValueSelected()){
                JOptionPane.showMessageDialog(this, "Error Occured. Can not write file. one or more Database value are missing ","Error",JOptionPane.ERROR_MESSAGE);
                log.append("!!! Error Occured. Can not write file. one or more Database value are missing\n");
            }
            else
                SavepmFile(pmfilechooser.getSelectedFile());
            
        }
    }//GEN-LAST:event_SavebtnActionPerformed
    
    /**
     * 
     * To check if the value is selected in the
     * database column of the elements table.
     * 
     * @return True - False
     * 
    */
    private boolean IsDatabaseValueSelected()
    {
        for(int i = 0; i < ElementsTable.getRowCount(); i++)
        {
            if(ElementsTable.getValueAt(i, 6) == null)
                        return false;
        }
        
        return true;
    }
    
    /**
     * 
     * Algorithm to calculate the center of mass for molecular 
     * structure and sets it as the center of mouse rotation.
     * 
     * @see CreateGraphics#setCenterOfMass(float, float, float)
     * 
    */
    private void centerofmassbtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_centerofmassbtnActionPerformed
        // TODO add your handling code here:
        float totalmass = 0.0f,Xcm = 0.0f,Ycm = 0.0f,Zcm = 0.0f;
        boolean flag = false;
        for(int i=0;i<this.ElementsTable.getRowCount();i++)
        {
            if(this.ElementsTable.getValueAt(i, 9) != null)
            {
            float mass = Float.parseFloat(this.ElementsTable.getValueAt(i, 9).toString());
            float x    = Float.parseFloat(this.ElementsTable.getValueAt(i, 1).toString());
            float y    = Float.parseFloat(this.ElementsTable.getValueAt(i, 2).toString());
            float z    = Float.parseFloat(this.ElementsTable.getValueAt(i, 3).toString());
            if(max != 0.0)
            {
                x = x/max;
                y = y/max;
                z = z/max;
            }
            else
                if(!(x > -1 && x < 1) || !(y > -1 && y < 1) || !(z > -1 && z < 1))
                   {        
                       float maxi = Math.max(Math.abs(x), Math.abs(y));
                       maxi = Math.max(maxi, Math.abs(z));
                       maxi = maxi + 1;
                       // System.out.println("max value : "+maxi);
                       x = x/maxi;
                       y = y/maxi;
                       z = z/maxi;
                                 
                      // System.out.println("x ="+x+"  y ="+y+"  z ="+z);
                    }
            totalmass += mass;
//            MandRmul += 0.035 * mass;
            Xcm += x * mass;
            Ycm += y * mass;
            Zcm += z * mass;
            }
            else
                flag = true;
               
        }
        if(flag)
        {
           JOptionPane.showMessageDialog(null,"One of the ID has no mass. Please select Database...");
           log.append("!!! One of the ID has no mass. Please select Database...\n");
        }
        else
        {
        Xcm = Xcm / totalmass;
        Ycm = Ycm / totalmass;
        Zcm = Zcm / totalmass;
        //JOptionPane.showMessageDialog(null, Xcm+","+Ycm+","+Zcm);
        CreateGraphics.setCenterOfMass(Xcm, Ycm, Zcm);
        log.append(">>> Center of mass has been changed ["+Xcm+" , "+Ycm+" , "+Zcm+"]\n");
        }
    }//GEN-LAST:event_centerofmassbtnActionPerformed

    /**
     * 
     * event handler for Geometry button:
     * Calls BADExtraction class to extract all the
     * bonds, angles and dihedrals
     * Call MolDyn class to calculate geometry from these parameters
     * Remove existing molecular structure from the scenegraph
     * Calculates new parameters and redraw structure with the help
     * of createGraphics class.
     * 
     * @see BADExtraction
     * @see MolDyn
     * @see CreateGraphics
     * 
    */
    private void CalculateGeometryActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_CalculateGeometryActionPerformed
        
        if(ElementsTable.getRowCount() > 1 && IsDatabaseValueSelected())
        {
            log.append(">>> Calculating...\n");
            int RowCount = ElementsTable.getRowCount();
            double[][] cord = new double[RowCount][3];
            double[] mass = new double[RowCount];
            float Xcm = 0.0f, Ycm = 0.0f, Zcm = 0.0f, totalMass = 0.0f;
            max = 0.0f;
            int[][] ex = PrepareArrayForMatrix();
          
            BADExtraction BADex = new BADExtraction();
        
            try {
                /*
                Extract2bonds,Extract3bonds and Extract4bonds functions are used to
                match bonds, angles and dihedrals respectively with the database...
                */
                StringBuffer Buffer2bonds = Extract2Bonds(BADex.getBondsList(ex)); 
                StringBuffer Buffer3bonds = Extract3Bonds(BADex.getAnglesList(ex),"CalGeo");
                StringBuffer Buffer4bonds = Extract4Bonds(BADex.getDihedralList(ex),"CalGeo");
                
                // Special condition diherals....
                if(RowCount == 4 || RowCount == 5)
                {
                    for(int i=0;i<RowCount;i++)
                    {
                        String linkColumn = ElementsTable.getValueAt(i, 4).toString();
                        String[] links = linkColumn.split(",");
                        if((links.length == 3 || links.length == 4) && Buffer4bonds == null)
                        {
                            Buffer4bonds = Buffer3bonds;
                        }
                    }
                }
                // Special condition  on bonds, angles and dihedrals...
                if(Buffer2bonds != null && (Buffer3bonds != null || RowCount <= 2) && (Buffer4bonds != null || RowCount <= 3))
                {
                    // Get coordinates from elements table and calculate totalmas from it...
                    for(int i = 0; i< RowCount;i++)
                    {
                        
                        //System.out.println("else called");
                        cord[i][0] = Double.parseDouble(ElementsTable.getValueAt(i, 1).toString().trim());
                        cord[i][1] = Double.parseDouble(ElementsTable.getValueAt(i, 2).toString().trim());
                        cord[i][2] = Double.parseDouble(ElementsTable.getValueAt(i, 3).toString().trim());

                        mass[i] = Double.parseDouble(ElementsTable.getValueAt(i, 9).toString().trim());
                        
                         totalMass += mass[i];
                    }
               // }
                
                    // Get center of mass....
                    double[] CoM= GetCenterOfMass(cord,mass,totalMass);
                    
                    // Subtract center of mass from original coordinates to use in geometry calculations...
                    for(int i = 0;i < cord.length;i++)
                    {
                            cord[i][0] = cord[i][0] - CoM[0];
                            cord[i][1] = cord[i][1] - CoM[1];
                            cord[i][2] = cord[i][2] - CoM[2];
                    }
                    
                    MolDyn mol = new MolDyn();
                    mol.SetBondPot(bondlist, bondlength, bondconst);//hchgfch
                    mol.SetAnglePot(anglelist, angle, angleconst);
                    mol.SetDihedPot(dihedlist, dihed, dihedconst);
                    
                   double[][] NewCord = mol.CalculateGeometry(cord, mass);
                   
                   for(int i=0;i < RowCount; i++)
                   {
                       float x = (float) NewCord[i][0];
                       float y = (float) NewCord[i][1];
                       float z = (float) NewCord[i][2];
                       
//                        max = 0.0f;
                        float Gr = Math.max(max, Math.abs(x));
                        Gr = Math.max(Gr, Math.abs(y));
                        max = Math.max(Gr, Math.abs(z));
                        
                       System.out.println("Atom ID : " + i + " [ "+x+" ] [ "+y+" ] [ "+z+" ]");
                       ElementsTable.setValueAt(x, i, 1);
                       ElementsTable.setValueAt(y, i, 2);
                       ElementsTable.setValueAt(z, i, 3);
                   }
                            max += 1.0;
                       System.out.println("Max value  : "+max);
                       log.append(">>>Max value by which coordinates are divided : ["+max+"]\n");
                       log.append(">>> values of atom on the scene graph...\n");
                       
                       for(int i=0;i < RowCount; i++)
                       {
                           // in order to resize structure to fit inside the scenegrap subtract all
                           // coordinates from the maximum...
                       float x = (float) (NewCord[i][0]/max);
                       float y = (float) (NewCord[i][1]/max);
                       float z = (float) (NewCord[i][2]/max);
                       
                       double massf = (float) mass[i];
                        
                       Xcm += x * massf;
                       Ycm += y * massf;
                       Zcm += z * massf;
                       
                       // new coordinates and values to draw structure on the scenegraph
                       // Note: it is different from the elements table parameters...
                       System.out.println("\t(" + i + ") [ "+x+" ] [ "+y+" ] [ "+z+" ]");
                       log.append("\t(" + i + ") [ "+x+" ] [ "+y+" ] [ "+z+" ]\n");
                       String links = ElementsTable.getValueAt(i, 4).toString();
                       CreateGraphics.changeVertex(i, x, y, z, links);
                       
                   }
                   
                    Xcm = Xcm / totalMass;
                    Ycm = Ycm / totalMass;
                    Zcm = Zcm / totalMass;
                   
                   CreateGraphics.setCenterOfMass(Xcm, Ycm, Zcm);
                   
                   BADError ec = new BADError(bondlist,anglelist,dihedlist,bondlength,angle,dihed,dihedconst,cord);
                   int[] error = ec.getGeometryError();
                   errorlbl.setText("::[ BE: "+error[0]+"% , AE: "+error[1]+"% , DE: "+error[2]+"% ]");
                   //errortxt.setEditable(true);
                   log.append(">>> Error calculated...\n");
                   log.append(">>> Geometry successfully calculated...\n");
                   
                   
                }
                else{
                    JOptionPane.showMessageDialog(this, "An Error Occured. Check your inputs. Database should also match for Angle and Dihedral ","Error",JOptionPane.ERROR_MESSAGE);
                    log.append("!!! Database match not found for Angles and Dihedrals\n");
                }
            } catch (IOException ex1) {
                Logger.getLogger(molecularModeling.class.getName()).log(Level.SEVERE, null, ex1);
                log.append("!!! Caught Exception while calculating geometry\n");
            }
               
            }
            else
                JOptionPane.showMessageDialog(this, "Error Occured. one or more Database value are missing / atoms are less than 2 ","Error",JOptionPane.ERROR_MESSAGE);
        
    }//GEN-LAST:event_CalculateGeometryActionPerformed
    
    /**
     * Calculate the center of mass of a structure.
     * 
     * @see molecularModeling#CalculateGeometryActionPerformed(java.awt.event.ActionEvent)
     * 
     * @param Coord Array of Coordinates of all atoms
     * @param mass Array of mass values of all atoms
     * @param totalMass Total mass value of all atoms
     * 
     * @return Array of center of mass
     */
    private double[] GetCenterOfMass(double[][] Coord,double[] mass,double totalMass)
    {
        int RowCount = Coord.length;
        double[] CoM = new double[3];
        double Xcm = 0.0f, Ycm = 0.0f, Zcm = 0.0f;
        
         for(int i = 0; i< RowCount;i++)
            {    
                    Xcm += Coord[i][0] * mass[i];
                    Ycm += Coord[i][1] * mass[i];
                    Zcm += Coord[i][2] * mass[i];
                            
            }
         
                   CoM[0] = Xcm / totalMass;
                   CoM[1] = Ycm / totalMass;
                   CoM[2] = Zcm / totalMass;
                   
                   return CoM;
    }
    
    /**
     * 
     * event handler: used to refresh database values inside the dropmenu of the elements table.
     * 
     * @see molecularModeling#getDBsNames(java.lang.String)
    */
    private void RefreshDBbtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_RefreshDBbtnActionPerformed
        try {
            DefaultTableModel dtm = (DefaultTableModel) ElementsTable.getModel();
                //dtm.fireTableDataChanged();
            int rows = ElementsTable.getRowCount();
            String[] DBNames = null;
            for(int i = 0; i < rows; i++)
            {
                if(ElementsTable.getValueAt(i, 5) == null)
                {
                    DBNames = getDBsNames("");
                }
                else
                {
                    String Value = ElementsTable.getValueAt(i, 5).toString();
                    DBNames = getDBsNames(Value);
                }
                
                Jbox.set(i, new JComboBox(DBNames));
                rowEditor.cancelCellEditing();
                rowEditor.setEditorAt(i, new DefaultCellEditor(Jbox.get(i)));
               
            }   
                log.append(">>> Databse field successfully refreshed...\n");
        } catch (IOException ex) {
            System.out.println(ex.getStackTrace());
            log.append("!!! Caught Exception while refreshing the database field\n");
        }
    }//GEN-LAST:event_RefreshDBbtnActionPerformed
    
    /**
     * 
     * event handler: used to restart the program
     * using command line method for both operating systems.
     * 
    */
    private void RestartbtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_RestartbtnActionPerformed
        String OperatingSystem = System.getProperty("os.name");
                String command = String.format("java -jar \""+UserDirectory+""+File.separator+"ms2potmod.jar\"");
                //JOptionPane.showMessageDialog(null, command);
                ProcessBuilder builder = null;
                if (OperatingSystem.toLowerCase().contains("windows")) {
                        
			builder = new ProcessBuilder("cmd", "/c", command);
		} else {
			// In case of Linux/Ubuntu run command using /bin/bash
			builder = new ProcessBuilder("/bin/bash", "-c", command);
		}
        
        try {
            builder.start();
        } catch (IOException ex) {
            Logger.getLogger(molecularModeling.class.getName()).log(Level.SEVERE, null, ex);
        }
        System.exit(0);
    }//GEN-LAST:event_RestartbtnActionPerformed
    
    /**
     * 
     * open button event handler
     * Used to Show open file dialogbox
     * Calls read file method.
     * 
     * @see molecularModeling#ReadPMFile(java.lang.String)
    */
    private void OpenBtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_OpenBtnActionPerformed
//        ISGeom = true;
        JFileChooser fileChooser = new JFileChooser();
        FileFilter fileFilter = new FileNameExtensionFilter("PM Files","pm");
        fileChooser.setFileFilter(fileFilter);
        if(OpenSelectedDir != null)
                fileChooser.setCurrentDirectory(OpenSelectedDir);
                //fileChooser.setCurrentDirectory(new File(System.getProperty("user.home")));
        
        int result = fileChooser.showOpenDialog(this);
        if (result == JFileChooser.APPROVE_OPTION) {
                File selectedFile = fileChooser.getSelectedFile();
                OpenSelectedDir = selectedFile.getParentFile();
                log.append(">>> Open file Selected:     "+selectedFile + "\n");
                
               // System.out.println("Selected file: " + selectedFile.getAbsolutePath());
            try {
                log.append(">>> reading file...\n");
                ReadPMFile(selectedFile.getAbsolutePath());
                log.append(">>> Done...\n");
            } catch (FileNotFoundException ex) {
                Logger.getLogger(molecularModeling.class.getName()).log(Level.SEVERE, null, ex);
                log.append("!!! Caught Exception while reading selected PM file\n");
            }
                
}
    }//GEN-LAST:event_OpenBtnActionPerformed

    private void formComponentResized(java.awt.event.ComponentEvent evt) {//GEN-FIRST:event_formComponentResized
       
       
    }//GEN-LAST:event_formComponentResized
    
    /**
     * 
     * help button event handler: used to show the help file.
     * 
    */
    private void helpbtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_helpbtnActionPerformed
       
    try{
         String filePath = UserDirectory+""+File.separator+"help"+File.separator+"ms2potmod_help_format.xlsx";
         Desktop.getDesktop().open(new File(filePath));
    }catch (Exception e){
      JOptionPane.showMessageDialog(null,"Exception caught ="+e.getMessage());
      log.append("!!! Caught Exception while fetching help file\n");
    }
    }//GEN-LAST:event_helpbtnActionPerformed
    
    /**
     * 
     * LJ-Site button event handler
     * Calls ljSite function and Show LJsite window.
     * 
     * @see molecularModeling#IsDatabaseValueSelected() 
    */
    private void ljsitebtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ljsitebtnActionPerformed
        if(IsDatabaseValueSelected())
        {
        String[][] tContent = new String[ElementsTable.getRowCount()][3];
        for(int i=0;i<ElementsTable.getRowCount();i++)
        {
            tContent[i][0] = ElementsTable.getValueAt(i, 0).toString();
            tContent[i][1] = ElementsTable.getValueAt(i, 6).toString();
            tContent[i][2] =  String.valueOf(ElementsTable.getValueAt(i, 4).toString().trim().split(",").length);
            
        }
        
       
        ljSite.GenerateLJTypeTable(tContent, this);
        ljSite.setVisible(true);
        }
        else{
            JOptionPane.showMessageDialog(this, "Select databases first!", "Info", JOptionPane.INFORMATION_MESSAGE);
        }
    }//GEN-LAST:event_ljsitebtnActionPerformed

    /**
     * 
     * add DB entry button event handler:
     * Show the DB entry window.
     * 
     * @see DBEntry
    */
    private void addDBentryActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_addDBentryActionPerformed
        this.setEnabled(false);
        DBEntry dbEntry = new DBEntry(UserDirectory, helpMap, this);
        dbEntry.setVisible(true);
    }//GEN-LAST:event_addDBentryActionPerformed
    
    /**
     * 
     * Set the frame size according to the width and height of the screen.
     * 
     */
    private void SetFrame()
    {
        
        Toolkit toolkit = Toolkit.getDefaultToolkit(); 
        int width= (int)(toolkit.getScreenSize().width * 2/3); 
        int height= (int)(toolkit.getScreenSize().height * 2/3); 
        //this.setSize(1000, 5000000); 
        this.setPreferredSize(new Dimension(width,height));
        this.pack();
        
        log.append("\n>>> Application frame adjusted.\n");
        
    }
    
    /**
     * 
     * Clear elements table
     * clear scenegraph
     * read provided file - all the required parameters from the file
     * Draw new structures.
     * 
     * @see CreateGraphics
     * @see molecularModeling#ReadDbNamesForReadFile(java.lang.String[][], int)
     * 
     * @param File Path of the file to read
    */
    private void ReadPMFile(String File) throws FileNotFoundException
    {
     try{
         DefaultTableModel dtm = (DefaultTableModel) ElementsTable.getModel();
                    for(int i = dtm.getRowCount() - 1 ; i >= 0  ; i--)
                    {
                        dtm.removeRow(i);
                    }
          Jbox.clear();
         CreateGraphics.CleanSceneGraph();
        FileReader fr = new FileReader(File);
        BufferedReader reader = new BufferedReader(fr);
        String line = null;
        int counter = 0,index = 0, IDcounter = 0;
        float maxi = 0.0f;
        boolean IsOutofRange = false;
        while ((line = reader.readLine()) != null) {
            if(line.contains("NSites"))
            {
                String[] Splitvalue = line.split("=");
                IDcounter = Integer.parseInt(Splitvalue[1].trim().toString());
               // System.out.print("counter : "+IDcounter);
                
                break;
            }
        }
        
        String[][] param = new String[IDcounter + 1][7];
        int[] apprd = new int[IDcounter + 1];
        float[][] CvCord = new float[IDcounter + 1][3];
        
        while ((line = reader.readLine()) != null) {
        //System.out.println(line);
            
            // read SiteID to get all the coordinates for the single ID...
            if(line.contains("SiteID"))
            {
                //String[] values = new String[7];
                String[] Splitvalue = line.split("=");
                int ID = Integer.parseInt(Splitvalue[1].trim().toString());
               
                for(int i=0; i < 6; i++)
                {
                line = reader.readLine();
                Splitvalue = line.split("=");
               // values[i] = Splitvalue[1];
                param[ID][i] = Splitvalue[1];
                    if(i < 3)
                    {
                        float val = Float.parseFloat(param[ID][i]);
                        maxi = Math.max(maxi, Math.abs(val));
                        CvCord[ID][i] = val;
                        
                        if(!(val > -1 && val < 1))
                         {
                             IsOutofRange = true;
                         }
                    }
                }
            }
           
            /*
            - only require bond information in order to recreate the structure
            */
            if(line.startsWith("Bond"))
            {
                if(IsOutofRange)
                {
                    maxi = maxi + 1.0f;
                    System.out.println("Max Value : "+maxi);
                    log.append(">>> Max value selected : ["+maxi+"]");
                    log.append(">>> values of atoms on scene graph are... \n");
                    for(int i= 1; i <=IDcounter;i++)
                    {
                        CvCord[i][0] = Float.parseFloat(param[i][0].trim())/maxi;
                        CvCord[i][1] = Float.parseFloat(param[i][1].trim())/maxi;
                        CvCord[i][2] = Float.parseFloat(param[i][2].trim())/maxi;
                        System.out.println("\t(" + i + ") [ "+CvCord[i][0]+" ] [ "+CvCord[i][1]+" ] [ "+CvCord[i][2]+" ]");
                        log.append("\t(" + i + ") [ "+CvCord[i][0]+" ] [ "+CvCord[i][1]+" ] [ "+CvCord[i][2]+" ]");
                      //  JOptionPane.showMessageDialog(null, i);
                    }
                    
                    IsOutofRange = false;
                }
                String[] FS = line.split("=");
                String[] SS = FS[1].split("\\s+");

                int a = Integer.parseInt(SS[1].trim());
                int b = Integer.parseInt(SS[2].trim());
                             float ax = CvCord[a][0];
                             float ay = CvCord[a][1];
                             float az = CvCord[a][2];
                
                 // First atom in a bond:
                // if not appeared before then add to the scenegraph and elements table
                // if appeared then add to the links coloumn of the elements table
                if(!Search(apprd,a))
                {
                             apprd[counter] = a;
                             dtm.addRow(new Object[]{a-1,param[a][0],param[a][1],param[a][2],b-1,null,null,param[a][3],param[a][4],param[a][5]});
                             //dtm.setValueAt("CH3-OPLS-mado", dtm.getRowCount() - 1, 6);
                             CreateGraphics.RedrawSphere(a-1, ax, ay, az);
                             
                }
                else
                {
                    String link = dtm.getValueAt(a-1, 4).toString()+","+(b-1);
                    dtm.setValueAt(link, a-1, 4);
                    counter--;
                }
                             float bx = CvCord[b][0];
                             float by = CvCord[b][1];
                             float bz = CvCord[b][2];
                             
               // Second atom in a bond:
                // if not appeared before then add to the scenegraph and elements table
                // if appeared then add to the links coloumn of the elements table
               if(!Search(apprd,b))
                {
                             apprd[counter + 1] = b;
                             dtm.addRow(new Object[]{b-1,param[b][0],param[b][1],param[b][2],a-1,null,null,param[b][3],param[b][4],param[b][5]});
                            // dtm.setValueAt("CH3-OPLS-mado", dtm.getRowCount() - 1, 6);
                             CreateGraphics.RedrawSphere(b-1, bx, by, bz);
                }
                else
                {
                    String link = dtm.getValueAt(b-1, 4).toString()+","+(a-1);
                    dtm.setValueAt(link, b-1, 4);
                }
                
               // then connects these two bonded atoms together on the scenegraph by
               // calling createGraphics class...
                Point3f CuRrPoint = new Point3f(bx, by, bz);
                Point3f PrEvpoint = new Point3f(ax,ay,az);
                CreateGraphics.DrawLinE(PrEvpoint, CuRrPoint, index);
                CreateGraphics.rEstScene();
                CreateGraphics.SetLineArray(index, a-1, b-1);
                counter+=2; index++;
            }
        }
        reader.close();
        fr.close();
        
        // Sets the Discription of the atoms on the sceneGraph by
        // providing links information which is important!
        // it shows which atom is connected to which...
        for(int i=0;i<dtm.getRowCount();i++)
        {
            String[] link = dtm.getValueAt(i, 4).toString().split(",");
            if(link.length < 2)
                CreateGraphics.SetDiscription(i, String.valueOf(link.length - 1));
            else
               CreateGraphics.SetDiscription(i, String.valueOf(link.length)); 
        }
        
        ReadDbNamesForReadFile(param, IDcounter);
        
        }
        catch(Exception ex)
                {
                    ex.printStackTrace();
                }
        
    }
    
    /**
     * 
     * Read following parameters from the Database file:
     * Site_Type,
     * Epsilon,
     * Sigma,
     * Mass.
     * 
     * @see molecularModeling#ReadPMFile(java.lang.String) 
     * 
     * @param ParamSEM Array of Parameters of the main table
     * @param counter
     * 
    */
    private void ReadDbNamesForReadFile(String[][] ParamSEM, int counter) throws IOException
    {
        
        try {
            for(int i=1; i<=counter; i++)
            {
            String DBName = null;
            BufferedReader in = new BufferedReader(new FileReader(UserDirectory+""+File.separator+"Database"+File.separator+"potparam.txt"));
            while(in.ready())
            {
                String line = in.readLine().trim();
                if(line.startsWith("site_type"))
                {
                   String[] DBn = line.split("=");
                   String DName = DBn[1].replace("\\par", "");
                   DName = DName.trim();
                   
                   line = in.readLine().trim();
                   DBn = line.split("=");
                   String Sigma = DBn[1].replace("\\par", "");
                   Sigma = Sigma.trim();
                   
                   line = in.readLine().trim();
                   DBn  = line.split("=");
                   String Epsilon = DBn[1].replace("\\par", "");
                   Epsilon = Epsilon.trim();
                   
                   line = in.readLine().trim();
                   DBn  = line.split("=");
                   String Mass = DBn[1].replace("\\par", "");
                   Mass = Mass.trim();
                  
                   if(Sigma.equals(ParamSEM[i][3].trim()) && Epsilon.equals(ParamSEM[i][4].trim()) && Mass.equals(ParamSEM[i][5].trim()))
                   {
                    DBName = DName;
                    break;
                   }
            }
            }
            in.close();
                if(DBName == null) {
                    JOptionPane.showMessageDialog(null, "Match not found in the Database file for row: " + (i - 1));
                    ElementsTable.setValueAt("", i-1, 7);
                    ElementsTable.setValueAt("", i-1, 8);
                    ElementsTable.setValueAt("", i-1, 9);
                }
            // add Site_Type values to the dropdown menus....
            JComboBox Jcomb = new JComboBox(DBvalues);
            Jcomb.setSelectedItem(DBName);
            ElementsTable.setValueAt(DBName, i-1, 6);
            Jbox.add(Jcomb);
            rowEditor.setEditorAt(i-1, new DefaultCellEditor(Jbox.get(i-1)));
            
            // add type search values inside the dropdown menu...
            String[] types = {"CH","CH2","CH3","C","O"};
            JComboBox jb = new JComboBox(types);
            jb.setEditable(true);
            TJbox.add(jb);
            TrowEditor.setEditorAt(i-1, new DefaultCellEditor(TJbox.get(i-1)));
            }
          
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        }
    }

    /**
     *
     * Check if the value exist or not
     * 
     * @param arr Array of values
     * @param targetValue Vlaue to check
     * 
     * @return True - False
     */
    public boolean Search(int[] arr, int targetValue) {
	for(int s: arr){
		if(s == targetValue)
			return true;
	}
	return false;
}
    /**
     * 
     * Start scenegrpah window frame.
     * 
    */
    private void StartSceneGraphWindow()
    {
        BasicInternalFrameUI ui = (BasicInternalFrameUI) CanvasInternalFrame.getUI();
        Container north = (Container)ui.getNorthPane();
        north.remove(0);
        north.validate();
        north.repaint();
        CanvasInternalFrame.setLayout(new BorderLayout());
        CanvasInternalFrame.add(new CreateGraphics().CreateGraphics(this));
        this.pack();
        
        log.append(">>> Scene graph started...\n");
    }
    
    /**
     * 
     * Prepare array for calculating the Matrix.
     * 
    */
    private int[][] PrepareArrayForMatrix()
    {
        int[][] ex = new int[ElementsTable.getRowCount()][];
        for(int i=0;i<ElementsTable.getRowCount();i++)
        {
            String value = ElementsTable.getValueAt(i, 4).toString();
            String[] bond = value.split(",");
            int id = Integer.parseInt(ElementsTable.getValueAt(i, 0).toString());
            ex[i] = new int[bond.length + 1];
            int j=0;
            while(j<bond.length + 1)
            {
                if(j==0)
                {
                    ex[i][j] = id + 1;
                   // System.out.println(ex[i][j] +" ");
                }
                else
                {
                    int ind = Integer.parseInt(bond[j - 1]);
                    ex[i][j] = ind + 1;
                   // System.out.println(ex[i][j] +" ");
                }
 
                j++;
                
            }
         // System.out.println("\n");
        }
        return ex;
    }
    
    /**
     * 
     * used to write the structures in the form of .par files.
     * 
     * @see molecularModeling#WriteDBPropertiesInPmFile(java.io.BufferedWriter) 
     * 
     * @param file Reference of the file to be saved
     * 
    */
    private void SavepmFile(File file)
    {
        try {
			if (!file.exists()) {
                            
				file.createNewFile();
			}
 
			FileWriter fw = new FileWriter(file.getAbsoluteFile(),false);
			BufferedWriter bw = new BufferedWriter(fw);
                        ArrayList<String> Charge = ljSite.GetLJType();
                        ArrayList<Integer> benzenelist = new ArrayList<Integer>();
                        
                        int count = 0;
                        if(Charge != null) {
                            for(int i = 0; i < Charge.size(); i++) {
                                String[] index = Charge.get(i).split("-");
                                if(index[1].equals("OH")) {
                                        count++;
                                }
                                else
                                if(index[1].equals("benzene")){
                                    int id = Integer.parseInt(index[0].trim());
                                    benzenelist.add(id);
                                }
                                
                            }
                        }
                                    
                        bw.append("NSiteTypes =  1\n" +
                                    "\n" +
                                    "SiteType   =  LJ126\n" +
                                    "NSites     =  "+(ElementsTable.getRowCount() + count)+"\n\n");
                        int siteID = 0;
                        for(int i=0;i<this.ElementsTable.getRowCount();i++)
                        {
                            
                            //siteID = Integer.parseInt(ElementsTable.getValueAt(i, 0).toString());
                            if(!benzenelist.contains(i)) {
                                bw.append("SiteID     =  " + (++siteID) + "\n" +
                                            "x          =  " + ElementsTable.getValueAt(i, 1) + "\n" +
                                            "y          =  " + ElementsTable.getValueAt(i, 2) + "\n" +
                                            "z          =  " + ElementsTable.getValueAt(i, 3) + "\n" +
                                            "sigma      = " + ElementsTable.getValueAt(i, 7) + "\n" +
                                            "epsilon    = " + ElementsTable.getValueAt(i, 8) + "\n" +
                                            "mass       = "+ ElementsTable.getValueAt(i, 9) + "\n\n");
                            }
                            
                        }
                        
                    
                        //ArrayList<Integer> _constraints1;
                        //ArrayList<Integer> _constraints2;
                        ArrayList<Integer> temp = new ArrayList<Integer>();
                        ArrayList<Integer> temp_x = new ArrayList<Integer>();
                        ArrayList<ArrayList<Integer>> constraints = new ArrayList<ArrayList<Integer>>();
                        ArrayList<ArrayList<Integer>> Bond = new ArrayList<ArrayList<Integer>>();
                        ArrayList<ArrayList<Integer>> Angle = new ArrayList<ArrayList<Integer>>();
                        ArrayList<ArrayList<Integer>> Dihedral = new ArrayList<ArrayList<Integer>>();
                    if(Charge != null){
                        int limit = siteID;
                        Vector3f[] O_cord = new Vector3f[Charge.size()];
                        Vector3f[] H_cord = new Vector3f[Charge.size()];
                        boolean[] branched = new boolean[Charge.size()];
                        boolean isOHexist = false;
                        int ch = 0;
                        
                            for(int i=0,k=0; i < Charge.size() ;i++)
                            {
                                String[] index = Charge.get(i).split("-");
                                k= Integer.parseInt(index[0]);
                                String links = ElementsTable.getValueAt(k, 4).toString();
                                String[] Clinks = links.split(",");
                                if(index[1].equals("charge")) {
                                    ch++;
                                }
                                if(index[1].equals("benzene") && Clinks.length < 2) {
                                   Vector3f KVec = getCoordinatesofAnID(k);
                                   limit = benzeneCHtemplate(KVec, bw, limit, false);
                                }
                             if(index[1].equals("OH")) {
                                 isOHexist = true;
                                
                                ///////////////////////////////////
                                
                                if(Clinks.length < 3) {
                                    float x_pos = Float.parseFloat(ElementsTable.getValueAt(k, 1).toString())/10.0f;
                                    float y_pos = Float.parseFloat(ElementsTable.getValueAt(k, 2).toString())/10.0f;
                                    float z_pos = Float.parseFloat(ElementsTable.getValueAt(k, 3).toString())/10.0f;
                                    if(Clinks.length == 0) {
                                        O_cord[i] = new Vector3f(x_pos + 0.15f, y_pos, z_pos);
                                        H_cord[i] = new Vector3f(O_cord[i].x + 0.15f, O_cord[i].y + 0.15f, O_cord[i].z);
                                    }
                                    else if(Clinks.length == 1) {
                                        int neigh = Integer.parseInt(Clinks[0]);
                                        Vector3f neighVec = getCoordinatesofAnID(neigh);
                                        Vector3f KVec = getCoordinatesofAnID(k);
                                        temp_x.add(neigh);
                                        x_pos = neighVec.x; y_pos = neighVec.y; z_pos = neighVec.z;
                                        if(neighVec.x == KVec.x && neighVec.y != KVec.y)
                                            y_pos = -1f * neighVec.y;
                                        else
                                          x_pos = -1f * neighVec.x;
                                       /* if(neighVec.y != KVec.y)
                                            y_pos = -1f * neighVec.y;
                                        if(neighVec.z != KVec.z)
                                            z_pos = -1f * neighVec.z;*/
                                        
                                        O_cord[i] = new Vector3f(x_pos, y_pos, z_pos);
                                        H_cord[i] = new Vector3f(O_cord[i].x + 0.15f, O_cord[i].y + 0.15f, O_cord[i].z);
                                        
                                    }
                                    else
                                    if(Clinks.length == 2) {
                                        int neighR = Integer.parseInt(Clinks[0]);
                                        int neighL = Integer.parseInt(Clinks[1]);
                                        
                                        Vector3f neighRvec = getCoordinatesofAnID(neighR);
                                        Vector3f neighLvec = getCoordinatesofAnID(neighL);
                                        
                                        x_pos = (neighRvec.y * neighLvec.z) - (neighRvec.z * neighLvec.y);
                                        y_pos = (neighRvec.z * neighLvec.x) - (neighRvec.x * neighLvec.z);
                                        z_pos = (neighRvec.x * neighLvec.y) - (neighRvec.y * neighLvec.x);
                                        
                                        //Vector3f newVec = new Vector3f(neighRvec.x + neighLvec.x, neighRvec.y + neighLvec.y, neighRvec.z + neighLvec.z);
                                        //float mag =  (float) Math.sqrt((newVec.x * newVec.x) + (newVec.y * newVec.y) + (newVec.z * newVec.z));
                                        
                                        //x_pos = newVec.x/mag; y_pos = newVec.y/mag; z_pos = newVec.y/mag;
                                        
                                        O_cord[i] = new Vector3f(x_pos, y_pos, z_pos);
                                        H_cord[i] = new Vector3f(O_cord[i].x + 0.15f, O_cord[i].y + 0.15f, O_cord[i].z);

                                        temp_x.add(neighR);
                                        branched[i] = true;
                                    }
                                }
                                //////////////////////////////////
                                bw.append("#O\n");
                                bw.append("SiteID     =  " + (++limit) + "\n" +
                                            "x          =  " + (O_cord[i].x) + "\n" +
                                            "y          =  " + (O_cord[i].y)+ "\n" +
                                            "z          =  " + (O_cord[i].z)+ "\n" +
                                            "sigma      =   3.149559\n" +
                                            "epsilon    =  85.053449\n" +
                                            "mass       =  15.999\n\n" );
                                temp.add(limit);
                              }
                            }
                            if(ch > 0 || isOHexist) {
                            bw.append("SiteType   =   CHARGE\n");
                            bw.append("NSites     =   "+ (((Charge.size() - ch) * 3 ) + ch)+"\n\n");
                            }
                            
                            for(int i = 0, k = 0; i < Charge.size() ; i++)
                            {
                                
                                //_constraints1 = new ArrayList<Integer>();
                                //_constraints2 = new ArrayList<Integer>();
                                ArrayList<Integer> bo = new ArrayList<Integer>();
                                ArrayList<Integer> an = new ArrayList<Integer>();
                                ArrayList<Integer> di = new ArrayList<Integer>();
                                
                                String[] index = Charge.get(i).split("-");
                                k= Integer.parseInt(index[0]);
                                String links = ElementsTable.getValueAt(k, 4).toString();
                                String[] Clinks = links.split(",");
                                
                                if(index[1].equals("charge")) {
                                    
                                    bw.append("SiteID     =  " + (++limit) + "\n" +
                                            "x          =  " + ElementsTable.getValueAt(k, 1) + "\n" +
                                            "y          =  " + ElementsTable.getValueAt(k, 2) + "\n" +
                                            "z          =  " + ElementsTable.getValueAt(k, 3) + "\n" +
                                            "Charge     =  -1.0000\n" +
                                            "mass       =  0.0\n" +
                                            "shielding  =  0.71092\n\n" );
                                    
                                }
                              if(index[1].equals("OH")) {
                                //k= Integer.parseInt(index[0]);

                                String db = ElementsTable.getValueAt(k, 6).toString();
                                String namen = db.substring(1,db.indexOf("-"));
                                bw.append("#" + namen + "\n");
                                bw.append("SiteID     =  " + (++limit) + "\n" +
                                            "x          =  " + ElementsTable.getValueAt(k, 1) + "\n" +
                                            "y          =  " + ElementsTable.getValueAt(k, 2) + "\n" +
                                            "z          =  " + ElementsTable.getValueAt(k, 3) + "\n" +
                                            "Charge     =  +0.2556\n" +
                                            "mass       =  0.0\n" +
                                            "shielding  =  1.38448\n\n" );
                               // _constraints1.add(k + 1);
                               // _constraints1.add(limit);
                                bw.append("#" + index[1].charAt(0) + "\n");
                                bw.append("SiteID     =  " + (++limit) + "\n" +
                                            "x          =  " + (O_cord[i].x) + "\n" +
                                            "y          =  " + (O_cord[i].y) + "\n" +
                                            "z          =  " + (O_cord[i].z) + "\n" +
                                            "Charge     =  -0.697107\n" +
                                            "mass       =  0.0\n" +
                                            "shielding  =   1.25982\n\n" );
                               // _constraints2.add(temp.get(i));
                               // _constraints2.add(limit);
                                bw.append("#" + index[1].charAt(1) + "\n");
                                bw.append("SiteID     =  " + (++limit) + "\n" +
                                            "x          =  " + (H_cord[i].x) + "\n" +
                                            "y          =  " + (H_cord[i].y) + "\n" +
                                            "z          =  " + (H_cord[i].z) + "\n" +
                                            "Charge     =  +0.441507\n" +
                                            "mass       =  1.008\n" +
                                            "shielding  =  0.0\n\n" );
                               // _constraints2.add(limit);
                                //Bond
                                ///////bo.add(k + 1);
                                ///////bo.add(temp.get(i));
                                ///////Bond.add(bo);
                                //
                                if(!branched[i]) {
                                // Angle
                                    an.add(k + 1);
                                    an.add(temp.get(i));
                                    an.add(limit);
                                    Angle.add(an);
                                //
                                // dihedral
                                    di.add(temp_x.get(i) + 1);
                                    di.add(k + 1);
                                    di.add(temp.get(i));
                                    di.add(limit);
                                    Dihedral.add(di);
                                //
                                }
                                else
                                {
                                 //Angle 1
                                    an.add(k + 1);
                                    an.add(temp.get(i));
                                    an.add(limit);
                                    Angle.add(an);
                                 //
                                 //Angle 2
                                    an = new ArrayList<Integer>();
                                    an.add(temp_x.get(i) + 1);
                                    an.add(k + 1);
                                    an.add(temp.get(i));
                                    Angle.add(an);
                                    
                                    // dihedral
                                    di.add(temp_x.get(i) + 1);
                                    di.add(k + 1);
                                    di.add(temp.get(i));
                                    di.add(limit);
                                    Dihedral.add(di);
                                    
                                }
                               // constraints.add(_constraints1);
                               // constraints.add(_constraints2);
                                
                                
                            }
                            
                        if(index[1].equals("benzene") && Clinks.length < 2) {
                            bw.append("\nSiteType   =   Quadrupole\n");
                            bw.append("NSites     =   6\n\n");
                             Vector3f KVec = getCoordinatesofAnID(k);
                             limit = benzeneCHtemplate(KVec, bw, limit, true);
                          }
                            
                          }
                            ArrayList<Integer> _constraints1 = null;
                            for(int j = 0; j < Charge.size(); j++) {
                                _constraints1 = new ArrayList<Integer>();
                                for(int k = 1; k <= limit; k++) {
                                   _constraints1.add(k);
                                }
                                
                                constraints.add(_constraints1);
                            }
                                
                        }
                        
                     boolean IsOut = WriteDBPropertiesInPmFile(bw, constraints, Bond, Angle, Dihedral);
                     bw.close();
                     if(!IsOut)
                     {
                         file.delete();
                         JOptionPane.showMessageDialog(this, "Error Occured. Can not write file. Check your inputs and Try again","Error",JOptionPane.ERROR_MESSAGE);
                         log.append(">>> Error Occured. Can not write file. Check your inputs and Try again...\n");
                     }
                     else
                     {
                        JOptionPane.showMessageDialog(this, "The file has successfully saved","Saved",JOptionPane.INFORMATION_MESSAGE);
			log.append(">>> File save done...\n");
                     }
		} catch (IOException e) {
			e.printStackTrace();
                       // JOptionPane.showMessageDialog(this, "Caught an Exception. Please try again","Exception",JOptionPane.ERROR_MESSAGE);
                        log.append(">>> Caught an Exception. Please try again...\n");
		}
    }
    
    Vector3f getCoordinatesofAnID(int k) {
        float x_pos = Float.parseFloat(ElementsTable.getValueAt(k, 1).toString());
        float y_pos = Float.parseFloat(ElementsTable.getValueAt(k, 2).toString());
        float z_pos = Float.parseFloat(ElementsTable.getValueAt(k, 3).toString());
        
        return new Vector3f(x_pos, y_pos, z_pos);
    }
    
    int benzeneCHtemplate(Vector3f selected, BufferedWriter bw, int inc, boolean quad) throws IOException {
        float x, y, z;
        //extend selected vector
        Vector3f temp = new Vector3f(selected.x * 2, selected.y * 2, selected.z * 2);
        
        //claculate magnitude of selected vector and extended vector (temp)
        double mag_temp = Math.sqrt((temp.x * temp.x) + (temp.y * temp.y) + (temp.z * temp.z));
        double mag_sel = Math.sqrt((selected.x * selected.x) + (selected.y * selected.y) + (selected.z * selected.z));
        
        //calculate division of magnitudes of selected and extended vectors
        float temp_sel_div = (float)(mag_temp/mag_sel);
        
        //claculate new vector from selected vector
        Vector3f temp_now = new Vector3f(selected.x * temp_sel_div, selected.y * temp_sel_div, selected.z * temp_sel_div);
        // orignal vector 1 to 6 from template
        Vector3f temp_orig = new Vector3f(-2.4455f, 1.4119f, 0);
        
        //calculate transform vector by substracting orignal vector (1 to 6 | temp_orig) from new selected vector (temp_now)
        Vector3f temp_trans = new Vector3f(temp_now.x - temp_orig.x, temp_now.y - temp_orig.y, temp_now.z - temp_orig.z);
        
        //calculate magnitude of transform vector
        float temp_trans_mag = (float) Math.sqrt((temp_trans.x * temp_trans.x) + (temp_trans.y * temp_trans.y) + (temp_trans.z * temp_trans.z));
        
        // calculate the unit vector of transform vector
        Vector3f temp_trans_u = new Vector3f(temp_trans.x/temp_trans_mag, temp_trans.y/temp_trans_mag, temp_trans.z/temp_trans_mag);
        
        bw.append("SiteID     =  " + (++inc) + "\n" +
                "x          =  " + selected.x + "\n" +
                "y          =  " + selected.y + "\n" +
                "z          =  " + selected.z + "\n");
        if(!quad) {
        bw.append("sigma      =   3.446\n" +
                "epsilon    =   70.019\n" +
                "mass       =    13.019\n\n" );
        } else {
                bw.append("theta      =   0\n" +
                          "phi        =   0\n" +
                          "quadrupole =  -1.0435\n" +
                          "mass       =   0\n" +
                          "shielding  =   1\n\n");
        }
        
        // orignal 1-2 vector from the tempate
        Vector3f orig12 =new Vector3f(-1.6303f, 0, 0);
        x = (orig12.y * temp_trans_u.z) - (orig12.z * temp_trans_u.y);
        y = (orig12.z * temp_trans_u.x) - (orig12.x * temp_trans_u.z);
        z = (orig12.x * temp_trans_u.y) - (orig12.y * temp_trans_u.x);
        
        bw.append("SiteID     =  " + (++inc) + "\n" +
                "x          =  " + x + "\n" +
                "y          =  " + y + "\n" +
                "z          =  " + z + "\n");
        if(!quad) {
        bw.append("sigma      =   3.446\n" +
                "epsilon    =   70.019\n" +
                "mass       =    13.019\n\n" );
        } else {
                bw.append("theta      =   0\n" +
                          "phi        =   0\n" +
                          "quadrupole =  -1.0435\n" +
                          "mass       =   0\n" +
                          "shielding  =   1\n\n");
        }
        
        // orignal 1-3 vector from the tempate
        Vector3f orig13 =new Vector3f(-2.4455f, 1.4119f, 0);
        x = (orig13.y * temp_trans_u.z) - (orig13.z * temp_trans_u.y);
        y = (orig13.z * temp_trans_u.x) - (orig13.x * temp_trans_u.z);
        z = (orig13.x * temp_trans_u.y) - (orig13.y * temp_trans_u.x);
        
        bw.append("SiteID     =  " + (++inc) + "\n" +
                "x          =  " + x + "\n" +
                "y          =  " + y + "\n" +
                "z          =  " + z + "\n");
        if(!quad) {
        bw.append("sigma      =   3.446\n" +
                "epsilon    =   70.019\n" +
                "mass       =    13.019\n\n" );
        } else {
                bw.append("theta      =   0\n" +
                          "phi        =   0\n" +
                          "quadrupole =  -1.0435\n" +
                          "mass       =   0\n" +
                          "shielding  =   1\n\n");
        }
        
       // orignal 1-4 vector from the tempate
        Vector3f orig14 =new Vector3f(-1.6303f, 2.8238f, 0);
        x = (orig14.y * temp_trans_u.z) - (orig14.z * temp_trans_u.y);
        y = (orig14.z * temp_trans_u.x) - (orig14.x * temp_trans_u.z);
        z = (orig14.x * temp_trans_u.y) - (orig14.y * temp_trans_u.x);
        
        bw.append("SiteID     =  " + (++inc) + "\n" +
                "x          =  " + x + "\n" +
                "y          =  " + y + "\n" +
                "z          =  " + z + "\n");
        if(!quad) {
        bw.append("sigma      =   3.446\n" +
                "epsilon    =   70.019\n" +
                "mass       =    13.019\n\n" );
        } else {
                bw.append("theta      =   0\n" +
                          "phi        =   0\n" +
                          "quadrupole =  -1.0435\n" +
                          "mass       =   0\n" +
                          "shielding  =   1\n\n");
        }
        
        // orignal 1-5 vector from the tempate
        Vector3f orig15 =new Vector3f(0, 2.8238f, 0);
        x = (orig15.y * temp_trans_u.z) - (orig15.z * temp_trans_u.y);
        y = (orig15.z * temp_trans_u.x) - (orig15.x * temp_trans_u.z);
        z = (orig15.x * temp_trans_u.y) - (orig15.y * temp_trans_u.x);
        
        bw.append("SiteID     =  " + (++inc) + "\n" +
                "x          =  " + x + "\n" +
                "y          =  " + y + "\n" +
                "z          =  " + z + "\n");
        if(!quad) {
        bw.append("sigma      =   3.446\n" +
                "epsilon    =   70.019\n" +
                "mass       =    13.019\n\n" );
        } else {
                bw.append("theta      =   0\n" +
                          "phi        =   0\n" +
                          "quadrupole =  -1.0435\n" +
                          "mass       =   0\n" +
                          "shielding  =   1\n\n");
        }
        
        // orignal 1-6 vector from the tempate
        Vector3f orig16 =new Vector3f(0.8152f, 1.4119f, 0);
        x = (orig16.y * temp_trans_u.z) - (orig16.z * temp_trans_u.y);
        y = (orig16.z * temp_trans_u.x) - (orig16.x * temp_trans_u.z);
        z = (orig16.x * temp_trans_u.y) - (orig16.y * temp_trans_u.x);
        
        bw.append("SiteID     =  " + (++inc) + "\n" +
                "x          =  " + x + "\n" +
                "y          =  " + y + "\n" +
                "z          =  " + z + "\n");
        if(!quad) {
        bw.append("sigma      =   3.446\n" +
                "epsilon    =   70.019\n" +
                "mass       =    13.019\n\n" );
        } else {
                bw.append("theta      =   0\n" +
                          "phi        =   0\n" +
                          "quadrupole =  -1.0435\n" +
                          "mass       =   0\n" +
                          "shielding  =   1\n\n");
        }
        
        return inc;
        
    }
    /**
     * 
     * Used to get database values from elements table for provided bonds
     * match the provided list of bonds with combination inside the database file
     * Prepare a string buffer to write bonds parameters in the .par file.
     * 
     * @see molecularModeling#CalculateGeometryActionPerformed(java.awt.event.ActionEvent)
     * @see molecularModeling#WriteDBPropertiesInPmFile(java.io.BufferedWriter)
     * 
     * @param B2s Array of bond list
     * 
     * @return Formated string buffer after matching database to the bond list
    */
    private StringBuffer Extract2Bonds(int[][] B2s) throws IOException
    {
            log.append(">>> Bond database matching started...\n");
            StringBuffer sb =new StringBuffer();
            
            bondlist = B2s;
            bondlength = new double[B2s.length];
            bondconst = new double[B2s.length];
            Bond2count = 0;
            for(int j = 0; j< B2s.length;j++)
            {
                try {
                boolean PairExists = false;
                    System.out.println("2Bond=>"+B2s[j][0]+"-"+B2s[j][1]);
                    String bondType1 = ElementsTable.getValueAt(B2s[j][0]-1, 6).toString();
                                        System.out.println(bondType1);
                    String bondType2 = ElementsTable.getValueAt(B2s[j][1]-1, 6).toString();
            
              BufferedReader in = new BufferedReader(new FileReader(UserDirectory+""+File.separator+"Database"+File.separator+"potparam.txt"));
              
            while(in.ready())
            {
                String line = in.readLine().trim();
                if(line.startsWith("&bond"))
                {
                    
//                    System.out.println(bondType1);
//                    System.out.println(bondType2);
                   String Junk = in.readLine();
                   String[] Type1 = Junk.split("=");
                   Type1[1] = Type1[1].replace("\\par", "");
                   Junk = in.readLine();
                   String[] Type2 = Junk.split("=");
                   Type2[1] = Type2[1].replace("\\par", "");
//                   System.out.println("Type[2]=\\s"+Type2[1]);
                   if((Type1[1].trim().equals(bondType1.trim()) && Type2[1].trim().equals(bondType2.trim())) || 
                           (Type1[1].trim().equals(bondType2.trim()) && Type2[1].trim().equals(bondType1.trim())))
                   {
                       Bond2count++;
                       Junk = in.readLine();
                       //Junk = in.readLine();
                       sb.append("Bond       =  "+B2s[j][0]+" "+B2s[j][1]+"\n");
                       System.out.println("Bond\t=\t"+B2s[j][0]+"\\s"+B2s[j][1]);
                       log.append(">>> match found: "+(B2s[j][0]-1)+" "+(B2s[j][1]-1)+"\n");
                       if(Junk.contains("bond_length"))
                       {
                           String[] BondLength = Junk.split("=");
                           BondLength[1] = BondLength[1].replace("\\par", "");
                           sb.append("R0      \t= "+BondLength[1]+"\n");
                           bondlength[j] = Double.parseDouble(BondLength[1].trim()); ///////// save bondlength for geomatery
                           System.out.println("R0      \t=\t"+BondLength[1]);
                       }
                       Junk = in.readLine();
                       if(Junk.contains("bonding_constant"))
                       {
                           String[] ForConst = Junk.split("=");
                           ForConst[1] = ForConst[1].replace("\\par", "");
                           sb.append("ForConst\t=  "+ForConst[1]+"\n\n");
                           bondconst[j] = Double.parseDouble(ForConst[1].trim()); ///////// save Forconstant for geomatery 
                           System.out.println("ForConst\t=\t"+ForConst[1]);
                       }
                       else
                       {
                           sb.append("ForConst   =  0\n\n");
                           System.out.println("ForConst\t=\t0");
                       }
                       PairExists = true;
                       break;
                   }
                       
                  
                }
            }
            in.close();
            if(!PairExists)
            {
                int a = B2s[j][0]-1, b = B2s[j][1]-1;
                JOptionPane.showMessageDialog(this,"Bond match not found: "+a+" and "+b+" for forcefields "+bondType1+" and "+bondType2+" respectively",
                        "Match Not Found",JOptionPane.ERROR_MESSAGE);
               log.append("!!! mactch not found: "+a+" and "+b+" for forcefields "+bondType1+" and "+bondType2+" respectively...\n");
               log.append(">>> process not complete...Exit\n"); 
               sb = null;
                break;
            }
            
                log.append(">>> Bonds matching process complete...\n");
         }
          catch (FileNotFoundException ex) {
            ex.printStackTrace();
            
            log.append("!!! Caught Exception while matching Bonds to database\n");
        }
            
    }
            
            return sb;
    }
    
    /**
     * 
     * Used to get database values from elements table for provided angles
     * Used to match the list of angles with combination inside the database file
     * Prepare a string buffer to write angles parameters in the .par file.
     * 
     * @see molecularModeling#CalculateGeometryActionPerformed(java.awt.event.ActionEvent)
     * @see molecularModeling#WriteDBPropertiesInPmFile(java.io.BufferedWriter)
     * 
     * @param B2s Array of Angle list
     * @param OP Operation to be carried out (calculate geometry or writing into file)
     * 
     * @return Formated string buffer after matching database to the Angle list
    */
    private StringBuffer Extract3Bonds(int[][] B2s, String  OP) throws IOException
    {
            log.append(">>> Angle database matching started...\n");
            StringBuffer sb = new StringBuffer();
            anglelist = B2s;
            angle = new double[B2s.length];
            angleconst = new double[B2s.length];
            Bond3count = 0;
            for(int j = 0; j< B2s.length ;j++)
            {
             try {
                boolean PairExists = false;
                    
                    System.out.println("3Bond=>"+B2s[j][0]+"-"+B2s[j][1]+"-"+B2s[j][2]);
                    
                    String bondType1 = ElementsTable.getValueAt(B2s[j][0]-1, 6).toString();
                    String bondType2 = ElementsTable.getValueAt(B2s[j][1]-1, 6).toString();
                    String bondType3 = ElementsTable.getValueAt(B2s[j][2]-1, 6).toString();
            
             BufferedReader in = new BufferedReader(new FileReader(UserDirectory+""+File.separator+"Database"+File.separator+"potparam.txt"));
            while(in.ready())
            {
                String line = in.readLine().trim();
                if(line.startsWith("&bending"))
                {
                    
//                    System.out.println(bondType1);
//                    System.out.println(bondType2);
                   String Junk = in.readLine();
                   String[] Type1 = Junk.split("=");
                   Type1[1] = Type1[1].replace("\\par", "");
//                   System.out.println("Type[1]=\\s"+Type1[1]);
                   Junk = in.readLine();
                   String[] Type2 = Junk.split("=");
                   Type2[1] = Type2[1].replace("\\par", "");
//                   System.out.println("Type[2]=\\s"+Type2[1]);
                   Junk = in.readLine();
                   String[] Type3 = Junk.split("=");
                   Type3[1] = Type3[1].replace("\\par", "");
                   
                   if((Type1[1].trim().equals(bondType1.trim()) && Type2[1].trim().equals(bondType2.trim()) && Type3[1].trim().equals(bondType3.trim()))
                           || (Type1[1].trim().equals(bondType3.trim()) && Type2[1].trim().equals(bondType2.trim()) && Type3[1].trim().equals(bondType1.trim())))
                   {
                       Bond3count++;
                       Junk = in.readLine();
                       sb.append("Angle      =  "+B2s[j][0]+" "+B2s[j][1]+" "+B2s[j][2]+"\n");
                       System.out.println("Angle\t=\t"+B2s[j][0]+"\\s"+B2s[j][1]);
                       log.append("\tmatch found: "+(B2s[j][0]-1)+" "+(B2s[j][1]-1)+" "+(B2s[j][2]-1)+"\n");
                       if(Junk.contains("bending_angle"))
                       {
                           String[] BendAngle = Junk.split("=");
                           BendAngle[1] = BendAngle[1].replace("\\par", "");
                           sb.append("Angle0     = "+BendAngle[1]+"\n");
                           angle[j] = Double.parseDouble(BendAngle[1].trim());// save angle0 for geomatery
                           System.out.println("Angle0  \t=\t"+BendAngle[1]);
                       }
                       Junk = in.readLine();
                       if(Junk.contains("bending_constant"))
                       {
                           String[] ForConst = Junk.split("=");
                           ForConst[1] = ForConst[1].replace("\\par", "");
                           sb.append("ForConst   = "+ForConst[1]+"\n\n");
                           angleconst[j] = Double.parseDouble(ForConst[1].trim());// save angle constant for geomatery
                           System.out.println("ForConst\t=\t"+ForConst[1]);
                       }
                       PairExists = true;
                       break;
                   }
                  
                }
            }
            in.close();
            if(!PairExists)
            {
                int a = B2s[j][0]-1, b = B2s[j][1]-1,c = B2s[j][2]-1;
                JOptionPane.showMessageDialog(this,"Angle match not found: "+a+", "+b+" and "+c+" for forcefields "
                        +bondType1+", "+bondType2+" and "+bondType3+" respectively","Match Not Found",JOptionPane.ERROR_MESSAGE);
                log.append("!!! match not found: "+a+", "+b+" and "+c+" for forcefields "+bondType1+", "+bondType2+" and "+bondType3+" respectively...\n");
               if(OP.equals("CalGeo"))
               {
                   log.append(">>> process not complete...Exit\n"); 
                   //sb = null;
                   break;
               }
            }
            
            log.append(">>> Angles matching process complete...\n");
            
                }
             
             catch (FileNotFoundException ex) {
                    ex.printStackTrace();
                    log.append("!!! Caught Exception while matching Angles to database\n");
                }
             
         }
            
          if(sb.toString().isEmpty())
                                   return null;
          else
             return sb;
    }
    
    /**
     * 
     * Used to get database values from elements table for provided dihedrals
     * Used to match the list of dihedrals with combination inside the database file
     * Prepare a string buffer to write dihedrals parameters in the .par file.
     * 
     * @see molecularModeling#CalculateGeometryActionPerformed(java.awt.event.ActionEvent)
     * @see molecularModeling#WriteDBPropertiesInPmFile(java.io.BufferedWriter)
     * 
     * @param B2s Array of Dihedral list
     * @param OP Operation to be carried out (calculate geometry or writing file)
     * 
     * @return Formated string buffer after matching database to the dihedral list
    */
    private StringBuffer Extract4Bonds(int[][] B2s, String OP) throws IOException
    {
        log.append(">>> Dihedrals database matching started...\n");
       StringBuffer sb = new StringBuffer();
       StringBuffer dihedrals = new StringBuffer();
            
            dihedlist = B2s;
            dihed = new double[B2s.length][];
            dihedconst = new double[B2s.length][];
            Bond4count = 0;
            
            for(int j = 0; j< B2s.length;j++)
            {
                try {
                boolean PairExists = false;
                    
                    System.out.println("4Bond=>"+B2s[j][0]+"-"+B2s[j][1]+"-"+B2s[j][2]+"-"+B2s[j][3]);
                    
                    String bondType1 = ElementsTable.getValueAt(B2s[j][0]-1, 6).toString();
                    String bondType2 = ElementsTable.getValueAt(B2s[j][1]-1, 6).toString();
                    String bondType3 = ElementsTable.getValueAt(B2s[j][2]-1, 6).toString();
                    String bondType4 = ElementsTable.getValueAt(B2s[j][3]-1, 6).toString();
            
              BufferedReader in = new BufferedReader(new FileReader(UserDirectory+""+File.separator+"Database"+File.separator+"potparam.txt"));
            while(in.ready())
            {
                String line = in.readLine().trim();
                if(line.startsWith("&torsion"))
                {
                    
                   String Junk = in.readLine();
                   String[] Type1 = Junk.split("=");
                   Type1[1] = Type1[1].replace("\\par", "");
//                   System.out.println("Type[1]=\\s"+Type1[1]);
                   Junk = in.readLine();
                   String[] Type2 = Junk.split("=");
                   Type2[1] = Type2[1].replace("\\par", "");
//                   System.out.println("Type[2]=\\s"+Type2[1]);
                   Junk = in.readLine();
                   String[] Type3 = Junk.split("=");
                   Type3[1] = Type3[1].replace("\\par", "");
                   
                   Junk = in.readLine();
                   String[] Type4 = Junk.split("=");
                   Type4[1] = Type4[1].replace("\\par", "");
                   
                   if((Type1[1].trim().equals(bondType1.trim()) && Type2[1].trim().equals(bondType2.trim())
                           && Type3[1].trim().equals(bondType3.trim()) && Type4[1].trim().equals(bondType4.trim())) ||
                          (Type1[1].trim().equals(bondType4.trim()) && Type2[1].trim().equals(bondType3.trim())
                           && Type3[1].trim().equals(bondType2.trim()) && Type4[1].trim().equals(bondType1.trim())))
                   {
                       Bond4count++;
                       dihedrals.append("Dihedral   =  "+B2s[j][0]+" "+B2s[j][1]+" "+B2s[j][2]+" "+B2s[j][3]+"\n");
                      // System.out.println("Dihedral\t=\t"+B2s[j][0]+"\\s"+B2s[j][1]);
                       log.append("\t match found: "+(B2s[j][0]-1)+" "+(B2s[j][1]-1)+" "+(B2s[j][2]-1)+" "+(B2s[j][3]-1)+"\n");
                       Junk = in.readLine();
                       //int nmax = 0;
                       sb.setLength(0);
                       ArrayList<Double> t_dihed = new ArrayList<Double>();
                       ArrayList<Double> t_dihed_const = new ArrayList<Double>();
                       while(!Junk.contains("/"))
                       {
                           
                       if(Junk.contains("torsion_c"))
                       {
                           String[] PotBarrier = Junk.split("=");
                           PotBarrier[1] = PotBarrier[1].replace("\\par", "");
                           
                           t_dihed_const.add(Double.parseDouble(PotBarrier[1].trim()));
                           System.out.println("PotBarrier\t=\t"+Double.parseDouble(PotBarrier[1].trim()));
                       }
                       
                       Junk = in.readLine();
                       
                       if(Junk.contains("gamma_c"))
                       {
                           String[] gamma = Junk.split("=");
                           gamma[1] = gamma[1].replace("\\par", "");
                          
                           t_dihed.add(Double.parseDouble(gamma[1].trim()));
                           System.out.println("gamma\t=\t"+gamma[1]);
                           Junk = in.readLine();
                       }
                       else
                       {
                          
                           t_dihed.add(0.0);
                           System.out.println("gamma\t=\t0");
                           
                       }
//                           
                       }
                           
                         for(int i = t_dihed_const.size() - 1; i >=0 ; i--)
                           {
                               if(t_dihed_const.get(i) == 0.0)
                               {
                                   t_dihed_const.remove(i);
                                   t_dihed.remove(i);
                               }
                               else
                                   break;
                           }
                           dihed[j] = new double[t_dihed.size()];
                           dihedconst[j] = new double[t_dihed_const.size()];
                           String num1; 
                           String num2;
                           for(int i=0;i<t_dihed_const.size();i++)
                           {
                               
                               dihed[j][i] = t_dihed.get(i);
                               dihedconst[j][i] = t_dihed_const.get(i);
                               num1 = String.format("%d", i);
                               num2 = String.format("%02d", i);
                               sb.append("ForConst"+num1+"  = "+t_dihed_const.get(i).toString()+"\n");
                               sb.append("gamma"+num2+"    = "+t_dihed.get(i).toString()+"\n");
                               
                           }
                           
                           sb.append("ScaleLJ14  =  0.0625\n");
                           sb.append("ScaleEl14  =  0.5\n\n");
                           int nmax = t_dihed.size() - 1;
                           dihedrals.append("nmax       =  "+nmax+"\n");
                           dihedrals.append(sb);
                       PairExists = true;
                       break;
                   }
                  
                }
            }
            in.close();
            if(!PairExists)
            {
                int a = B2s[j][0]-1, b = B2s[j][1]-1,c = B2s[j][2]-1, d = B2s[j][3]-1;
                JOptionPane.showMessageDialog(this,"Dihedral match not found: "+a+", "+b+", "+c+" and "+d+" for forcefields "
                        +bondType1+", "+bondType2+", "+bondType3+" and "+bondType4+" respectively","Match Not Found",JOptionPane.ERROR_MESSAGE);
                log.append("!!! match not found: "+a+", "+b+", "+c+" and "+d+" for forcefields "
                        +bondType1+", "+bondType2+", "+bondType3+" and "+bondType4+" respectively\n");
                if(OP.equals("CalGeo"))
                {
                    log.append(">>> process not complete...Exit\n");
                    //dihedrals = null;
                    break;
                }
            }
            log.append(">>> Dihedrals matching process complete...\n");
           
            } catch (FileNotFoundException ex) {
            ex.printStackTrace();
            log.append("!!! Caught Exception while matching Dihedrals to database\n");
        }
                
         }
            if(dihedrals.toString().isEmpty())
                                   return null;
          else
            return dihedrals;
    }
    
    /**
     * 
     * calls Extract2bonds, Extract3bonds and Extract4bonds
     * Write the prepared string buffers by them into the .par file
     * And rest of the required parameters for the par file.
     * 
     * @see molecularModeling#Extract2Bonds(int[][])
     * @see molecularModeling#Extract3Bonds(int[][], java.lang.String)
     * @see molecularModeling#Extract4Bonds(int[][], java.lang.String)
     * 
     * @param bw Reference of the Buffer Writer
     * 
     * @return True - False
     * 
    */
    private boolean WriteDBPropertiesInPmFile(BufferedWriter bw, ArrayList<ArrayList<Integer>> cons, ArrayList<ArrayList<Integer>> bond,
            ArrayList<ArrayList<Integer>> angle, ArrayList<ArrayList<Integer>> dihed) throws IOException
    {
        
        boolean bondsOK = false;
        int count = 1;
        if(ElementsTable.getRowCount() > 1)
        {
        int[][] ex = PrepareArrayForMatrix();
        
        BADExtraction BADEx = new BADExtraction();
        
                StringBuffer Bufferbonds = Extract2Bonds(BADEx.getBondsList(ex));
                if(Bufferbonds == null)
                            bondsOK = false;
                else{
               StringBuffer BufferAngle = Extract3Bonds(BADEx.getAnglesList(ex),"");
               if(BufferAngle != null)
                             count++;
               
               StringBuffer BufferDiheds = Extract4Bonds(BADEx.getDihedralList(ex),"");
               
               if(BufferDiheds != null)
                              count++;
               
                bw.append("NRotAxes   =  auto\n" +
                                    "\n" +
                                    "NIdfTypes  =  "+count+"\n\n"
                                +  "IdfType    =  Bond\n"
                                + "NIdfs      =  "+(Bond2count + bond.size())+"\n\n");
               bw.append(Bufferbonds);
               if(bond != null) {
                   for(int i = 0; i < bond.size(); i++) {
                       String _ids = "";
                       for(int j = 0; j < bond.get(i).size(); j++) {
                           _ids += " " + bond.get(i).get(j);
                       }
                       bw.append("Bond       = "+ _ids +"\n");
                       bw.append("R0         =  1.43\n");
                       bw.append("ForConst   =  0\n\n");
                   }
                   //bw.append("\n");
               }
               
               if(BufferAngle != null)
               {
                bw.append("IdfType    =  Angle\n" +
                        "NIdfs      =  "+(Bond3count + angle.size())+"\n\n");
               bw.append(BufferAngle);
               }
               if(angle != null) {
                   for(int i = 0; i < angle.size(); i++) {
                       String _ids = "";
                       for(int j = 0; j < angle.get(i).size(); j++) {
                           _ids += " " + angle.get(i).get(j);
                       }
                       bw.append("Angle      = "+ _ids +"\n");
                       bw.append("Angle0     =  108.5\n");
                       bw.append("ForConst   =  27700.0\n\n");
                   }
                   //bw.append("\n");
               }
               if(BufferDiheds != null || dihed != null)
               {
               bw.append("IdfType    =  Dihedral\n" +
                        "NIdfs      =  "+(Bond4count + dihed.size())+"\n\n");
               bw.append(BufferDiheds);
               for(int i = 0; i < dihed.size(); i++) {
                       String _ids = "";
                       for(int j = 0; j < dihed.get(i).size(); j++) {
                           _ids += " " + dihed.get(i).get(j);
                       }
                       bw.append("Dihedral   ="+ _ids +"\n");
                       bw.append("nmax       =   3\n" +
                                "ForConst0  = -29.17\n" +
                                "gamma00    =   0\n" +
                                "ForConst1  = 209.82\n" +
                                "gamma01    =   0\n" +
                                "ForConst2  =  29.17\n" +
                                "gamma02    =   0\n" +
                                "ForConst3  = 187.93\n" +
                                "gamma03    =   0\n" +
                                "ScaleLJ14  =   0.0625\n" +
                                "ScaleEl14  =   0.5\n\n");
                   }
               }
               
               if(cons != null) {
                   bw.append("\nNConstrU   =  "+ cons.size() +"\n");
                   for(int i = 0; i < cons.size(); i++) {
                       bw.append("\nConstraint = "+ cons.get(i).size() +"\n");
                       String _ids = "";
                       for(int j = 0; j < cons.get(i).size(); j++) {
                           _ids += " " + cons.get(i).get(j);
                       }
                       bw.append("SiteIDs    ="+ _ids +"\n");
                       bw.append("NRotAxes   =  auto\n");
                   }
               }
               else {
                    bw.append("\nNConstrU   =  0");
               }
               
               bw.append("\nNRotAxes   =  auto");
               bondsOK = true;
               log.append(">>> Parameters ready and file is ready to save...\n");
                }
        }
        
        return bondsOK;
    }
   
    /**
     * 
     * add dropdown menu to the elements table
     * create arraylist of dropdowns for reusing or later updates.
     * 
    */
    private void addComboTotable() throws IOException
    {
        int row = ElementsTable.getRowCount() - 1;
        Jbox.add(new JComboBox(DBvalues));
        rowEditor.setEditorAt(row, new DefaultCellEditor(Jbox.get(Jbox.size() - 1)));
        
        String[] types = {"CH","CH2","CH3","C","O"};
        JComboBox jb = new JComboBox(types);
        jb.setEditable(true);
        TJbox.add(jb);
             TrowEditor.setEditorAt(row, new DefaultCellEditor(TJbox.get(TJbox.size() - 1)));
    }
    
    /**
     * 
     * Get Site_types for the dropdown menus inside elements table.
     * 
     * @see molecularModeling#addTablelistener()
     * @see molecularModeling#RefreshDBbtnActionPerformed(java.awt.event.ActionEvent)
     * @see molecularModeling#molecularModeling()
     * 
     * @param Filter Sring to be filter
     * 
     * @return Array of filtered strings
    */
    private String[] getDBsNames(String Filter) throws IOException
    {
        StringBuffer nameBuffer = new StringBuffer();
        try {
            BufferedReader in = new BufferedReader(new FileReader(UserDirectory+""+File.separator+"Database"+File.separator+"potparam.txt"));
            while(in.ready())
            {
                String line = in.readLine().trim();
                if(line.startsWith("site_type"))
                {
                   String[] DBn = line.split("=");
                   DBn[1] = DBn[1].replace("\\par", "");
       
                   if(!Filter.isEmpty())
                   {
                       String Nvalue = DBn[1].trim();
                       String substr = Nvalue.substring(1,Nvalue.indexOf("-")).toUpperCase();
                       if(substr.equals(Filter.toUpperCase()))
                                    nameBuffer.append(DBn[1] + " ");
                   }
                   else
                      nameBuffer.append(DBn[1] + " "); 
                   
               
            }
            }
            
            in.close();
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
            log.append("!!! Caught Exception while reading Database types\n");
        }
        
        StringTokenizer st = new StringTokenizer(nameBuffer.toString());
       String[] valuesNam = new String[st.countTokens()];
       for(int i = 0; i<valuesNam.length;i++)
       {
           valuesNam[i] = st.nextToken();
       }
       
       return valuesNam;
    }
    
    /**
     * 
     * add all the required listener to the elements table
     * Click to highlight the atom shape
     * change values inside x,y,x coordinates and immediatly effect on scenegraph
     * Select Database name value and get values for epsilon, sigma and mass
     * Search type related database names.
     * 
     * @see molecularModeling#CreateElementsTable() 
    */
    private void addTablelistener()
    {        
        ElementsTable.getModel().addTableModelListener(new TableModelListener(){

            @Override
            public void tableChanged(TableModelEvent e) {
                
                int index = e.getColumn();
    switch (e.getType()) {
        
    case TableModelEvent.INSERT:
        
      break;
    case TableModelEvent.UPDATE:
        // for coordinates....
        if(ElementsTable.getSelectedRow() > -1)
        {
        if(index == 1 || index == 2 || index == 3)
        {
        int id = Integer.parseInt(String.valueOf(ElementsTable.getValueAt(ElementsTable.getSelectedRow(), 0)));
        float x = Float.parseFloat(String.valueOf(ElementsTable.getValueAt(ElementsTable.getSelectedRow(), 1)));
        float y = Float.parseFloat(String.valueOf(ElementsTable.getValueAt(ElementsTable.getSelectedRow(), 2)));
        float z = Float.parseFloat(String.valueOf(ElementsTable.getValueAt(ElementsTable.getSelectedRow(), 3)));
        String links = String.valueOf(String.valueOf(ElementsTable.getValueAt(ElementsTable.getSelectedRow(), 4)));
        
        if(!(x > -1 && x < 1) || !(y > -1 && y < 1) || !(z > -1 && z < 1))
                    {
                                        float maxi = 0.0f;
                                        for(int i=0;i < ElementsTable.getRowCount(); i++)
                                            {
                                                x = Float.parseFloat(String.valueOf(ElementsTable.getValueAt(i, 1)));
                                                y = Float.parseFloat(String.valueOf(ElementsTable.getValueAt(i, 2)));
                                                z = Float.parseFloat(String.valueOf(ElementsTable.getValueAt(i, 3)));

                                                 float Gr = Math.max(maxi, Math.abs(x));
                                                 Gr = Math.max(Gr, Math.abs(y));
                                                 maxi = Math.max(Gr, Math.abs(z));

                                                
                                            }
                                 
                                maxi = maxi + 1.0f;
                                
                                //System.out.println("max value : "+maxi);
                                for(int i=0;i < ElementsTable.getRowCount(); i++)
                                {
                                    x = Float.parseFloat(String.valueOf(ElementsTable.getValueAt(i, 1)));
                                    y = Float.parseFloat(String.valueOf(ElementsTable.getValueAt(i, 2)));
                                    z = Float.parseFloat(String.valueOf(ElementsTable.getValueAt(i, 3)));
                                    links = String.valueOf(String.valueOf(ElementsTable.getValueAt(i, 4)));
                                    
                                    x = x/maxi;
                                    y = y/maxi;
                                    z = z/maxi;
                                 
                                //System.out.println("divided value : x ="+x+"  y ="+y+"  z ="+z);
                                CreateGraphics.changeVertex(i,x,y,z,links);
                                }
                            }
                            else
                                CreateGraphics.changeVertex(id,x,y,z,links);
        
        CreateGraphics.HighlightSphere(id);
        
        //log.append(">>> "+id+" atom value changed\n");
        
        }
        // for databse name value and epsilon, sigma and mass...
        if(index == 6)
        {
            
            String value = String.valueOf(ElementsTable.getValueAt(ElementsTable.getSelectedRow(), 6));
                    try {
                        setDBProperties(value);
                        int id = Integer.parseInt(String.valueOf(ElementsTable.getValueAt(ElementsTable.getSelectedRow(), 0)));
                        log.append(">>> "+id+" atom selected database ["+value+"]\n");
                    } catch (IOException ex) {
                        Logger.getLogger(molecularModeling.class.getName()).log(Level.SEVERE, null, ex);
                    }
        }
        // for type search database names....
        if(index == 5)
        {
            if(ElementsTable.getSelectedRow() > -1)
            {
            String value = String.valueOf(ElementsTable.getValueAt(ElementsTable.getSelectedRow(), 5));
            int id = Integer.parseInt(String.valueOf(ElementsTable.getValueAt(ElementsTable.getSelectedRow(), 0)));
                    try {
                        //System.out.println("searched>>>>>>"+value);
                        String[] DBsNames = getDBsNames(value);
                        //System.out.println("values>>>>>>"+DBsNames.length);
                        Jbox.set(id, new JComboBox(DBsNames));
                       // System.out.println("searched>>>>>>"+id);
                        rowEditor.setEditorAt(id, new DefaultCellEditor(Jbox.get(id)));
                    } catch (IOException ex) {
                        Logger.getLogger(molecularModeling.class.getName()).log(Level.SEVERE, null, ex);
                    }
            }
        }
        }
      break;
    case TableModelEvent.DELETE:
      break;
    }
               // throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
            }
            
        });
        
        // to highlight the selected atom....
        ElementsTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            @Override
            public void valueChanged(ListSelectionEvent event) {
               
                if (ElementsTable.getSelectedRow() > -1) {
                    // print first column value from selected row
                    //System.out.println(ElementsTable.getValueAt(ElementsTable.getSelectedRow(), 1).toString());
                    int id = -1;
                    try{
                    id = Integer.parseInt(String.valueOf(ElementsTable.getValueAt(ElementsTable.getSelectedRow(), 0)));
                    CreateGraphics.HighlightSphere(id);
                    }
                    catch(Exception ex)
                    {
                         log.append("!!! An unexpected exception occured while selecting atom "+id+"\n");
                    }
                    //log.append(">>> "+id+" atom is selected\n");
                }
            }
});
    }
    
    /**
     *  get epsilon , sigma and mass when DB selected.
     * 
     * @param Value Name of the database value
     * 
     * @throws IOException 
     */
    private void setDBProperties(String Value) throws IOException
    {
        try {
            BufferedReader in = new BufferedReader(new FileReader(UserDirectory+""+File.separator+"Database"+File.separator+"potparam.txt"));
            while(in.ready())
            {
                String line = in.readLine().trim();
                if(line.startsWith("site_type") && line.contains(Value))
                {
                    int i=7;
                    while(i<10)
                    {
                   String nextline = in.readLine();
                   String[] DBn = nextline.split("=");
                   DBn[1] = DBn[1].replace("\\par", "");
                   this.ElementsTable.setValueAt(DBn[1], ElementsTable.getSelectedRow(), i);
                   i++;
                    }
                }
               
            }
            
            in.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(molecularModeling.class.getName()).log(Level.SEVERE, null, ex);
            log.append("!!! Caught Exception while setting Database types \n");
        }
        
    }

    /**
     * 
     * used to call from CreateGraphics class to add
     * row in elements table when user click to add new atom
     * also in readfile function to add rows and their specified information.
     * 
     * @see CreateGraphics#CreateSphere(com.sun.j3d.utils.geometry.Primitive)
     * 
     * @param mSphEre Name of current atom
     * @param LspHEre Name of last atom
     * @param mSphcord Old coordinates
     * @param LSphcord New coordinates
     * 
     */
    
    public void addTableRows(String mSphEre,String LspHEre,Vector3f mSphcord,Vector3f LSphcord)
    {
        if(mSphEre !=null && LspHEre !=null)
        {
             
           // JOptionPane.showMessageDialog(null, mSphEre +" / "+LspHEre);
        String[] NSo = mSphEre.split("-");
        String[] NSn = LspHEre.split("-");
        boolean main = false;
        DefaultTableModel tm = (DefaultTableModel) ElementsTable.getModel();
        int rowCount = tm.getRowCount();
       for(int i=0;i<rowCount;i++)
       {
         String row = tm.getValueAt(i, 0).toString();
         if(NSo[1].equals(row))
         {
            // JOptionPane.showMessageDialog(null, row);
             int r = Integer.parseInt(row);
            //JOptionPane.showMessageDialog(null,r); 
             String Links = tm.getValueAt(i,4).toString();
             if(Links.equals(""))
                 Links = NSn[1];
             else
                 Links = Links + "," + NSn[1];
           // JOptionPane.showMessageDialog(null,Links);
            //tm.removeRow(r);
            //tm.insertRow(r, new Object[]{row,mSphcord.x,mSphcord.y,mSphcord.z,Links});
             tm.setValueAt(Links, r, 4);
             main = true;
         }
       }
       if(!main)
       {
           float x = (float) mSphcord.x * 10.0f;
           float y = (float) mSphcord.y * 10.0f;
           float z = (float) mSphcord.z * 10.0f;
           
           //tm.addRow(new Object[]{NSo[1],mSphcord.x,mSphcord.y,mSphcord.z,NSn[1]});
           tm.addRow(new Object[]{NSo[1],x,y,z,NSn[1]});
       }
       main = false;
           float x = (float) LSphcord.x * 10.0f;
           float y = (float) LSphcord.y * 10.0f;
           float z = (float) LSphcord.z * 10.0f;
          
          // tm.addRow(new Object[]{NSn[1],LSphcord.x,LSphcord.y,LSphcord.z,NSo[1]});
          tm.addRow(new Object[]{NSn[1],x,y,z,NSo[1]});
          log.append(">>> atom successfully added...\n");
       undobtn.setEnabled(true);
        }
        
        try {
                 
            addComboTotable();
            
        } catch (IOException ex) {
            Logger.getLogger(molecularModeling.class.getName()).log(Level.SEVERE, null, ex);
            log.append("!!! Caught Exception while adding combo box to table \n");
        }
        
    }

    /**
     * 
     * used to update the "links" coloumn fo the elements table
     * when there is new atom added or deleted by user
     * 
     * @param i_id1 Id of one atom
     * @param i_id2 Id of second atom
     */
    public void EditTableRows(int i_id1,int i_id2)
    {
        DefaultTableModel tm = (DefaultTableModel) ElementsTable.getModel();
        int rowCount = tm.getRowCount();
       for(int i=0;i<rowCount;i++)
       {
           if(i_id1 == i)
           {
               String link = tm.getValueAt(i, 4).toString();
               link += ","+i_id2;
               tm.setValueAt(link, i, 4);
           }
           else
               if(i_id2 == i)
               {
               String link = tm.getValueAt(i, 4).toString();
               link += ","+i_id1;
               tm.setValueAt(link, i, 4);
               }
       }
    }
    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(molecularModeling.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(molecularModeling.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(molecularModeling.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(molecularModeling.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new molecularModeling().setVisible(true);
            }
        });
    }
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton CalculateGeometry;
    private javax.swing.JInternalFrame CanvasInternalFrame;
    private javax.swing.JPopupMenu ElementTablePopupmenu;
    private javax.swing.JButton OpenBtn;
    private javax.swing.JButton RefreshDBbtn;
    private javax.swing.JButton Restartbtn;
    private javax.swing.JButton Savebtn;
    private javax.swing.JButton addDBentry;
    private javax.swing.JToolBar canvasmenu;
    private javax.swing.JButton centerofmassbtn;
    private javax.swing.JLabel error;
    private javax.swing.JLabel errorlbl;
    private javax.swing.JButton helpbtn;
    private javax.swing.JMenu jMenu3;
    private javax.swing.JPopupMenu jPopupMenu1;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JScrollPane jScrollPane2;
    public javax.swing.JPanel leftPanelM;
    private javax.swing.JButton ljsitebtn;
    public static javax.swing.JTextArea log;
    private javax.swing.JPanel logpanel;
    private javax.swing.JToolBar mainmenu;
    private javax.swing.JButton resetbtn;
    private javax.swing.JButton undobtn;
    // End of variables declaration//GEN-END:variables
}
