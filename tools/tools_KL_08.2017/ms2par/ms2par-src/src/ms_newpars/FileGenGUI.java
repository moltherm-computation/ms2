/*
 * FileGenGUI.java
 *
 * Created on 20. Oktober 2006, 11:06
 */

package ms_newpars;

import java.awt.Component;
import java.awt.Desktop;
import java.awt.Event;
import java.io.*;
import java.util.*;
import javax.swing.*;
import javax.swing.filechooser.FileFilter;
import java.awt.GraphicsDevice;
import java.awt.GraphicsEnvironment;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.util.Iterator;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

/**
 * The GUI for MS parameter files.
 * <P>
 * 0.1a GUI Design since 10/2006
 * 0.1.1a ComboBoxes
 * 0.1.2a EnsemlePanel -> external File
 * 0.1.3a GeFrames added
 * 0.1.4a class for PmTable
 * 0.1.5a FileIO test added
 * 0.1.6a Table Pm, Button GE, etaxi
 * 0.1.7a Datastructure classes
 * 0.1.8a complete code review
 * 0.1.9a new read and write methods
 * 0.2a first alpha release 2/2007
 * 0.2.1a added MCorSteps and VapDensity
 * 0.2.2a corrected unit calc reduced units
 * 0.2.3a clearButton and values corrected
 * 0.2.4a selection values corrected MD MC MCOR
 * 0.2.5a added copy for *.pm files
 * 0.2.6a added EtaXi
 * 0.2.7a added Dialog for existing pm files
 * 0.2.8a EtaXi Tableeditorqu
 * 0.2.9a PmList Editor MolFrac
 * 0.3a 2nd alpha first release 4/2007
 * 0.3.1a added chemPotMethod GradIns
 * 0.3.2a changed GeData validate Input
 * 0.3.3a added Button and Frame for pmEditor
 * 0.3.4a added res2ge
 * 0.3.5a corrected partmolvol in res2ge
 * 0.3.6a corrected molFrac added 0.0 as possible value
 * 0.3.7a added read chemPotMethod from parFile
 * 0.3.8a removed vapdensity for ge / partmolvol only reduced
 * 0.3.9a read and print gi parameter
 * 0.4a editor for Potential Models
 * 0.4.1a corrected pmOverwrite name Problem
 * 0.4.2a fixed pmbutton problem with savedialog
 * 0.5a import resultfiles -- 3rd release 5/2007
 * 0.5.1a autocorrection for molfrac sum=1
 * 0.5.2a corrected NTPSteps Output only for NPT
 * 0.5.3a removed NTest with code
 * 0.5.4a moved type of ensemble and simulation in GUI and Output
 * 0.5.5a correced some combobox typeof enabled
 * 0.5.6a corrected GE output for PartMolVol (enabled <2)
 * 0.6a added SVC simulation -- 4th release 6/2007
 * 0.6.1a corrected output eta/xi lower case
 * 0.6.2a fixed bug multiple _ge adding in save dialog
 * 0.6.3a fixed system of units combo box -> gdata
 * 0.6.4a added initial value for cuttoff
 * 0.6.5a fixed calcmaxcutoff for reduced systemofunits
 * 0.6.6a removed chempotmethod for GE
 * 0.6.7a fixed chempotmethod button for svc
 * 0.6.8a fixed reading mueVTSteps into NPT
 * 0.6.9a fixed cpmrow in chempotmethod 21.8.07
 * 0.7b   5th release 1st beta 8/2007
 * 0.7.1b fixed EtaXi empty lines in par
 *
 * 0.8.1a find solution for save dialog
 *
 * </P>
 * @author  Christian Berreth
 * @version 0.7b
 */
public class FileGenGUI extends javax.swing.JFrame {
    private File loadFile = new File(System.getProperty("user.dir"));  //user.home
    private File loadDir = new File(System.getProperty("user.dir"));   //user.home

    /**
     *
     */
    public File saveFile = null;
    private int importRes = 0;
    private String UserDirectory = null;
    private String lineSep = System.getProperty("line.separator");
    private String fileSep = System.getProperty("file.separator");

    /**
     *
     */
    public HashMap<String,String> helpMap = new HashMap<String,String>();
    /** Creates new form FileGenGUI */
    public FileGenGUI() {
        
        try {
            // Set cross-platform Java L&F
        UIManager.setLookAndFeel(
            UIManager.getCrossPlatformLookAndFeelClassName());
    } 
    catch (UnsupportedLookAndFeelException e) {
       // handle exception
    }
    catch (ClassNotFoundException e) {
       // handle exception
    }
    catch (InstantiationException e) {
       // handle exception
    }
    catch (IllegalAccessException e) {
       // handle exception
    }
        initComponents();
        
        initPanel();
        readHelpxlsxfile();
        addEnsemble();
        showConsoleMessage();
        setwindowsize();
        settooltipforlabels();
        addShortCutsForButtonActions();
        
    }
    
    /**
    * Get all the label and buttons from the UI
    * Assign tool tip from the populated hashMap (helpmap)
    * Call addlinebreaks function to make tool tip in paragraph form
    * use html tags inside tool tip to formate paragrah form.
    */
    private void settooltipforlabels()
    {
        for(Component cm: generalPanel.getComponents())
        {
            if(cm instanceof JLabel)
            {
                JLabel lbl = (JLabel) cm;
                if(helpMap.containsKey(lbl.getText().toLowerCase())) {
                String tooltip = addLinebreaks(helpMap.get(lbl.getText().toLowerCase()),30);
                lbl.setToolTipText("<html><body><p>"+tooltip+"</p></body></html>");
                }
            }
            else
            if(cm instanceof JButton)
            {
                JButton btn = (JButton) cm;
                if(helpMap.containsKey(btn.getText().toLowerCase())) {
                String tooltip = addLinebreaks(helpMap.get(btn.getText().toLowerCase()),30);
                btn.setToolTipText("<html><body><p>"+tooltip+"</p></body></html>");
                }
            }
            else
                if(cm instanceof JPanel)
                {
                    JPanel jp = (JPanel) cm;
                    for(Component cmp: jp.getComponents())
                    {
                        if(cmp instanceof JLabel)
                        {
                            JLabel lbl = (JLabel) cmp;
                            if(helpMap.containsKey(lbl.getText().toLowerCase())) {
                            String tooltip = addLinebreaks(helpMap.get(lbl.getText().toLowerCase()),30);
                            lbl.setToolTipText("<html><body><p>"+tooltip+"</p></body></html>");
                            }
                        }
                        
                        if(cmp instanceof JButton)
                        {
                            JButton btn = (JButton) cmp;
                            if(helpMap.containsKey(btn.getText().toLowerCase())) {
                            String tooltip = addLinebreaks(helpMap.get(btn.getText().toLowerCase()),30);
                            btn.setToolTipText("<html><body><p>"+tooltip+"</p></body></html>");
                            }
                        }
                    }
                }
        }
    }

    /**
     * To add line breaks after every specified characters length
     * Input: String which needed to be break in line
     * maxLineLength: maximum length specified to break the line.
     * 
     * @param input Help text
     * @param maxLineLength Length of the line
     * 
     * @return formatted help text
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
     * Read help file from the help folder
     * Find the "Name" & "Meaning" coloumn inside the help xlsx file (we need only these two)
     * populate the hashmap (helpMap) to save all the names and meanings
            Key: Text from the "Name" column
            Value: Text from the "Meaning"" column.
     */
    private void readHelpxlsxfile()
    {
        try 
        {
        // Get the workbook instance for XLS file
        String filePath = UserDirectory+""+File.separator+"help"+File.separator+"ms2par_help_format.xlsx";
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
                   helpMap.put(cell.getStringCellValue().toLowerCase(), mcell.getStringCellValue());
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
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        EnsembleAddDialog = new javax.swing.JDialog();
        mainjscroll = new javax.swing.JScrollPane();
        mainjpanel = new javax.swing.JPanel();
        topjscroll = new javax.swing.JScrollPane();
        topjpanel = new javax.swing.JPanel();
        loadButton = new javax.swing.JButton();
        saveButton = new javax.swing.JButton();
        clearButton = new javax.swing.JButton();
        jLabel19 = new javax.swing.JLabel();
        jLabel20 = new javax.swing.JLabel();
        ensemblesaddButton = new javax.swing.JButton();
        ensemblesremoveButton = new javax.swing.JButton();
        helpButton = new javax.swing.JButton();
        filenameLabel = new javax.swing.JLabel();
        leftjscroll = new javax.swing.JScrollPane();
        generalPanel = new javax.swing.JPanel();
        jLabel1 = new javax.swing.JLabel();
        systemofunitsComboBox = new javax.swing.JComboBox();
        jLabel2 = new javax.swing.JLabel();
        lengthTextField = new javax.swing.JTextField();
        jLabel3 = new javax.swing.JLabel();
        jLabel4 = new javax.swing.JLabel();
        energyTextField = new javax.swing.JTextField();
        jLabel6 = new javax.swing.JLabel();
        jLabel5 = new javax.swing.JLabel();
        massTextField = new javax.swing.JTextField();
        jLabel7 = new javax.swing.JLabel();
        simPanel = new javax.swing.JPanel();
        jLabel9 = new javax.swing.JLabel();
        typeofsimulationComboBox = new javax.swing.JComboBox();
        jLabel10 = new javax.swing.JLabel();
        integratorComboBox = new javax.swing.JComboBox();
        jLabel11 = new javax.swing.JLabel();
        timestepTextField = new javax.swing.JTextField();
        jLabel12 = new javax.swing.JLabel();
        acceptanceTextField = new javax.swing.JTextField();
        jLabel22 = new javax.swing.JLabel();
        mcorStepsTextField = new javax.swing.JTextField();
        nOrientTextField = new javax.swing.JTextField();
        jLabel23 = new javax.swing.JLabel();
        jLabel24 = new javax.swing.JLabel();
        rStepsTextField = new javax.swing.JTextField();
        jLabel26 = new javax.swing.JLabel();
        rMinRadiusTextField = new javax.swing.JTextField();
        optPressureComboBox = new javax.swing.JComboBox();
        optpressurelabel = new javax.swing.JLabel();
        commonequiComboBox = new javax.swing.JComboBox();
        commonequilabel = new javax.swing.JLabel();
        jLabel25 = new javax.swing.JLabel();
        rMaxRadiusTextField = new javax.swing.JTextField();
        timestepunitlabel = new javax.swing.JLabel();
        jLabel13 = new javax.swing.JLabel();
        nvtstepsTextField = new javax.swing.JTextField();
        jLabel14 = new javax.swing.JLabel();
        nptstepsTextField = new javax.swing.JTextField();
        jLabel15 = new javax.swing.JLabel();
        runstepsTextField = new javax.swing.JTextField();
        jLabel16 = new javax.swing.JLabel();
        resultfreqTextField = new javax.swing.JTextField();
        jLabel17 = new javax.swing.JLabel();
        errorfreqTextField = new javax.swing.JTextField();
        jLabel18 = new javax.swing.JLabel();
        visualfreqTextField = new javax.swing.JTextField();
        jLabel21 = new javax.swing.JLabel();
        cutofftypeComboBox = new javax.swing.JComboBox();
        jLabel8 = new javax.swing.JLabel();
        typeofensembleComboBox = new javax.swing.JComboBox();
        jLabel28 = new javax.swing.JLabel();
        longRangeCorrComboBox = new javax.swing.JComboBox();
        jLabel29 = new javax.swing.JLabel();
        kappaTextField = new javax.swing.JTextField();
        jLabel30 = new javax.swing.JLabel();
        nVecMaxTextField = new javax.swing.JTextField();
        jLabel31 = new javax.swing.JLabel();
        nsqMaxTextField = new javax.swing.JTextField();
        jLabel32 = new javax.swing.JLabel();
        nMaxTextField = new javax.swing.JTextField();
        jLabel33 = new javax.swing.JLabel();
        intDegFreedComboBox = new javax.swing.JComboBox();
        jLabel34 = new javax.swing.JLabel();
        intraLjElComboBox = new javax.swing.JComboBox();
        jLabel35 = new javax.swing.JLabel();
        ljElComboBox = new javax.swing.JComboBox();
        jLabel36 = new javax.swing.JLabel();
        qshaker = new javax.swing.JComboBox();
        jLabel37 = new javax.swing.JLabel();
        printidf = new javax.swing.JComboBox();
        jLabel38 = new javax.swing.JLabel();
        jLabel39 = new javax.swing.JLabel();
        tolerance = new javax.swing.JTextField();
        scale = new javax.swing.JTextField();
        mvtlabel = new javax.swing.JLabel();
        mvttxtbox = new javax.swing.JTextField();
        rdflabel = new javax.swing.JLabel();
        rdftxtbox = new javax.swing.JTextField();
        numshelllabel = new javax.swing.JLabel();
        numshelltxtbox = new javax.swing.JTextField();
        jLabel27 = new javax.swing.JLabel();
        versiontxt = new javax.swing.JTextField();
        jLabel40 = new javax.swing.JLabel();
        jLabel41 = new javax.swing.JLabel();
        ndbtxt = new javax.swing.JTextField();
        viewparbtn = new javax.swing.JButton();
        viewresbtn = new javax.swing.JButton();
        jLabel42 = new javax.swing.JLabel();
        walltimetxt = new javax.swing.JTextField();
        rightjscroll = new javax.swing.JScrollPane();
        ensemblesTabbedPane = new javax.swing.JTabbedPane();

        EnsembleAddDialog.setDefaultCloseOperation(javax.swing.WindowConstants.DO_NOTHING_ON_CLOSE);

        org.jdesktop.layout.GroupLayout EnsembleAddDialogLayout = new org.jdesktop.layout.GroupLayout(EnsembleAddDialog.getContentPane());
        EnsembleAddDialog.getContentPane().setLayout(EnsembleAddDialogLayout);
        EnsembleAddDialogLayout.setHorizontalGroup(
            EnsembleAddDialogLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(0, 400, Short.MAX_VALUE)
        );
        EnsembleAddDialogLayout.setVerticalGroup(
            EnsembleAddDialogLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(0, 300, Short.MAX_VALUE)
        );

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setTitle("ms2par - Editor for MS Parameter Files");

        mainjscroll.setBorder(javax.swing.BorderFactory.createBevelBorder(javax.swing.border.BevelBorder.RAISED));

        mainjpanel.setPreferredSize(new java.awt.Dimension(1178, 1118));

        topjscroll.setBorder(null);
        topjscroll.setHorizontalScrollBarPolicy(javax.swing.ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
        topjscroll.setVerticalScrollBarPolicy(javax.swing.ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER);

        topjpanel.setBorder(javax.swing.BorderFactory.createBevelBorder(javax.swing.border.BevelBorder.LOWERED));
        topjpanel.setPreferredSize(new java.awt.Dimension(1100, 50));

        loadButton.setText("Load");
        loadButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                loadButtonActionPerformed(evt);
            }
        });

        saveButton.setText("Save");
        saveButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                saveButtonActionPerformed(evt);
            }
        });

        clearButton.setText("Clear");
        clearButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                clearButtonActionPerformed(evt);
            }
        });

        jLabel19.setText("Number of Ensembles:");

        jLabel20.setText("n");

        ensemblesaddButton.setText("Add ensemble");
        ensemblesaddButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ensemblesaddButtonActionPerformed(evt);
            }
        });

        ensemblesremoveButton.setText("Remove ensemble");
        ensemblesremoveButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ensemblesremoveButtonActionPerformed(evt);
            }
        });

        helpButton.setText("Help");
        helpButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                helpButtonActionPerformed(evt);
            }
        });

        filenameLabel.setText("jLabel42");

        org.jdesktop.layout.GroupLayout topjpanelLayout = new org.jdesktop.layout.GroupLayout(topjpanel);
        topjpanel.setLayout(topjpanelLayout);
        topjpanelLayout.setHorizontalGroup(
            topjpanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(topjpanelLayout.createSequentialGroup()
                .addContainerGap()
                .add(loadButton, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 84, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(saveButton, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 77, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .add(18, 18, 18)
                .add(filenameLabel)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED, 569, Short.MAX_VALUE)
                .add(clearButton, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 90, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .add(18, 18, 18)
                .add(jLabel19)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.UNRELATED)
                .add(jLabel20, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 19, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .add(33, 33, 33)
                .add(ensemblesaddButton)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(ensemblesremoveButton)
                .add(115, 115, 115)
                .add(helpButton)
                .addContainerGap())
        );
        topjpanelLayout.setVerticalGroup(
            topjpanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(topjpanelLayout.createSequentialGroup()
                .add(13, 13, 13)
                .add(topjpanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(loadButton)
                    .add(saveButton, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(clearButton, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(jLabel19)
                    .add(jLabel20)
                    .add(filenameLabel))
                .addContainerGap(org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
            .add(topjpanelLayout.createSequentialGroup()
                .add(12, 12, 12)
                .add(topjpanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(ensemblesaddButton)
                    .add(ensemblesremoveButton, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(helpButton, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap(org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        topjscroll.setViewportView(topjpanel);

        leftjscroll.setBorder(null);
        leftjscroll.setHorizontalScrollBarPolicy(javax.swing.ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
        leftjscroll.setVerticalScrollBarPolicy(javax.swing.ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER);

        generalPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("General settings"));
        generalPanel.setAutoscrolls(true);
        generalPanel.setMinimumSize(new java.awt.Dimension(100, 100));
        generalPanel.setPreferredSize(new java.awt.Dimension(450, 1030));

        jLabel1.setText("System of Units");

        systemofunitsComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "SI", "Reduced" }));
        systemofunitsComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                systemofunitsComboBoxActionPerformed(evt);
            }
        });

        jLabel2.setText("Sigma");

        lengthTextField.setText("0.0");
        lengthTextField.setMinimumSize(new java.awt.Dimension(69, 19));
        lengthTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                FileGenGUI.this.focusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                lengthTextFieldFocusLost(evt);
            }
        });
        lengthTextField.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                lengthTextFieldActionPerformed(evt);
            }
        });

        jLabel3.setText("[Å]");

        jLabel4.setText("Epsilon");

        energyTextField.setText("0.0");
        energyTextField.setMinimumSize(new java.awt.Dimension(69, 19));
        energyTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                FileGenGUI.this.focusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                energyTextFieldFocusLost(evt);
            }
        });

        jLabel6.setText("[K]");

        jLabel5.setText("Mass");

        massTextField.setText("0.0");
        massTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                FileGenGUI.this.focusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                massTextFieldFocusLost(evt);
            }
        });

        jLabel7.setText("[u]");

        jLabel9.setText("Type of Simulation");

        typeofsimulationComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "MD", "MC", "SVC" }));
        typeofsimulationComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                typeofsimulationComboBoxActionPerformed(evt);
            }
        });

        jLabel10.setText("Integrator");

        integratorComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Gear", "Leapfrog" }));
        integratorComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                integratorComboBoxActionPerformed(evt);
            }
        });

        jLabel11.setText("Time Step");

        timestepTextField.setText("2.0");
        timestepTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                FileGenGUI.this.focusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                timestepTextFieldFocusLost(evt);
            }
        });

        jLabel12.setText("Acceptance");
        jLabel12.setEnabled(false);

        acceptanceTextField.setText("0.0");
        acceptanceTextField.setEnabled(false);
        acceptanceTextField.setPreferredSize(new java.awt.Dimension(69, 19));
        acceptanceTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                FileGenGUI.this.focusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                acceptanceTextFieldFocusLost(evt);
            }
        });

        jLabel22.setText("MCORSteps");

        mcorStepsTextField.setText("0");
        mcorStepsTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                FileGenGUI.this.focusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                mcorStepsTextFieldFocusLost(evt);
            }
        });

        nOrientTextField.setText("1000");
        nOrientTextField.setPreferredSize(new java.awt.Dimension(69, 19));
        nOrientTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                FileGenGUI.this.focusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                nOrientTextFieldFocusLost(evt);
            }
        });

        jLabel23.setText("NOrient");
        jLabel23.setMaximumSize(new java.awt.Dimension(66, 13));
        jLabel23.setMinimumSize(new java.awt.Dimension(66, 13));
        jLabel23.setPreferredSize(new java.awt.Dimension(66, 13));

        jLabel24.setText("RSteps");
        jLabel24.setMaximumSize(new java.awt.Dimension(56, 13));
        jLabel24.setMinimumSize(new java.awt.Dimension(56, 13));

        rStepsTextField.setText("100");
        rStepsTextField.setPreferredSize(new java.awt.Dimension(69, 19));
        rStepsTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                FileGenGUI.this.focusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                rStepsTextFieldFocusLost(evt);
            }
        });

        jLabel26.setText("RMin");
        jLabel26.setMaximumSize(new java.awt.Dimension(56, 13));
        jLabel26.setMinimumSize(new java.awt.Dimension(56, 13));

        rMinRadiusTextField.setText("1.8");
        rMinRadiusTextField.setPreferredSize(new java.awt.Dimension(69, 19));
        rMinRadiusTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                FileGenGUI.this.focusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                rMinRadiusTextFieldFocusLost(evt);
            }
        });

        optPressureComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));
        optPressureComboBox.setEnabled(false);
        optPressureComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                optPressureComboBoxActionPerformed(evt);
            }
        });

        optpressurelabel.setText("Explicit Pressure");
        optpressurelabel.setEnabled(false);

        commonequiComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Yes", "No" }));
        commonequiComboBox.setEnabled(false);
        commonequiComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                commonequiComboBoxActionPerformed(evt);
            }
        });

        commonequilabel.setText("Common Equilibration");
        commonequilabel.setEnabled(false);

        jLabel25.setText("RMax");
        jLabel25.setMaximumSize(new java.awt.Dimension(56, 13));
        jLabel25.setMinimumSize(new java.awt.Dimension(56, 13));

        rMaxRadiusTextField.setText("6.0");
        rMaxRadiusTextField.setPreferredSize(new java.awt.Dimension(69, 19));
        rMaxRadiusTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                FileGenGUI.this.focusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                rMaxRadiusTextFieldFocusLost(evt);
            }
        });

        timestepunitlabel.setText("femtosec");

        org.jdesktop.layout.GroupLayout simPanelLayout = new org.jdesktop.layout.GroupLayout(simPanel);
        simPanel.setLayout(simPanelLayout);
        simPanelLayout.setHorizontalGroup(
            simPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(simPanelLayout.createSequentialGroup()
                .add(simPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING)
                    .add(jLabel9, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 228, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(simPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING, false)
                        .add(org.jdesktop.layout.GroupLayout.TRAILING, jLabel10, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .add(org.jdesktop.layout.GroupLayout.TRAILING, jLabel11, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .add(org.jdesktop.layout.GroupLayout.TRAILING, jLabel22, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .add(org.jdesktop.layout.GroupLayout.TRAILING, jLabel12, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .add(org.jdesktop.layout.GroupLayout.TRAILING, optpressurelabel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .add(org.jdesktop.layout.GroupLayout.TRAILING, commonequilabel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .add(org.jdesktop.layout.GroupLayout.TRAILING, jLabel23, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .add(org.jdesktop.layout.GroupLayout.TRAILING, jLabel24, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .add(org.jdesktop.layout.GroupLayout.TRAILING, jLabel26, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .add(org.jdesktop.layout.GroupLayout.TRAILING, jLabel25, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                .add(27, 27, 27)
                .add(simPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, nOrientTextField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, commonequiComboBox, 0, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, acceptanceTextField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, mcorStepsTextField)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, timestepTextField)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, integratorComboBox, 0, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, typeofsimulationComboBox, 0, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(optPressureComboBox, 0, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, rStepsTextField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(rMinRadiusTextField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(rMaxRadiusTextField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .add(27, 27, 27)
                .add(timestepunitlabel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 98, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .add(115, 115, 115))
        );
        simPanelLayout.setVerticalGroup(
            simPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(simPanelLayout.createSequentialGroup()
                .addContainerGap()
                .add(simPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel9)
                    .add(typeofsimulationComboBox, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(simPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel10)
                    .add(integratorComboBox, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(simPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(timestepTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(jLabel11)
                    .add(timestepunitlabel))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(simPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel22)
                    .add(mcorStepsTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(simPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel12)
                    .add(acceptanceTextField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .add(simPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(optPressureComboBox, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(optpressurelabel))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(simPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(commonequilabel)
                    .add(commonequiComboBox, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(simPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel23, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 26, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(nOrientTextField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(simPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel24, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(rStepsTextField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(simPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel26, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(rMinRadiusTextField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(simPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel25, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(rMaxRadiusTextField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap(org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        jLabel13.setText("NVT Steps");

        nvtstepsTextField.setText("jTextField5");
        nvtstepsTextField.setMinimumSize(new java.awt.Dimension(69, 19));
        nvtstepsTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                FileGenGUI.this.focusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                nvtstepsTextFieldFocusLost(evt);
            }
        });
        nvtstepsTextField.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                nvtstepsTextFieldActionPerformed(evt);
            }
        });

        jLabel14.setText("NPT Steps");
        jLabel14.setEnabled(false);

        nptstepsTextField.setText("jTextField6");
        nptstepsTextField.setEnabled(false);
        nptstepsTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                FileGenGUI.this.focusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                nptstepsTextFieldFocusLost(evt);
            }
        });

        jLabel15.setText("Run Steps");

        runstepsTextField.setText("jTextField7");
        runstepsTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                FileGenGUI.this.focusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                runstepsTextFieldFocusLost(evt);
            }
        });

        jLabel16.setText("Average Result Frequency");

        resultfreqTextField.setText("jTextField8");
        resultfreqTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                FileGenGUI.this.focusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                resultfreqTextFieldFocusLost(evt);
            }
        });

        jLabel17.setText("Error Frequency");

        errorfreqTextField.setText("jTextField9");
        errorfreqTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                FileGenGUI.this.focusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                errorfreqTextFieldFocusLost(evt);
            }
        });

        jLabel18.setText("Visual Capture Frequency");

        visualfreqTextField.setText("jTextField10");
        visualfreqTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                FileGenGUI.this.focusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                visualfreqTextFieldFocusLost(evt);
            }
        });

        jLabel21.setText("Cutoff-mode");

        cutofftypeComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "COM", "Site" }));
        cutofftypeComboBox.setEnabled(false);
        cutofftypeComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                cutofftypeComboBoxActionPerformed(evt);
            }
        });

        jLabel8.setText("Ensemble");

        typeofensembleComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "NVE", "NVT", "NPH", "NPT", "GE" }));
        typeofensembleComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                typeofensembleComboBoxActionPerformed(evt);
            }
        });

        jLabel28.setText("Long Range Correction");

        longRangeCorrComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "ReactionField", "Ewald", "Rodgers" }));
        longRangeCorrComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                longRangeCorrComboBoxActionPerformed(evt);
            }
        });

        jLabel29.setText("Kappa");

        kappaTextField.setText("jTextField1");
        kappaTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                FileGenGUI.this.focusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                kappaTextFieldFocusLost(evt);
            }
        });

        jLabel30.setText("NVecMax");

        nVecMaxTextField.setText("jTextField1");
        nVecMaxTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                FileGenGUI.this.focusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                nVecMaxTextFieldFocusLost(evt);
            }
        });

        jLabel31.setText("NsqMax");

        nsqMaxTextField.setText("jTextField1");
        nsqMaxTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                FileGenGUI.this.focusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                nsqMaxTextFieldFocusLost(evt);
            }
        });

        jLabel32.setText("Nmax");

        nMaxTextField.setText("jTextField1");
        nMaxTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                FileGenGUI.this.focusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                nMaxTextFieldFocusLost(evt);
            }
        });

        jLabel33.setText("Internal Degrees of Freedom");

        intDegFreedComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "On", "Off" }));
        intDegFreedComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                intDegFreedComboBoxActionPerformed(evt);
            }
        });

        jLabel34.setText("Intramolecular 1-5-Interactions");

        intraLjElComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "On", "Off" }));
        intraLjElComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                intraLjElComboBoxActionPerformed(evt);
            }
        });

        jLabel35.setText("Intramolecular 1-4-Interactions");

        ljElComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "On", "Off" }));
        ljElComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ljElComboBoxActionPerformed(evt);
            }
        });

        jLabel36.setText("Print IDF-Energies");

        qshaker.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Off", "On" }));
        qshaker.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                qshakerActionPerformed(evt);
            }
        });

        jLabel37.setText("Qshake");

        printidf.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "On", "Off" }));
        printidf.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                printidfActionPerformed(evt);
            }
        });

        jLabel38.setText("Tolerance");

        jLabel39.setText("Scale");

        tolerance.setText("0.0001");
        tolerance.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                tolerancefocusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                toleranceFocusLost(evt);
            }
        });
        tolerance.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                toleranceActionPerformed(evt);
            }
        });

        scale.setText("0");
        scale.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                scalefocusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                scaleFocusLost(evt);
            }
        });

        mvtlabel.setText("μVT-Steps");
        mvtlabel.setEnabled(false);

        mvttxtbox.setText("20000");
        mvttxtbox.setEnabled(false);
        mvttxtbox.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                mvttxtboxfocusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                mvttxtboxFocusLost(evt);
            }
        });

        rdflabel.setText("RDF");

        rdftxtbox.setText("0");
        rdftxtbox.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                rdftxtboxfocusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                rdftxtboxFocusLost(evt);
            }
        });

        numshelllabel.setText("Num Shells");

        numshelltxtbox.setText("200");
        numshelltxtbox.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                numshelltxtboxfocusGained(evt);
            }
            public void focusLost(java.awt.event.FocusEvent evt) {
                numshelltxtboxFocusLost(evt);
            }
        });

        jLabel27.setText("Version : ");

        versiontxt.setText("2.0");
        versiontxt.setEnabled(false);
        versiontxt.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                versiontxtFocusLost(evt);
            }
        });
        versiontxt.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                versiontxtActionPerformed(evt);
            }
        });

        jLabel40.setText("Special Features:");

        jLabel41.setText("Number of density bins");

        ndbtxt.setText("jTextField1");
        ndbtxt.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                ndbtxtFocusLost(evt);
            }
        });

        viewparbtn.setText("View par");
        viewparbtn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                viewparbtnActionPerformed(evt);
            }
        });

        viewresbtn.setText("View res");
        viewresbtn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                viewresbtnActionPerformed(evt);
            }
        });

        jLabel42.setText("WallTime");

        walltimetxt.setText("20160");
        walltimetxt.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                walltimetxtFocusLost(evt);
            }
        });
        walltimetxt.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                walltimetxtActionPerformed(evt);
            }
        });

        org.jdesktop.layout.GroupLayout generalPanelLayout = new org.jdesktop.layout.GroupLayout(generalPanel);
        generalPanel.setLayout(generalPanelLayout);
        generalPanelLayout.setHorizontalGroup(
            generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(generalPanelLayout.createSequentialGroup()
                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(jLabel40)
                    .add(generalPanelLayout.createSequentialGroup()
                        .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING, false)
                            .add(org.jdesktop.layout.GroupLayout.LEADING, generalPanelLayout.createSequentialGroup()
                                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                                    .add(jLabel33, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 221, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                                    .add(generalPanelLayout.createSequentialGroup()
                                        .add(84, 84, 84)
                                        .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING, false)
                                            .add(jLabel35, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                            .add(jLabel34, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                            .add(jLabel36, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                            .add(org.jdesktop.layout.GroupLayout.TRAILING, jLabel37, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))))
                                .add(27, 27, 27)
                                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING, false)
                                    .add(intDegFreedComboBox, 0, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .add(org.jdesktop.layout.GroupLayout.LEADING, ljElComboBox, 0, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .add(intraLjElComboBox, 0, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .add(qshaker, 0, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .add(printidf, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)))
                            .add(org.jdesktop.layout.GroupLayout.LEADING, generalPanelLayout.createSequentialGroup()
                                .add(jLabel41)
                                .add(73, 73, 73)
                                .add(ndbtxt)))
                        .add(18, 18, 18)
                        .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING, false)
                            .add(jLabel38, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .add(jLabel39, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING, false)
                            .add(tolerance, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 76, Short.MAX_VALUE)
                            .add(scale)))
                    .add(generalPanelLayout.createSequentialGroup()
                        .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                            .add(jLabel8)
                            .add(generalPanelLayout.createSequentialGroup()
                                .add(12, 12, 12)
                                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                                    .add(jLabel21, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 97, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                                    .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING, false)
                                        .add(org.jdesktop.layout.GroupLayout.LEADING, jLabel15, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .add(org.jdesktop.layout.GroupLayout.LEADING, mvtlabel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .add(org.jdesktop.layout.GroupLayout.LEADING, jLabel14, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .add(org.jdesktop.layout.GroupLayout.LEADING, jLabel13, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 100, Short.MAX_VALUE)))))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING, false)
                            .add(nvtstepsTextField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 93, Short.MAX_VALUE)
                            .add(nptstepsTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 1, Short.MAX_VALUE)
                            .add(org.jdesktop.layout.GroupLayout.LEADING, runstepsTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 1, Short.MAX_VALUE)
                            .add(mvttxtbox)
                            .add(typeofensembleComboBox, 0, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .add(org.jdesktop.layout.GroupLayout.LEADING, cutofftypeComboBox, 0, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .add(32, 32, 32)
                        .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                            .add(generalPanelLayout.createSequentialGroup()
                                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING, false)
                                    .add(org.jdesktop.layout.GroupLayout.LEADING, jLabel16, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .add(jLabel17, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .add(org.jdesktop.layout.GroupLayout.LEADING, jLabel18, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .add(rdflabel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .add(org.jdesktop.layout.GroupLayout.LEADING, numshelllabel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING, false)
                                    .add(org.jdesktop.layout.GroupLayout.LEADING, errorfreqTextField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 126, Short.MAX_VALUE)
                                    .add(org.jdesktop.layout.GroupLayout.LEADING, visualfreqTextField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 126, Short.MAX_VALUE)
                                    .add(org.jdesktop.layout.GroupLayout.LEADING, rdftxtbox, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 126, Short.MAX_VALUE)
                                    .add(org.jdesktop.layout.GroupLayout.LEADING, resultfreqTextField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 126, Short.MAX_VALUE)
                                    .add(org.jdesktop.layout.GroupLayout.LEADING, numshelltxtbox)))
                            .add(generalPanelLayout.createSequentialGroup()
                                .add(viewparbtn)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                .add(viewresbtn))))
                    .add(simPanel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(generalPanelLayout.createSequentialGroup()
                        .add(jLabel28)
                        .add(18, 18, 18)
                        .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING, false)
                            .add(generalPanelLayout.createSequentialGroup()
                                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                                    .add(jLabel32)
                                    .add(jLabel30, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 85, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                                .add(18, 18, 18)
                                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING, false)
                                    .add(nMaxTextField)
                                    .add(nVecMaxTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)))
                            .add(longRangeCorrComboBox, 0, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                    .add(generalPanelLayout.createSequentialGroup()
                        .addContainerGap()
                        .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                            .add(generalPanelLayout.createSequentialGroup()
                                .add(jLabel29)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.UNRELATED)
                                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                                    .add(kappaTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 69, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                                    .add(nsqMaxTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 69, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)))
                            .add(jLabel31, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 109, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)))
                    .add(generalPanelLayout.createSequentialGroup()
                        .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                            .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                                .add(jLabel1)
                                .add(jLabel27)
                                .add(org.jdesktop.layout.GroupLayout.TRAILING, jLabel2)
                                .add(org.jdesktop.layout.GroupLayout.TRAILING, jLabel4)
                                .add(org.jdesktop.layout.GroupLayout.TRAILING, jLabel5))
                            .add(jLabel42))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING)
                            .add(org.jdesktop.layout.GroupLayout.LEADING, walltimetxt, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 99, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                            .add(org.jdesktop.layout.GroupLayout.LEADING, generalPanelLayout.createSequentialGroup()
                                .add(massTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 69, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                .add(jLabel7))
                            .add(org.jdesktop.layout.GroupLayout.LEADING, generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING, false)
                                .add(generalPanelLayout.createSequentialGroup()
                                    .add(lengthTextField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                    .add(jLabel3))
                                .add(org.jdesktop.layout.GroupLayout.LEADING, systemofunitsComboBox, 0, 99, Short.MAX_VALUE)
                                .add(org.jdesktop.layout.GroupLayout.LEADING, versiontxt)
                                .add(org.jdesktop.layout.GroupLayout.LEADING, generalPanelLayout.createSequentialGroup()
                                    .add(energyTextField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 72, Short.MAX_VALUE)
                                    .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                    .add(jLabel6))))))
                .addContainerGap(198, Short.MAX_VALUE))
        );
        generalPanelLayout.setVerticalGroup(
            generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(generalPanelLayout.createSequentialGroup()
                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(versiontxt, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(jLabel27, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 17, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .add(8, 8, 8)
                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(walltimetxt, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(jLabel42, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 17, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(systemofunitsComboBox, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(jLabel1))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel2)
                    .add(lengthTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(jLabel3))
                .add(9, 9, 9)
                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel4)
                    .add(energyTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(jLabel6))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel5)
                    .add(massTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(jLabel7))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(simPanel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 203, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel8)
                    .add(typeofensembleComboBox, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(viewparbtn, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(viewresbtn, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(generalPanelLayout.createSequentialGroup()
                        .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                            .add(jLabel13)
                            .add(nvtstepsTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                            .add(jLabel14)
                            .add(nptstepsTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                            .add(jLabel17))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                            .add(jLabel18)
                            .add(mvttxtbox, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                            .add(mvtlabel))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                            .add(runstepsTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                            .add(jLabel15)))
                    .add(generalPanelLayout.createSequentialGroup()
                        .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                            .add(resultfreqTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                            .add(jLabel16))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(errorfreqTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(visualfreqTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                            .add(rdftxtbox, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                            .add(rdflabel))))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(numshelltxtbox, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                        .add(numshelllabel)
                        .add(cutofftypeComboBox, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                        .add(jLabel21)))
                .add(25, 25, 25)
                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel28)
                    .add(longRangeCorrComboBox, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(generalPanelLayout.createSequentialGroup()
                        .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                            .add(jLabel29)
                            .add(kappaTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                            .add(jLabel31)
                            .add(nsqMaxTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)))
                    .add(generalPanelLayout.createSequentialGroup()
                        .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                            .add(nMaxTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                            .add(jLabel32))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                            .add(nVecMaxTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                            .add(jLabel30))))
                .add(27, 27, 27)
                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel33)
                    .add(intDegFreedComboBox, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel36)
                    .add(printidf, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel37)
                    .add(qshaker, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(jLabel38)
                    .add(tolerance, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(intraLjElComboBox, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(jLabel34))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                        .add(ljElComboBox, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                        .add(jLabel39)
                        .add(scale, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                    .add(jLabel35))
                .add(18, 18, 18)
                .add(jLabel40)
                .add(18, 18, 18)
                .add(generalPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel41)
                    .add(ndbtxt, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(161, Short.MAX_VALUE))
        );

        leftjscroll.setViewportView(generalPanel);

        rightjscroll.setBorder(null);
        rightjscroll.setHorizontalScrollBarPolicy(javax.swing.ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
        rightjscroll.setVerticalScrollBarPolicy(javax.swing.ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER);
        rightjscroll.setPreferredSize(new java.awt.Dimension(600, 1030));

        ensemblesTabbedPane.setAutoscrolls(true);
        ensemblesTabbedPane.setMinimumSize(new java.awt.Dimension(100, 100));
        ensemblesTabbedPane.setPreferredSize(new java.awt.Dimension(600, 1030));
        rightjscroll.setViewportView(ensemblesTabbedPane);
        ensemblesTabbedPane.getAccessibleContext().setAccessibleParent(mainjpanel);

        org.jdesktop.layout.GroupLayout mainjpanelLayout = new org.jdesktop.layout.GroupLayout(mainjpanel);
        mainjpanel.setLayout(mainjpanelLayout);
        mainjpanelLayout.setHorizontalGroup(
            mainjpanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(mainjpanelLayout.createSequentialGroup()
                .addContainerGap()
                .add(mainjpanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(topjscroll)
                    .add(mainjpanelLayout.createSequentialGroup()
                        .add(leftjscroll, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 790, Short.MAX_VALUE)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(rightjscroll, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 857, Short.MAX_VALUE)))
                .addContainerGap())
        );
        mainjpanelLayout.setVerticalGroup(
            mainjpanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(mainjpanelLayout.createSequentialGroup()
                .addContainerGap()
                .add(topjscroll, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 56, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(mainjpanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING, false)
                    .add(leftjscroll, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 1211, Short.MAX_VALUE)
                    .add(rightjscroll, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap(org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        topjscroll.getAccessibleContext().setAccessibleParent(mainjscroll);
        leftjscroll.getAccessibleContext().setAccessibleParent(mainjscroll);
        rightjscroll.getAccessibleContext().setAccessibleParent(mainjscroll);

        mainjscroll.setViewportView(mainjpanel);

        org.jdesktop.layout.GroupLayout layout = new org.jdesktop.layout.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(org.jdesktop.layout.GroupLayout.TRAILING, mainjscroll, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 1705, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(mainjscroll, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 1315, Short.MAX_VALUE)
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents
    
    /**
     * Set the windows size according to the actual screen.
     */
    private void setwindowsize()
    {
        GraphicsDevice Gd = GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice();
        int screenwidth = Gd.getDisplayMode().getWidth();
        int screenheight = Gd.getDisplayMode().getHeight();
        JFrame parentframe = (JFrame) SwingUtilities.getRoot(this);
        //parentframe.setResizable(false);
        int Pheight = parentframe.getHeight();
        int Pwidth = parentframe.getWidth();
        
        if(screenwidth <= Pwidth && screenheight <= Pheight)
        {
         parentframe.setSize(screenwidth - 10, screenheight - 50);
         
        }
        else
        if(screenwidth <= Pwidth)
            {
             parentframe.setSize(screenwidth - 10, Pheight);
            
            
            }
            else
            if(screenheight < Pheight)
            {
               parentframe.setSize(Pwidth, screenheight - 50);
   
            }
        mainjscroll.getVerticalScrollBar().setUnitIncrement(16);
        leftjscroll.addMouseWheelListener(new MouseWheelListener() {
           
        @Override
        public void mouseWheelMoved(MouseWheelEvent e) {
        leftjscroll.getParent().dispatchEvent(e);
                }
        });
        rightjscroll.addMouseWheelListener(new MouseWheelListener() {
           
        @Override
        public void mouseWheelMoved(MouseWheelEvent e) {
        rightjscroll.getParent().dispatchEvent(e);
                }
        });
        topjscroll.addMouseWheelListener(new MouseWheelListener() {
           
        @Override
        public void mouseWheelMoved(MouseWheelEvent e) {
        topjscroll.getParent().dispatchEvent(e);
                }
        });
    }
    
    /**
     * 
     */
    private void intDegFreedChanged() {
        if (intDegFreedComboBox.getSelectedItem().toString().equals("On")){
            jLabel34.setEnabled(true);
            intraLjElComboBox.setEnabled(true);
            jLabel35.setEnabled(true);
            ljElComboBox.setEnabled(true);
            jLabel36.setEnabled(true);
            printidf.setEnabled(true);
            jLabel37.setEnabled(true);
            qshaker.setEnabled(true);
            jLabel38.setEnabled(true);
            tolerance.setEnabled(true);
            jLabel39.setEnabled(true);
            scale.setEnabled(true);
        }else{
            jLabel34.setEnabled(false);
            intraLjElComboBox.setEnabled(false);
            jLabel35.setEnabled(false);
            ljElComboBox.setEnabled(false);
            jLabel36.setEnabled(false);
            printidf.setEnabled(false);
            jLabel37.setEnabled(false);
            qshaker.setEnabled(false);
            jLabel38.setEnabled(false);
            tolerance.setEnabled(false);
            jLabel39.setEnabled(false);
            scale.setEnabled(false);
        }
        qshakerchanged();
        ljElchanged();
    }
    private void longRangeCorrChanged() {
        
        if (longRangeCorrComboBox.getSelectedItem().toString().equals("Ewald")){
            jLabel29.setEnabled(true);
            kappaTextField.setEnabled(true);
            jLabel30.setEnabled(true);
            nVecMaxTextField.setEnabled(true);
            jLabel31.setEnabled(true);
            nsqMaxTextField.setEnabled(true);
            jLabel32.setEnabled(true);
            nMaxTextField.setEnabled(true);
        }
        else if (longRangeCorrComboBox.getSelectedItem().toString().equals("Rodgers")){
            jLabel29.setEnabled(true);
            kappaTextField.setEnabled(true);
            jLabel30.setEnabled(false);
            nVecMaxTextField.setEnabled(false);
            jLabel31.setEnabled(false);
            nsqMaxTextField.setEnabled(false);
            jLabel32.setEnabled(false);
            nMaxTextField.setEnabled(false);
            
        }
        else{
            jLabel29.setEnabled(false);
            kappaTextField.setEnabled(false);
            jLabel30.setEnabled(false);
            nVecMaxTextField.setEnabled(false);
            jLabel31.setEnabled(false);
            nsqMaxTextField.setEnabled(false);
            jLabel32.setEnabled(false);
            nMaxTextField.setEnabled(false);
        }
    }
                                                               
    /**
     *
     * @param str
     * @return
     */
    public double validateInputDouble(String str){
        double ret = 0.0;
        if (str.length() > 0) {
            try {
                ret = Double.valueOf(str).doubleValue();
            } catch (Exception e){
                JOptionPane.showMessageDialog(new JFrame(), "The number you entered is not a valid decimal number.");
            }
        }
        return ret;
    }
    
    /**
     *
     * @param str
     * @return
     */
    public int validateInputInt(String str){
        int ret = 0;
        if (str.length() > 0) {
            try {
                ret = Integer.parseInt(str);
            } catch (Exception e){
                JOptionPane.showMessageDialog(new JFrame(), "The number you entered is not a valid number.");
            }
        }
        
        return ret;
    }
    
        
    private void saveButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_saveButtonActionPerformed

      
        checkMolFrac();
        String pmdirpathSelected;
        File pmfileSelected = null;
        File pmdirSelected;
        File pmdirDefault = null;
        
        class pmFilter extends FileFilter {
            public String getDescription() {
                return "Only Parameter Files";
            }
            
            public boolean accept(File f) {
                if (f == null)
                    return false;
                if (f.isDirectory())
                    return true;
                return f.getName().toLowerCase().endsWith(".par");
            }
        }
        
        JFileChooser pmfileChooser = new JFileChooser();
        pmfileChooser.addChoosableFileFilter(new pmFilter());
        pmfileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
        pmfileChooser.setCurrentDirectory(this.saveFile);
        //TODO File Temperature
        
        if (this.saveFile==null){
            /*String temperature = String.valueOf(((EnsemblePanel)this.ensemblesTabbedPane.getSelectedComponent()).eData.getTemperature());
            temperature = temperature.substring(0,temperature.indexOf("."));*/
            String f_name = this.filenameLabel.getText().substring(this.filenameLabel.getText().lastIndexOf(":") + 2);
            String tempFile = this.loadDir.toString();
            tempFile = tempFile+fileSep+f_name;
            this.saveFile = new File(tempFile);
        }
        if (this.importRes ==1 || this.typeofensembleComboBox.getSelectedItem().toString().equals("GE")) {
            String geFileStr = this.saveFile.toString().substring(0,this.saveFile.toString().length()-4);
            if (!geFileStr.contains("_ge"))
                geFileStr += "_ge";
            geFileStr += ".par";
            this.saveFile = new File(geFileStr);
        }
        pmfileChooser.setSelectedFile(this.saveFile);
        int returnVal = pmfileChooser.showSaveDialog(this);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            pmdirSelected = pmfileChooser.getCurrentDirectory();
            pmfileSelected = pmfileChooser.getSelectedFile();
            pmdirpathSelected = pmdirSelected.getAbsolutePath();
            
            try{
                /* Copy PmFiles */
                String copyPath = pmfileChooser.getSelectedFile().getAbsolutePath();
                copyPath = copyPath.substring(0,copyPath.lastIndexOf(this.fileSep));
                copyPmFiles(copyPath+this.fileSep,pmdirSelected);
                
                writeParFile(pmfileSelected);
                
                this.saveFile = pmfileSelected;
                JOptionPane.showMessageDialog(new JFrame(), pmfileChooser.getSelectedFile().toString() + " has been saved.");
                this.filenameLabel.setText(pmfileSelected.toString());
                
            } catch (Exception e) {
               JOptionPane.showMessageDialog(new JFrame(), "An unknown error occurred");
               // e.printStackTrace();
            }
            ((EnsemblePanel) this.ensemblesTabbedPane.getSelectedComponent()).updatePmButtons();
        } 
    }//GEN-LAST:event_saveButtonActionPerformed
    
    private void loadButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_loadButtonActionPerformed
        
        File pmdirDefault = null;
        
        class pmFilter extends FileFilter {
            public String getDescription() {
                return "Only Parameter Files";
            }
            
            public boolean accept(File f) {
                if (f == null)
                    return false;
                if (f.isDirectory())
                    return true;
                return f.getName().toLowerCase().endsWith(".par");
            }
        }
        
        JFileChooser pmfileChooser = new JFileChooser();
        pmfileChooser.addChoosableFileFilter(new pmFilter());
        pmfileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
        pmfileChooser.setCurrentDirectory(this.loadFile);
        int returnVal = pmfileChooser.showOpenDialog(this);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            
            this.loadFile = pmfileChooser.getSelectedFile();
            this.loadDir = pmfileChooser.getCurrentDirectory();
            this.filenameLabel.setText(this.loadFile.toString());
            
            try{
                readParFile(pmfileChooser.getSelectedFile());
                this.saveFile = this.loadFile;
                /*remove comment if res is finished */
                if (typeofensembleComboBox.getSelectedItem().toString().equals("GE")){
                    this.viewparbtn.setEnabled(true);
                    this.viewresbtn.setEnabled(true);
                }
                /* */
                
            } catch (IOException e) {
                JOptionPane.showMessageDialog(new JFrame(), "An unknown error occurred");
            }
        }
        
    }//GEN-LAST:event_loadButtonActionPerformed
    
    private void ensemblesremoveButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ensemblesremoveButtonActionPerformed
        
        removeEnsemble();
        
    }//GEN-LAST:event_ensemblesremoveButtonActionPerformed
    
    private void ensemblesaddButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ensemblesaddButtonActionPerformed
        
        addEnsemble();
        
    }//GEN-LAST:event_ensemblesaddButtonActionPerformed
                    
    private void integratorComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_integratorComboBoxActionPerformed
        integratorChanged();
    }//GEN-LAST:event_integratorComboBoxActionPerformed
    
    private void clearButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_clearButtonActionPerformed
        
        this.clearPanel();
        
        this.ensemblesTabbedPane.removeAll();
        this.ensemblesTabbedPane.add("Ensemble 1", new EnsemblePanel(this));
        this.importRes = 0;
        this.saveFile = null;
        this.loadFile = null;
        this.filenameLabel.setText("File: newname.par");
        
    }//GEN-LAST:event_clearButtonActionPerformed

    private void helpButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_helpButtonActionPerformed

    try{
         String filePath = UserDirectory+""+File.separator+"help"+File.separator+"ms2par_help_format.xlsx";
                Desktop.getDesktop().open(new File(filePath));
     // ViewFiles view = new ViewFiles("Help");
     // view.setVisible(true);
       // Desktop.getDesktop().open(new File(getClass().getClassLoader().getResource("resources\\ms2tools_help_format.xlsx").getPath()));
    }catch (Exception e){
      JOptionPane.showMessageDialog(null,"Exception caught ="+e.getMessage());
    }
        
        
}//GEN-LAST:event_helpButtonActionPerformed

    private void scaleFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_scaleFocusLost
        double value;
        try{
            value = Double.parseDouble(scale.getText());
        }
        catch(Exception e)
        {
            JOptionPane.showMessageDialog(null, "Scale should have the value between 0 and 1");
            value = 0;
            scale.setText(Double.toString(value));

        }
        if(value < 0 || value > 1)
        {
            JOptionPane.showMessageDialog(null, "Scale should have the value between 0 and 1");
            value = 0;
            scale.setText(Double.toString(value));
        }
    }//GEN-LAST:event_scaleFocusLost

    private void scalefocusGained(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_scalefocusGained
        scale.selectAll();
    }//GEN-LAST:event_scalefocusGained

    private void toleranceFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_toleranceFocusLost
        double value;
        try{
            value = Double.parseDouble(tolerance.getText());
        }
        catch(Exception e)
        {
            JOptionPane.showMessageDialog(null, "Tolerance should have the value equal or less than 0.0001");
            value = 0.0001;
            tolerance.setText(Double.toString(value));

        }
        if(value > 0.0001)
        {
            JOptionPane.showMessageDialog(null, "Tolerance should have the value equal or less than 0.0001");
            value = 0.0001;
            tolerance.setText(Double.toString(value));

        }
    }//GEN-LAST:event_toleranceFocusLost

    private void tolerancefocusGained(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_tolerancefocusGained
        tolerance.selectAll();
    }//GEN-LAST:event_tolerancefocusGained

    private void toleranceActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_toleranceActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_toleranceActionPerformed

    private void printidfActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_printidfActionPerformed

    }//GEN-LAST:event_printidfActionPerformed

    private void qshakerActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_qshakerActionPerformed
        qshakerchanged();
    }//GEN-LAST:event_qshakerActionPerformed

    private void ljElComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ljElComboBoxActionPerformed
        ljElchanged();
    }//GEN-LAST:event_ljElComboBoxActionPerformed

    private void intraLjElComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_intraLjElComboBoxActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_intraLjElComboBoxActionPerformed

    private void intDegFreedComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_intDegFreedComboBoxActionPerformed
        intDegFreedChanged();
    }//GEN-LAST:event_intDegFreedComboBoxActionPerformed

    private void nMaxTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_nMaxTextFieldFocusLost
        this.gData.setNMax(this.validateInputInt(this.nMaxTextField.getText()));
        this.nMaxTextField.setText(String.valueOf(this.gData.getNMax()));
        calcNVecMax();
    }//GEN-LAST:event_nMaxTextFieldFocusLost

    private void focusGained(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_focusGained
        JTextField jtf = (JTextField) evt.getSource();
        jtf.selectAll();
    }//GEN-LAST:event_focusGained

    private void nsqMaxTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_nsqMaxTextFieldFocusLost
        this.gData.setNsqMax(this.validateInputInt(this.nsqMaxTextField.getText()));
        this.nsqMaxTextField.setText(String.valueOf(this.gData.getNsqMax()));
        calcNVecMax();
    }//GEN-LAST:event_nsqMaxTextFieldFocusLost

    private void nVecMaxTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_nVecMaxTextFieldFocusLost
        this.gData.setNVecMax(this.validateInputInt(this.nVecMaxTextField.getText()));
        this.nVecMaxTextField.setText(String.valueOf(this.gData.getNVecMax()));
    }//GEN-LAST:event_nVecMaxTextFieldFocusLost

    private void kappaTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_kappaTextFieldFocusLost
        this.gData.setKappa(this.validateInputDouble(this.kappaTextField.getText()));
        this.kappaTextField.setText(String.valueOf(this.gData.getKappa()));
    }//GEN-LAST:event_kappaTextFieldFocusLost

    private void longRangeCorrComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_longRangeCorrComboBoxActionPerformed
        longRangeCorrChanged();
    }//GEN-LAST:event_longRangeCorrComboBoxActionPerformed

    private void typeofensembleComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_typeofensembleComboBoxActionPerformed
        typeOfEnsembleChanged();
        //this.simPanel.setPreferredSize(new Dimension(320,347));
    }//GEN-LAST:event_typeofensembleComboBoxActionPerformed

    private void cutofftypeComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cutofftypeComboBoxActionPerformed
        cuttoffTypeChanged();
    }//GEN-LAST:event_cutofftypeComboBoxActionPerformed

    private void visualfreqTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_visualfreqTextFieldFocusLost
        this.gData.setVisualFreq(this.validateInputInt(this.visualfreqTextField.getText()));
        this.visualfreqTextField.setText(String.valueOf(this.gData.getVisualFreq()));
    }//GEN-LAST:event_visualfreqTextFieldFocusLost

    private void errorfreqTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_errorfreqTextFieldFocusLost
        this.gData.setErrorFreq(this.validateInputInt(this.errorfreqTextField.getText()));
        this.errorfreqTextField.setText(String.valueOf(this.gData.getErrorFreq()));
    }//GEN-LAST:event_errorfreqTextFieldFocusLost

    private void resultfreqTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_resultfreqTextFieldFocusLost
        this.gData.setResultFreq(this.validateInputInt(this.resultfreqTextField.getText()));
        this.resultfreqTextField.setText(String.valueOf(this.gData.getResultFreq()));
    }//GEN-LAST:event_resultfreqTextFieldFocusLost

    private void runstepsTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_runstepsTextFieldFocusLost
        this.gData.setRunSteps(this.validateInputInt(this.runstepsTextField.getText()));
        this.runstepsTextField.setText(String.valueOf(this.gData.getRunSteps()));
    }//GEN-LAST:event_runstepsTextFieldFocusLost

    private void nptstepsTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_nptstepsTextFieldFocusLost
        this.gData.setNptSteps(this.validateInputInt(this.nptstepsTextField.getText()));
        this.nptstepsTextField.setText(String.valueOf(this.gData.getNptSteps()));
    }//GEN-LAST:event_nptstepsTextFieldFocusLost

    private void nvtstepsTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_nvtstepsTextFieldFocusLost
        this.gData.setNvtSteps(this.validateInputInt(this.nvtstepsTextField.getText()));
        this.nvtstepsTextField.setText(String.valueOf(this.gData.getNvtSteps()));
    }//GEN-LAST:event_nvtstepsTextFieldFocusLost

    private void rMinRadiusTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_rMinRadiusTextFieldFocusLost
        this.gData.setRMinRadius(this.validateInputDouble(this.rMinRadiusTextField.getText()));
        this.rMinRadiusTextField.setText(String.valueOf(this.gData.getRMinRadius()));
    }//GEN-LAST:event_rMinRadiusTextFieldFocusLost

    private void rMaxRadiusTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_rMaxRadiusTextFieldFocusLost
        this.gData.setRMaxRadius(this.validateInputDouble(this.rMaxRadiusTextField.getText()));
        this.rMaxRadiusTextField.setText(String.valueOf(this.gData.getRMaxRadius()));
    }//GEN-LAST:event_rMaxRadiusTextFieldFocusLost

    private void rStepsTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_rStepsTextFieldFocusLost
        this.gData.setRSteps(this.validateInputInt(this.rStepsTextField.getText()));
        this.rStepsTextField.setText(String.valueOf(this.gData.getRSteps()));
    }//GEN-LAST:event_rStepsTextFieldFocusLost

    private void nOrientTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_nOrientTextFieldFocusLost
        this.gData.setNOrient(this.validateInputInt(this.nOrientTextField.getText()));
        this.nOrientTextField.setText(String.valueOf(this.gData.getNOrient()));
    }//GEN-LAST:event_nOrientTextFieldFocusLost

    private void mcorStepsTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_mcorStepsTextFieldFocusLost
        this.gData.setMcorSteps(this.validateInputInt(this.mcorStepsTextField.getText()));
        this.mcorStepsTextField.setText(String.valueOf(this.gData.getMcorSteps()));
        if (this.gData.getMcorSteps()>0){
            this.acceptanceTextField.setEnabled(true);
            this.jLabel12.setEnabled(true);
        } else {
            this.acceptanceTextField.setEnabled(false);
            this.jLabel12.setEnabled(false);
        }
    }//GEN-LAST:event_mcorStepsTextFieldFocusLost

    private void acceptanceTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_acceptanceTextFieldFocusLost
        this.gData.setAcceptance(this.validateInputDouble(this.acceptanceTextField.getText()));
        this.acceptanceTextField.setText(String.valueOf(this.gData.getAcceptance()));
    }//GEN-LAST:event_acceptanceTextFieldFocusLost

    private void timestepTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_timestepTextFieldFocusLost
        this.gData.setTimeStep(this.validateInputDouble(this.timestepTextField.getText()));
        if(systemofunitsComboBox.getSelectedItem().toString().equals("Reduced"))
            this.timestepTextField.setText(String.valueOf(this.gData.getReducedTimeStep(this.gData.getTimeStep())));
        else
           this.timestepTextField.setText(String.valueOf(this.gData.getTimeStep())); 
    }//GEN-LAST:event_timestepTextFieldFocusLost

    private void typeofsimulationComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_typeofsimulationComboBoxActionPerformed
        typeOfSimulationChanged();
        //this.simPanel.setPreferredSize(new Dimension(320,347));
    }//GEN-LAST:event_typeofsimulationComboBoxActionPerformed

    private void massTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_massTextFieldFocusLost
        this.gData.setMassUnit(this.validateInputDouble(this.massTextField.getText()));
        this.massTextField.setText(String.valueOf(this.gData.getMassUnit()));
    }//GEN-LAST:event_massTextFieldFocusLost

    private void energyTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_energyTextFieldFocusLost
        this.gData.setEnergyUnit(this.validateInputDouble(this.energyTextField.getText()));
        this.energyTextField.setText(String.valueOf(this.gData.getEnergyUnit()));
    }//GEN-LAST:event_energyTextFieldFocusLost

    private void lengthTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_lengthTextFieldFocusLost
        this.gData.setLengthUnit(this.validateInputDouble(this.lengthTextField.getText()));
        this.lengthTextField.setText(String.valueOf(this.gData.getLengthUnit()));
    }//GEN-LAST:event_lengthTextFieldFocusLost

    private void lengthTextFieldActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_lengthTextFieldActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_lengthTextFieldActionPerformed

    private void systemofunitsComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_systemofunitsComboBoxActionPerformed
        systemOfUnitsChanged();
    }//GEN-LAST:event_systemofunitsComboBoxActionPerformed

    private void mvttxtboxfocusGained(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_mvttxtboxfocusGained
        // TODO add your handling code here:
    }//GEN-LAST:event_mvttxtboxfocusGained

    private void mvttxtboxFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_mvttxtboxFocusLost
        this.gData.setmuevttSteps(this.validateInputInt(this.mvttxtbox.getText()));
        this.mvttxtbox.setText(String.valueOf(this.gData.getmuevttSteps()));
    }//GEN-LAST:event_mvttxtboxFocusLost

    private void rdftxtboxfocusGained(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_rdftxtboxfocusGained
        // TODO add your handling code here:
    }//GEN-LAST:event_rdftxtboxfocusGained

    private void rdftxtboxFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_rdftxtboxFocusLost
        this.gData.setrdffreq(this.validateInputInt(this.rdftxtbox.getText()));
        this.rdftxtbox.setText(String.valueOf(this.gData.getrdffreq()));
    }//GEN-LAST:event_rdftxtboxFocusLost

    private void numshelltxtboxfocusGained(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_numshelltxtboxfocusGained
        // TODO add your handling code here:
    }//GEN-LAST:event_numshelltxtboxfocusGained

    private void numshelltxtboxFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_numshelltxtboxFocusLost
        this.gData.setnumshell(this.validateInputInt(this.numshelltxtbox.getText()));
        this.numshelltxtbox.setText(String.valueOf(this.gData.getnumshell()));
    }//GEN-LAST:event_numshelltxtboxFocusLost

    private void optPressureComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_optPressureComboBoxActionPerformed
        this.eData.setOptPressure(optPressureComboBox.getSelectedItem().toString());
    }//GEN-LAST:event_optPressureComboBoxActionPerformed

    private void commonequiComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_commonequiComboBoxActionPerformed
       this.eData.setcommonequi(commonequiComboBox.getSelectedItem().toString());
    }//GEN-LAST:event_commonequiComboBoxActionPerformed

    private void viewparbtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_viewparbtnActionPerformed
        // TODO add your handling code here:
        
       ViewFiles view = new ViewFiles(TrimFileName() + ".par");
       view.setVisible(true);
     
    }//GEN-LAST:event_viewparbtnActionPerformed

    private void viewresbtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_viewresbtnActionPerformed
        // TODO add your handling code here:
        File file = new File(TrimFileName() + ".res");
        if(file.exists())
        {
        ViewFiles view = new ViewFiles(TrimFileName() + ".res");
        view.setVisible(true);
        }
        else
        {
            JOptionPane.showMessageDialog(null,"Error: No .res file found");
        }
    }//GEN-LAST:event_viewresbtnActionPerformed

    private void versiontxtActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_versiontxtActionPerformed
        // TODO add your handling code here:
        
    }//GEN-LAST:event_versiontxtActionPerformed

    private void versiontxtFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_versiontxtFocusLost
         this.versiontxt.setText(String.valueOf(this.validateInputDouble(this.versiontxt.getText())));
    }//GEN-LAST:event_versiontxtFocusLost

    private void nvtstepsTextFieldActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_nvtstepsTextFieldActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_nvtstepsTextFieldActionPerformed

    private void ndbtxtFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_ndbtxtFocusLost
        this.gData.setNumDinBins(this.validateInputInt(this.ndbtxt.getText()));
        this.ndbtxt.setText(String.valueOf(this.gData.getNumDinBins()));
    }//GEN-LAST:event_ndbtxtFocusLost

    private void walltimetxtFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_walltimetxtFocusLost
        this.gData.setWallTime(this.validateInputInt(this.walltimetxt.getText()));
        this.walltimetxt.setText(String.valueOf(this.gData.getWallTime()));
    }//GEN-LAST:event_walltimetxtFocusLost

    private void walltimetxtActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_walltimetxtActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_walltimetxtActionPerformed

    private String TrimFileName()
    {
        String path = this.loadFile.toString();
        path = path.substring(0, path.length() - 4);
        return path;
    }
    private void ljElchanged()
    {
        if(ljElComboBox.getSelectedItem().toString().equals("On"))
        {
            jLabel39.setEnabled(true);
            scale.setEnabled(true);
        }
        else
        {
            jLabel39.setEnabled(false);
            scale.setEnabled(false);
        }
    }    
    private void qshakerchanged()
    {
        if(qshaker.getSelectedItem().toString().equals("On"))
        {
            jLabel38.setEnabled(true);
            tolerance.setEnabled(true);
            jLabel10.setEnabled(false);
             this.gData.setIntegrator("Leapfrog");
             this.integratorComboBox.setSelectedItem("Leapfrog");
            this.integratorComboBox.setEnabled(false);
        }
        else
        {
            jLabel38.setEnabled(false);
            tolerance.setEnabled(false);
            jLabel10.setEnabled(true);
            this.gData.setIntegrator("Gear");
            this.integratorComboBox.setEnabled(true);
        }
    }
    
    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new FileGenGUI().setVisible(true);
            }
        });
    }
    
    
    /**
     * call default parameters for form
     */
    private void initPanel() {
        
        UserDirectory = System.getProperty("user.dir");
        this.gData = new General();
        this.eData = new Ensemble();
        updatePanel();
        this.filenameLabel.setText("File: newname.par");
        viewparbtn.setEnabled(false);
        viewresbtn.setEnabled(false); 
 

        
}
/*
    Add key short cuts for button actions e.g save, help, load
    */
private void addShortCutsForButtonActions() {
    //Short cut for save button
    Action savebuttonAction = new AbstractAction("Save") {
        @Override
        public void actionPerformed(ActionEvent evt) {
            saveButtonActionPerformed(evt);
        }
    };
 
    this.saveButton.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(
        KeyStroke.getKeyStroke(KeyEvent.VK_S, KeyEvent.CTRL_DOWN_MASK), "Save");
    this.saveButton.getActionMap().put("Save", savebuttonAction);

    //Shortcut for load button
    Action loadbuttonAction = new AbstractAction("Load") {
        @Override
        public void actionPerformed(ActionEvent evt) {
           loadButtonActionPerformed(evt);
        }
    };
 
    this.loadButton.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(
        KeyStroke.getKeyStroke(KeyEvent.VK_A, KeyEvent.CTRL_DOWN_MASK), "Load");
    this.loadButton.getActionMap().put("Load", loadbuttonAction);
    
   //shortcut for help button
    Action helpbuttonAction = new AbstractAction("Help") {

        @Override
        public void actionPerformed(ActionEvent evt) {
            helpButtonActionPerformed(evt);
        }
    };

    this.helpButton.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(
            KeyStroke.getKeyStroke(KeyEvent.VK_F1, 0), "Help");
    this.helpButton.getActionMap().put("Help", helpbuttonAction);

}
     private void settimeunitsteps()
    {
        double version = Double.parseDouble(versiontxt.getText());
        if(this.systemofunitsComboBox.getSelectedItem().toString().equals("SI") && version < 2.0)
            {
                JOptionPane.showMessageDialog(this, "System of Unit 'SI' detected but version is less then 2.0.\n\n"
                        + "Now program will update values to a new version.","Version Warning", JOptionPane.WARNING_MESSAGE);
                this.versiontxt.setText("2.0");
                this.systemofunitsComboBox.setSelectedItem("SI");
                    this.gData.setSystemOfUnits("SI");
                this.timestepunitlabel.setText("femtosec");
                    this.gData.setTimeStepUnit("femtosec");
                double timestep_new = this.gData.getTimeStep() * 1096.6877 * Math.sqrt(this.gData.getMassUnit() / this.gData.getEnergyUnit()) * this.gData.getLengthUnit();
                timestep_new = Math.floor(timestep_new);
                this.gData.setTimeStep(timestep_new);
                timestepTextField.setText(String.valueOf(this.gData.getTimeStep()));// recalculate values
                java.util.List ComponentList = Arrays.asList(this.ensemblesTabbedPane.getComponents());
                Iterator CIt = ComponentList.iterator();
                EnsemblePanel curE;
                System.out.println("now ensemble");
                while (CIt.hasNext()) {
                    curE = (EnsemblePanel)CIt.next();
                    double pistonmass_new = curE.eData.getPistonMass() * 1.66054 * Math.pow(10, 13) * (this.gData.getMassUnit() / this.gData.getLengthUnit());
                    curE.eData.setPistonMass(pistonmass_new);
                    System.out.println("this is :"+String.valueOf(pistonmass_new));
                    curE.pistonmassTextField.setText(String.valueOf(pistonmass_new));
                }
            }
            else
            if(systemofunitsComboBox.getSelectedItem().toString().equals("SI") && version >= 2.0)
            {
                if(timestepunitlabel.getText().equals("reduced")) {
                    double timestep_new = this.gData.getTimeStep() * 1096.6877 * Math.sqrt(this.gData.getMassUnit() / this.gData.getEnergyUnit()) * this.gData.getLengthUnit();
                    timestep_new = Math.floor(timestep_new);
                    this.gData.setTimeStep(timestep_new);
                }
                
                timestepunitlabel.setText("femtosec");
                this.gData.setTimeStepUnit("femtosec");
                timestepTextField.setText(String.valueOf(this.gData.getTimeStep()));
            }
            else
            if(systemofunitsComboBox.getSelectedItem().toString().equals("Reduced")) {
                this.timestepunitlabel.setText("reduced");
                this.gData.setTimeStepUnit("reduced");
                timestepTextField.setText(String.valueOf(this.gData.getReducedTimeStep(this.gData.getTimeStep())));
            }
    }
    private void systemOfUnitsChanged() {
        
        EnsemblePanel curE;
        this.gData.setSystemOfUnits(this.systemofunitsComboBox.getSelectedItem().toString());
        if (systemofunitsComboBox.getSelectedItem().toString().equals("SI")){
            jLabel2.setEnabled(true);
            lengthTextField.setEnabled(true);
            jLabel4.setEnabled(true);
            energyTextField.setEnabled(true);
            jLabel5.setEnabled(true);
            massTextField.setEnabled(true);
            jLabel3.setEnabled(true);
            jLabel6.setEnabled(true);
            jLabel7.setEnabled(true);

            //timestepunitComboBox.setSelectedItem("femtosec");
            
            java.util.List ComponentList = Arrays.asList(this.ensemblesTabbedPane.getComponents());
            Iterator CIt = ComponentList.iterator();
            while (CIt.hasNext()) {
                curE = (EnsemblePanel)CIt.next();
                curE.jLabel26.setVisible(false);
                curE.temperatureunitComboBox.setVisible(true);
                curE.jLabel31.setVisible(false);
                curE.pressureunitComboBox.setVisible(true);
                curE.pressureunitComboBox.setEnabled(false);                
                curE.jLabel32.setVisible(false);
                curE.densityunitComboBox.setVisible(true);
                curE.jLabel6.setVisible(false);
                curE.vapdensityunitComboBox.setVisible(true);
                curE.vapdensityunitComboBox.setEnabled(true);
                curE.jLabel7.setVisible(false);
                curE.pistonmassunitComboBox.setVisible(true);
                curE.pistonmassunitComboBox.setEnabled(false);
                if (this.typeofensembleComboBox.getSelectedItem().equals("NPT")){
                       curE.pressureunitComboBox.setEnabled(true);
                       if (this.typeofsimulationComboBox.getSelectedItem().equals("MD")){
                          curE.pistonmassunitComboBox.setEnabled(true);
                       }
                       else {
                          curE.pistonmassunitComboBox.setEnabled(false);
                       }
                }
                else if (this.typeofensembleComboBox.getSelectedItem().equals("NPH")){
                       curE.pressureunitComboBox.setEnabled(true);
                       if (this.typeofsimulationComboBox.getSelectedItem().equals("MD")){
                          curE.pistonmassunitComboBox.setEnabled(true);
                       }
                       else {
                          curE.pistonmassunitComboBox.setEnabled(false);
                       }
                }

            }
        } else if (systemofunitsComboBox.getSelectedItem().toString().equals("Reduced")){
            jLabel2.setEnabled(true);
            lengthTextField.setEnabled(true);
            jLabel4.setEnabled(true);
            energyTextField.setEnabled(true);
            jLabel5.setEnabled(true);
            jLabel3.setEnabled(true);
            jLabel6.setEnabled(true);
            jLabel7.setEnabled(true);
            massTextField.setEnabled(true);
//            timestepunitComboBox.setSelectedItem("reduced");
            java.util.List ComponentList = Arrays.asList(this.ensemblesTabbedPane.getComponents());
            Iterator CIt = ComponentList.iterator();
            while (CIt.hasNext()) {
                curE = (EnsemblePanel)CIt.next();
                curE.jLabel26.setVisible(true);
                curE.temperatureunitComboBox.setVisible(false);
                curE.jLabel31.setVisible(true);
                curE.pressureunitComboBox.setVisible(false);
                curE.jLabel31.setEnabled(false);
                curE.jLabel32.setVisible(true);
                curE.densityunitComboBox.setVisible(false);
                curE.jLabel6.setVisible(true);
                curE.vapdensityunitComboBox.setVisible(false);
                curE.jLabel7.setVisible(true);
                curE.jLabel7.setEnabled(false);
                if (this.typeofensembleComboBox.getSelectedItem().equals("NPT")){
                       curE.jLabel31.setEnabled(true);
                       curE.jLabel7.setEnabled(false);
                       if (this.typeofsimulationComboBox.getSelectedItem().equals("MD")){
                           curE.jLabel7.setEnabled(true);
                       }
                }
                else if (this.typeofensembleComboBox.getSelectedItem().equals("NPH")){
                       curE.jLabel31.setEnabled(true);
                       curE.jLabel7.setEnabled(false);
                       if (this.typeofsimulationComboBox.getSelectedItem().equals("MD")){
                           curE.jLabel7.setEnabled(true);
                       }
                }
                curE.pistonmassunitComboBox.setVisible(false);
            }
        }
        
        settimeunitsteps();
        
    }
    
    private void typeOfEnsembleChanged() {
        /* NPT SELECTION */
        EnsemblePanel curE;
        if (typeofensembleComboBox.getSelectedItem().toString().equals("NPT")){
            jLabel14.setEnabled(true);
            nptstepsTextField.setEnabled(true);
            mvtlabel.setEnabled(false);
            mvttxtbox.setEnabled(false);
            
            
            java.util.List ComponentList = Arrays.asList(this.ensemblesTabbedPane.getComponents());
            Iterator CIt = ComponentList.iterator();
            while (CIt.hasNext()) {
                curE = (EnsemblePanel)CIt.next();
                if (typeofsimulationComboBox.getSelectedItem().toString().equals("MD")){
   /*ich*/      curE.jLabel46.setEnabled(true);
                curE.transportComboBox.setEnabled(true);  
            } else {
   /*ich*/      curE.jLabel46.setEnabled(false);
                curE.transportComboBox.setEnabled(false);
            }
                curE.jLabel23.setEnabled(true);      
                curE.pressureTextField.setEnabled(true);
                 if (this.systemofunitsComboBox.getSelectedItem().toString().equals("SI")){
                   curE.pressureunitComboBox.setEnabled(true);
                   curE.jLabel31.setEnabled(false);
                } else {
                   curE.jLabel31.setEnabled(true);
                   curE.pressureunitComboBox.setEnabled(false);
                }
                
               // optpressurelabel.setEnabled(true);
                //optPressureComboBox.setEnabled(true);
                //commonequilabel.setEnabled(true);
                //commonequiComboBox.setEditable(true);
                curE.hamiltonianlbl.setEnabled(false);
                curE.hamitoniantxt.setEnabled(false);
                curE.enthalphylbl.setEnabled(false);
                curE.enthalphytxt.setEnabled(false);
                
                if (curE.transportComboBox.getSelectedItem().toString().equals("On")){ 
                    curE.jLabel40.setEnabled(true);
                    curE.corrLengthCFTextField.setEnabled(true);
                    curE.jLabel41.setEnabled(true);
                    curE.resFreqCFTextField.setEnabled(true);
                    curE.jLabel42.setEnabled(true);
                    curE.spanFunCFTextField.setEnabled(true);
                    curE.jLabel44.setEnabled(true);
                    curE.stepcorrfuntxt.setEnabled(true);
                    curE.jLabel43.setEnabled(true);
                    curE.jLabel45.setEnabled(true);
                    curE.viewFunCFTextField.setEnabled(true);
                }else{
                    curE.jLabel40.setEnabled(false);
                    curE.corrLengthCFTextField.setEnabled(false);
                    curE.jLabel41.setEnabled(false);
                    curE.resFreqCFTextField.setEnabled(false);
                    curE.jLabel42.setEnabled(false);
                    curE.spanFunCFTextField.setEnabled(false);
                    curE.jLabel44.setEnabled(false);
                    curE.stepcorrfuntxt.setEnabled(false);
                    curE.jLabel43.setEnabled(false);
                    curE.jLabel45.setEnabled(false);
                    curE.viewFunCFTextField.setEnabled(false); 
                }
                    
                if (typeofsimulationComboBox.getSelectedItem().toString().equals("MC")){
                    curE.jLabel46.setEnabled(false);
                    curE.transportComboBox.setEnabled(false);
                    curE.jLabel27.setEnabled(false);
                    curE.pistonmassTextField.setEnabled(false);
                    curE.pistonmassunitComboBox.setEnabled(false);
                } else if (typeofsimulationComboBox.getSelectedItem().toString().equals("MD")){
                   /* optpressure.setEnabled(false);
                    curE.optPressureComboBox.setEnabled(false);
                    if (curE.optPressureComboBox.getSelectedItem().toString().equals("No")){
                        curE.eData.setOptPressure("Yes");
                        curE.optPressureComboBox.setSelectedItem(curE.eData.getOptPressure());
                       // curE.optPressureChanged();
                    }*/
                    
                    curE.jLabel46.setEnabled(true);
                    curE.transportComboBox.setEnabled(true);
                    curE.jLabel27.setEnabled(true);
                    curE.pistonmassTextField.setEnabled(true);
                    if (this.systemofunitsComboBox.getSelectedItem().toString().equals("SI")){
                        curE.pistonmassunitComboBox.setEnabled(true);
                        curE.jLabel7.setEnabled(false);
                    }
                    else {
                        curE.pistonmassunitComboBox.setEnabled(false);
                        curE.jLabel7.setEnabled(true);
                    }
                }
                else{
/*                    curE.jLabel8.setEnabled(false);
                    curE.optPressureComboBox.setEnabled(false);
                    if (curE.optPressureComboBox.getSelectedItem().toString().equals("No")){
                        curE.eData.setOptPressure("Yes");
                        curE.optPressureComboBox.setSelectedItem(curE.eData.getOptPressure());
                      //  curE.optPressureChanged();
                    } */
                }
            }
        }
        else if (typeofensembleComboBox.getSelectedItem().toString().equals("NPH")){
            jLabel14.setEnabled(true);
            nptstepsTextField.setEnabled(true);
             mvtlabel.setEnabled(false);
            mvttxtbox.setEnabled(false);
 /*ich*/     
          
            
            java.util.List ComponentList = Arrays.asList(this.ensemblesTabbedPane.getComponents());
            Iterator CIt = ComponentList.iterator();
            while (CIt.hasNext()) {
                curE = (EnsemblePanel)CIt.next();
                curE.jLabel46.setEnabled(false);
                curE.transportComboBox.setEnabled(false);
                curE.jLabel23.setEnabled(true);
                curE.pressureTextField.setEnabled(true);
                curE.enthalphylbl.setEnabled(true);
                curE.enthalphytxt.setEnabled(true);
                curE.hamiltonianlbl.setEnabled(false);
                curE.hamitoniantxt.setEnabled(false);
                if (this.systemofunitsComboBox.getSelectedItem().toString().equals("SI")){
                   curE.pressureunitComboBox.setEnabled(true);
                   curE.jLabel31.setEnabled(false);
                } else {
                   curE.jLabel31.setEnabled(true);
                   curE.pressureunitComboBox.setEnabled(false);
                }

                
/*                curE.jLabel8.setEnabled(false);
                curE.optPressureComboBox.setEnabled(false);
                if (curE.optPressureComboBox.getSelectedItem().toString().equals("No")){
                        curE.eData.setOptPressure("Yes");
                        curE.optPressureComboBox.setSelectedItem(curE.eData.getOptPressure());
                      //  curE.optPressureChanged();
                } */
                if (curE.transportComboBox.getSelectedItem().toString().equals("On")){
                curE.eData.setTransport("Off");
                curE.transportComboBox.setSelectedItem(curE.eData.getTransport());
                transportChanged(curE);
            }
                    curE.jLabel40.setEnabled(false);
                    curE.corrLengthCFTextField.setEnabled(false);
                   curE.jLabel41.setEnabled(false);
                    curE.resFreqCFTextField.setEnabled(false);
                    curE.jLabel42.setEnabled(false);
                    curE.spanFunCFTextField.setEnabled(false);
                    curE.jLabel43.setEnabled(false);
                    curE.jLabel44.setEnabled(false);
                    curE.stepcorrfuntxt.setEnabled(false);
                    curE.jLabel43.setEnabled(false);
                    curE.jLabel45.setEnabled(false);
                    curE.viewFunCFTextField.setEnabled(false);               
                
                if (typeofsimulationComboBox.getSelectedItem().toString().equals("MC")){
                    //this.gData.setTypeOfSimulation("MD");
                    //typeofsimulationComboBox.setSelectedItem(this.gData.getTypeOfSimulation());
                    //typeOfSimulationChanged();
                } else if (typeofsimulationComboBox.getSelectedItem().toString().equals("MD")){
                    curE.jLabel46.setEnabled(true);
                    curE.transportComboBox.setEnabled(true);
                    curE.jLabel27.setEnabled(true);
                    curE.pistonmassTextField.setEnabled(true);
                    if (this.systemofunitsComboBox.getSelectedItem().toString().equals("SI")){
                        curE.pistonmassunitComboBox.setEnabled(true);
                        curE.jLabel7.setEnabled(false);
                    }
                    else {
                        curE.pistonmassunitComboBox.setEnabled(false);
                        curE.jLabel7.setEnabled(true);
                    }

                }
            }
        }
//        else if (typeofensembleComboBox.getSelectedItem().toString().equals("NVE"))
//        {
//            java.util.List ComponentList = Arrays.asList(this.ensemblesTabbedPane.getComponents());
//            Iterator CIt = ComponentList.iterator();
//            while (CIt.hasNext()) {
//                curE = (EnsemblePanel)CIt.next();
//                curE.hamiltonianlbl.setEnabled(true);
//                curE.hamitoniantxt.setEnabled(true);
//            }
//        }
        else {
            jLabel14.setEnabled(false);
            nptstepsTextField.setEnabled(false);
            mvtlabel.setEnabled(false);
            mvttxtbox.setEnabled(false);
             
            
            if (typeofensembleComboBox.getSelectedItem().toString().equals("GE")){
                jLabel14.setEnabled(true);
                nptstepsTextField.setEnabled(true);
                jLabel14.setEnabled(true);
                nptstepsTextField.setEnabled(true);
     /*ich*/      
                mvtlabel.setEnabled(true);
                mvttxtbox.setEnabled(true);
                this.gData.setTypeOfSimulation("MC");
                typeofsimulationComboBox.setSelectedItem(this.gData.getTypeOfSimulation());
                typeOfSimulationChanged();
            }
            
            java.util.List ComponentList = Arrays.asList(this.ensemblesTabbedPane.getComponents());
            Iterator CIt = ComponentList.iterator();
            while (CIt.hasNext()) {
                curE = (EnsemblePanel)CIt.next();
                if (typeofensembleComboBox.getSelectedItem().toString().equals("GE")){
                    curE.jLabel46.setEnabled(false);
                    curE.transportComboBox.setEnabled(false);
                }
                if (typeofsimulationComboBox.getSelectedItem().toString().equals("MD")){
   /*ich*/      curE.jLabel46.setEnabled(true);
                curE.transportComboBox.setEnabled(true);  
            } else {
   /*ich*/      curE.jLabel46.setEnabled(false);
                curE.transportComboBox.setEnabled(false);
            }
                curE.jLabel23.setEnabled(false);
                curE.pressureTextField.setEnabled(false);
                curE.pressureunitComboBox.setEnabled(false);
                curE.jLabel31.setEnabled(false);
                curE.hamiltonianlbl.setEnabled(false);
                curE.hamitoniantxt.setEnabled(false);
                curE.enthalphylbl.setEnabled(false);
                curE.enthalphytxt.setEnabled(false);
                
                if (typeofensembleComboBox.getSelectedItem().toString().equals("NVE"))
                {
                        curE.hamiltonianlbl.setEnabled(true);
                        curE.hamitoniantxt.setEnabled(true);
                    
                }
/*                curE.jLabel8.setEnabled(false);
                curE.optPressureComboBox.setEnabled(false);
                if (curE.optPressureComboBox.getSelectedItem().toString().equals("No")){
                        curE.eData.setOptPressure("Yes");
                        curE.optPressureComboBox.setSelectedItem(curE.eData.getOptPressure());
                      //  curE.optPressureChanged();
                }*/
                
                if (curE.transportComboBox.getSelectedItem().toString().equals("On")){ 
                    curE.jLabel40.setEnabled(true);
                    curE.corrLengthCFTextField.setEnabled(true);
                    curE.jLabel41.setEnabled(true);
                    curE.resFreqCFTextField.setEnabled(true);
                    curE.jLabel42.setEnabled(true);
                    curE.spanFunCFTextField.setEnabled(true);
                    curE.jLabel43.setEnabled(true);
                    curE.jLabel44.setEnabled(true);
                    curE.stepcorrfuntxt.setEnabled(true);
                    curE.jLabel43.setEnabled(true);
                    curE.jLabel45.setEnabled(true);
                    curE.viewFunCFTextField.setEnabled(true);
                }else{
                    curE.jLabel40.setEnabled(false);
                    curE.corrLengthCFTextField.setEnabled(false);
                    curE.jLabel41.setEnabled(false);
                    curE.resFreqCFTextField.setEnabled(false);
                    curE.jLabel42.setEnabled(false);
                    curE.spanFunCFTextField.setEnabled(false);
                    curE.jLabel44.setEnabled(false);
                    curE.stepcorrfuntxt.setEnabled(false);
                    curE.jLabel43.setEnabled(false);
                    curE.jLabel45.setEnabled(false);
                    curE.viewFunCFTextField.setEnabled(false);
                }               
                
                curE.jLabel46.setEnabled(false);
                curE.transportComboBox.setEnabled(false);
                curE.jLabel27.setEnabled(false);
                curE.pistonmassTextField.setEnabled(false);
                curE.pistonmassunitComboBox.setEnabled(false);
                curE.jLabel7.setEnabled(false);
            }
        }
        /* GE SELECTION */
        java.util.List ComponentList = Arrays.asList(this.ensemblesTabbedPane.getComponents());
        Iterator CIt = ComponentList.iterator();
        while (CIt.hasNext()) {
            curE = (EnsemblePanel)CIt.next();
            curE.updatePmButtons();
            if (typeofensembleComboBox.getSelectedItem().toString().equals("GE")){
                curE.ensemblesgeButton.setEnabled(true);
                curE.vapDensityTextField.setEnabled(false);
                curE.jLabel1.setEnabled(false);
//                curE.chemPotMethodjButton.setEnabled(false); its now in the table...
                
                    curE.jLabel40.setEnabled(false);
                    curE.corrLengthCFTextField.setEnabled(false);
                    curE.jLabel41.setEnabled(false);
                    curE.resFreqCFTextField.setEnabled(false);
                    curE.jLabel42.setEnabled(false);
                    curE.spanFunCFTextField.setEnabled(false);
                    curE.jLabel43.setEnabled(false);
                    curE.jLabel44.setEnabled(false);
                    curE.stepcorrfuntxt.setEnabled(false);
                    curE.jLabel45.setEnabled(false);
                    curE.viewFunCFTextField.setEnabled(false);
/*                curE.jLabel8.setEnabled(false);
                curE.optPressureComboBox.setEnabled(false);
                if (curE.optPressureComboBox.getSelectedItem().toString().equals("No")){
                        curE.eData.setOptPressure("Yes");
                        curE.optPressureComboBox.setSelectedItem(curE.eData.getOptPressure());
                      //  curE.optPressureChanged();
                }*/
            } else {
                curE.ensemblesgeButton.setEnabled(false);
                curE.vapDensityTextField.setEnabled(true);
                curE.jLabel1.setEnabled(true);
//                curE.chemPotMethodjButton.setEnabled(false); its now int the potmod table
           
/*                curE.editPmFile.setEnabled(false);*/
            }
            
        }
    }
    
    private void typeOfSimulationChanged() {
        EnsemblePanel curE;
        /* MC */
        if (typeofsimulationComboBox.getSelectedItem().toString().equals("MC")){
            jLabel10.setVisible(false);
            integratorComboBox.setVisible(false);
            jLabel11.setVisible(false);
            timestepTextField.setVisible(false);
//            timestepunitComboBox.setVisible(false);
            timestepunitlabel.setVisible(false);
            jLabel22.setVisible(false);
            mcorStepsTextField.setVisible(false);
            
            jLabel12.setVisible(true);
            jLabel12.setEnabled(true);
            acceptanceTextField.setVisible(true);
            acceptanceTextField.setEnabled(true);
            optpressurelabel.setVisible(true);
            optpressurelabel.setEnabled(true);
            optPressureComboBox.setVisible(true);
            optPressureComboBox.setEnabled(true);
            commonequilabel.setVisible(true);
            commonequilabel.setEnabled(true);
            commonequiComboBox.setVisible(true);
            commonequiComboBox.setEnabled(true);
            jLabel23.setVisible(false);
            nOrientTextField.setVisible(false);
            jLabel24.setVisible(false);
            rStepsTextField.setVisible(false);
            jLabel25.setVisible(false);
            rMinRadiusTextField.setVisible(false);
            jLabel26.setVisible(false);
            rMaxRadiusTextField.setVisible(false);
            
/*ich*/     
            
            this.jLabel8.setEnabled(true);
            this.typeofensembleComboBox.setEnabled(true);
            if (this.typeofensembleComboBox.getSelectedItem().toString().equals("NPT") || this.typeofensembleComboBox.getSelectedItem().toString().equals("GE")){
                this.jLabel14.setEnabled(true);
                this.nptstepsTextField.setEnabled(true);

            } 
            else if (this.typeofensembleComboBox.getSelectedItem().toString().equals("NPH")){
                this.jLabel14.setEnabled(true);
                this.nptstepsTextField.setEnabled(true);

            } 
            else {
                this.jLabel14.setEnabled(false);
                this.nptstepsTextField.setEnabled(false);
            }
            this.jLabel13.setEnabled(true);
            this.nvtstepsTextField.setEnabled(true);
            this.jLabel15.setEnabled(true);
            this.runstepsTextField.setEnabled(true);
            this.jLabel16.setEnabled(true);
            this.resultfreqTextField.setEnabled(true);
            this.jLabel17.setEnabled(true);
            this.errorfreqTextField.setEnabled(true);
            this.jLabel18.setEnabled(true);
            this.visualfreqTextField.setEnabled(true);
            this.jLabel21.setEnabled(true);
            this.cutofftypeComboBox.setEnabled(true);
            this.jLabel28.setEnabled(true);
            this.longRangeCorrComboBox.setEnabled(true);
            longRangeCorrChanged();
            this.jLabel33.setEnabled(true);
            this.intDegFreedComboBox.setEnabled(true);
            intDegFreedChanged();
            
            java.util.List ComponentList = Arrays.asList(this.ensemblesTabbedPane.getComponents());
            Iterator CIt = ComponentList.iterator();
            
            while (CIt.hasNext()) {
                curE = (EnsemblePanel)CIt.next();
                curE.jLabel46.setEnabled(false);
            curE.transportComboBox.setEnabled(false);
                
                curE.jLabel25.setEnabled(true);
                curE.densityTextField.setEnabled(true);
                curE.densityunitComboBox.setEnabled(true);
                curE.jLabel1.setEnabled(true);
                curE.vapDensityTextField.setEnabled(true);
                curE.vapdensityunitComboBox.setEnabled(true);
                curE.jLabel28.setEnabled(true);
                curE.particlesTextField.setEnabled(true);
                
                    curE.jLabel40.setEnabled(false);
                    curE.corrLengthCFTextField.setEnabled(false);
                    curE.jLabel41.setEnabled(false);
                    curE.resFreqCFTextField.setEnabled(false);
                    curE.jLabel42.setEnabled(false);
                    curE.spanFunCFTextField.setEnabled(false);
                    curE.jLabel43.setEnabled(false);
                    curE.jLabel44.setEnabled(false);
                    curE.stepcorrfuntxt.setEnabled(false);
                    curE.jLabel45.setEnabled(false);
                    curE.viewFunCFTextField.setEnabled(false);             
                
                /*if (this.cutofftypeComboBox.getSelectedItem().toString().equals("COM")){
                  curE.jLabel24.setEnabled(true);
                  curE.cutoffTextField.setEnabled(true);
                }
                else {
                  curE.jLabel36.setEnabled(true);
                  curE.cutoffljTextField.setEnabled(true);
                  curE.jLabel37.setEnabled(true);
                  curE.cutoffddTextField.setEnabled(true);
                  curE.jLabel38.setEnabled(true);
                  curE.cutoffdqTextField.setEnabled(true);
                  curE.jLabel39.setEnabled(true);
                  curE.cutoffqqTextField.setEnabled(true);
                }
                */
                cuttoffTypeChanged();
                curE.jLabel40.setEnabled(true);
                curE.epsilonTextField.setEnabled(true);
                // curE.calcMaxCutoff.setEnabled(true);
                curE.updatePmButtons();
                
                
                if (this.typeofensembleComboBox.getSelectedItem().toString().equals("NPT")){
                    curE.jLabel23.setEnabled(true);
                    curE.pressureTextField.setEnabled(true);
                    curE.pressureunitComboBox.setEnabled(true);
                    curE.jLabel46.setEnabled(false);
                    curE.transportComboBox.setEnabled(false);
                    curE.jLabel27.setEnabled(false);
                    curE.pistonmassTextField.setEnabled(false);
                    curE.pistonmassunitComboBox.setEnabled(false);
//                    curE.jLabel8.setEnabled(true);
 //                   curE.optPressureComboBox.setEnabled(true);
                } 
                else if (this.typeofensembleComboBox.getSelectedItem().toString().equals("NPH")){
                    curE.jLabel23.setEnabled(true);
                    curE.pressureTextField.setEnabled(true);
                    curE.pressureunitComboBox.setEnabled(true);
                    curE.jLabel46.setEnabled(false);
                    curE.transportComboBox.setEnabled(false);
                    curE.jLabel27.setEnabled(false);
                    curE.pistonmassTextField.setEnabled(false);
                    curE.pistonmassunitComboBox.setEnabled(false);
/*                    curE.jLabel8.setEnabled(false);
                    curE.optPressureComboBox.setEnabled(false);
                    if (curE.optPressureComboBox.getSelectedItem().toString().equals("No")){
                        curE.eData.setOptPressure("Yes");
                        curE.optPressureComboBox.setSelectedItem(curE.eData.getOptPressure());
                      //  curE.optPressureChanged();
                    }*/
                } 
                else if (this.typeofensembleComboBox.getSelectedItem().toString().equals("GE")){
                    curE.jLabel1.setEnabled(false);
                    curE.vapDensityTextField.setEnabled(false);
                    curE.vapdensityunitComboBox.setEnabled(false);
/*                    curE.jLabel8.setEnabled(false);
                    curE.optPressureComboBox.setEnabled(false);
                    if (curE.optPressureComboBox.getSelectedItem().toString().equals("No")){
                        curE.eData.setOptPressure("Yes");
                        curE.optPressureComboBox.setSelectedItem(curE.eData.getOptPressure());
                      //  curE.optPressureChanged();
                    }*/
                } 
                else {
            /*        curE.jLabel8.setEnabled(false);
                    curE.optPressureComboBox.setEnabled(false);
                    if (curE.optPressureComboBox.getSelectedItem().toString().equals("No")){
                        curE.eData.setOptPressure("Yes");
                        curE.optPressureComboBox.setSelectedItem(curE.eData.getOptPressure());
                      //  curE.optPressureChanged();
                    }*/
                }
            }
            
            if (this.typeofensembleComboBox.getSelectedItem().toString().equals("NPH")){
                //this.gData.setTypeOfSimulation("MD");
                //typeofsimulationComboBox.setSelectedItem(this.gData.getTypeOfSimulation());
                //typeOfSimulationChanged();

            } 

            /* MD */
        } else if (typeofsimulationComboBox.getSelectedItem().toString().equals("MD")){
            
            jLabel10.setVisible(true);
            integratorComboBox.setVisible(true);
            if(qshaker.getSelectedItem().toString().equals("On")) {
                jLabel10.setEnabled(false);
                integratorComboBox.setEnabled(false);
            }
            else {
                jLabel10.setEnabled(true);
                integratorComboBox.setEnabled(true);
            }
            
            jLabel11.setVisible(true);
            timestepTextField.setVisible(true);
//            timestepunitComboBox.setVisible(true);
            timestepunitlabel.setVisible(true);
            jLabel22.setVisible(true);
            mcorStepsTextField.setVisible(true);
            
            jLabel12.setVisible(false);
            jLabel12.setEnabled(false);
            acceptanceTextField.setVisible(false);
            acceptanceTextField.setEnabled(false);
            optpressurelabel.setVisible(false);
            optpressurelabel.setEnabled(false);
            optPressureComboBox.setVisible(false);
            optPressureComboBox.setEnabled(false);
            commonequilabel.setVisible(false);
            commonequilabel.setEnabled(false);
            commonequiComboBox.setVisible(false);
            commonequiComboBox.setEnabled(false);
            if (this.gData.getMcorSteps()>0){
                this.acceptanceTextField.setEnabled(true);
                this.jLabel12.setEnabled(true);
            }
            
            jLabel23.setVisible(false);
            nOrientTextField.setVisible(false);
            jLabel24.setVisible(false);
            rStepsTextField.setVisible(false);
            jLabel25.setVisible(false);
            rMinRadiusTextField.setVisible(false);
            jLabel26.setVisible(false);
            rMaxRadiusTextField.setVisible(false);    
           
 /*ich*/      

            this.jLabel8.setEnabled(true);
            this.typeofensembleComboBox.setEnabled(true);

            this.jLabel13.setEnabled(true);
            this.nvtstepsTextField.setEnabled(true);
            if (this.typeofensembleComboBox.getSelectedItem().toString().equals("NPH")){
                this.jLabel14.setEnabled(true);
                this.nptstepsTextField.setEnabled(true);
    /* ich */                   
            }
            else if (this.typeofensembleComboBox.getSelectedItem().toString().equals("NPT")){
                this.jLabel14.setEnabled(true);
                this.nptstepsTextField.setEnabled(true); 
     /* ich */  
            }
            else {
                this.jLabel14.setEnabled(false);
                this.nptstepsTextField.setEnabled(false);
    /* ich */                                 
            }            
            
            this.jLabel15.setEnabled(true);
            this.runstepsTextField.setEnabled(true);
            this.jLabel16.setEnabled(true);
            this.resultfreqTextField.setEnabled(true);
            this.jLabel17.setEnabled(true);
            this.errorfreqTextField.setEnabled(true);
            this.jLabel18.setEnabled(true);
            this.visualfreqTextField.setEnabled(true);
            this.jLabel21.setEnabled(true);
            this.cutofftypeComboBox.setEnabled(true);
            this.jLabel28.setEnabled(true);
            this.longRangeCorrComboBox.setEnabled(true);
            longRangeCorrChanged();
            this.jLabel33.setEnabled(true);
            this.intDegFreedComboBox.setEnabled(true);
            intDegFreedChanged();
            
            EnsemblePanel ECur = (EnsemblePanel) this.ensemblesTabbedPane.getSelectedComponent();
            
            java.util.List ComponentList = Arrays.asList(this.ensemblesTabbedPane.getComponents());
            Iterator CIt = ComponentList.iterator();
            
            while (CIt.hasNext()) {
                curE = (EnsemblePanel)CIt.next();
                
                curE.jLabel46.setEnabled(true);
                curE.transportComboBox.setEnabled(true);
                
                if (this.typeofensembleComboBox.getSelectedItem().toString().equals("NPH")){
    /* ich */   curE.jLabel46.setEnabled(false);
                curE.transportComboBox.setEnabled(false);                 
            }
            else if (this.typeofensembleComboBox.getSelectedItem().toString().equals("NPT")){
     /* ich */  curE.jLabel46.setEnabled(true);
                curE.transportComboBox.setEnabled(true);
            }
            else {
    /* ich */  curE.jLabel46.setEnabled(true);
                curE.transportComboBox.setEnabled(true);                              
            } 
                curE.jLabel25.setEnabled(true);
                curE.densityTextField.setEnabled(true);
                curE.densityunitComboBox.setEnabled(true);
                curE.jLabel1.setEnabled(true);
                curE.vapDensityTextField.setEnabled(true);
                curE.vapdensityunitComboBox.setEnabled(true);
                curE.jLabel28.setEnabled(true);
                curE.particlesTextField.setEnabled(true);
                
/*                curE.jLabel8.setEnabled(false);
                curE.optPressureComboBox.setEnabled(false);
                if (curE.optPressureComboBox.getSelectedItem().toString().equals("No")){
                        curE.eData.setOptPressure("Yes");
                        curE.optPressureComboBox.setSelectedItem(curE.eData.getOptPressure());
                      //  curE.optPressureChanged();
                }*/
                
                /*
                if (this.cutofftypeComboBox.getSelectedItem().toString().equals("COM")){
                  curE.jLabel24.setEnabled(true);
                  curE.cutoffTextField.setEnabled(true);
                }
                else {
                  curE.jLabel36.setEnabled(true);
                  curE.cutoffljTextField.setEnabled(true);
                  curE.jLabel37.setEnabled(true);
                  curE.cutoffddTextField.setEnabled(true);
                  curE.jLabel38.setEnabled(true);
                  curE.cutoffdqTextField.setEnabled(true);
                  curE.jLabel39.setEnabled(true);
                  curE.cutoffqqTextField.setEnabled(true);
                }
                */
                cuttoffTypeChanged();

                if (curE.transportComboBox.getSelectedItem().toString().equals("Off")){  
 /*ich*/            curE.jLabel40.setEnabled(false);
                    curE.corrLengthCFTextField.setEnabled(false);
                    curE.jLabel41.setEnabled(false);
                    curE.resFreqCFTextField.setEnabled(false);
                    curE.jLabel42.setEnabled(false);
                    curE.spanFunCFTextField.setEnabled(false);
                    curE.jLabel43.setEnabled(false);
                    curE.jLabel44.setEnabled(false);
                    curE.stepcorrfuntxt.setEnabled(false);
                    curE.jLabel45.setEnabled(false);
                    curE.viewFunCFTextField.setEnabled(false);
                }
                else if (this.typeofensembleComboBox.getSelectedItem().toString().equals("NPH") || this.typeofensembleComboBox.getSelectedItem().toString().equals("GE")){
                    curE.jLabel40.setEnabled(false);
                    curE.corrLengthCFTextField.setEnabled(false);
                    curE.jLabel41.setEnabled(false);
                    curE.resFreqCFTextField.setEnabled(false);
                    curE.jLabel42.setEnabled(false);
                    curE.spanFunCFTextField.setEnabled(false);
                    curE.jLabel44.setEnabled(false);
                    curE.stepcorrfuntxt.setEnabled(false);
                    curE.jLabel43.setEnabled(false);
                    curE.jLabel45.setEnabled(false);
                    curE.viewFunCFTextField.setEnabled(false);
                  } 
                else{
 /*ich*/            curE.jLabel40.setEnabled(true);
                    curE.corrLengthCFTextField.setEnabled(true);
                    curE.jLabel41.setEnabled(true);
                    curE.resFreqCFTextField.setEnabled(true);
                    curE.jLabel42.setEnabled(true);
                    curE.spanFunCFTextField.setEnabled(true);
                    curE.jLabel44.setEnabled(true);
                    curE.stepcorrfuntxt.setEnabled(true);
                    curE.jLabel43.setEnabled(true);
                    curE.jLabel45.setEnabled(true);
                    curE.viewFunCFTextField.setEnabled(true);
                    }
                
                curE.jLabel40.setEnabled(true);
                curE.epsilonTextField.setEnabled(true);
                curE.calcMaxCutoff.setEnabled(true);
                curE.updatePmButtons();
                
                
                if (this.typeofensembleComboBox.getSelectedItem().toString().equals("NPT")){
                    curE.jLabel23.setEnabled(true);
                    curE.pressureTextField.setEnabled(true);
                    curE.pressureunitComboBox.setEnabled(true);
                    curE.jLabel46.setEnabled(true);
                    curE.transportComboBox.setEnabled(true);
                    curE.jLabel27.setEnabled(true);
                    curE.pistonmassTextField.setEnabled(true);
                    curE.pistonmassunitComboBox.setEnabled(true);
                }
                if (this.typeofensembleComboBox.getSelectedItem().toString().equals("NPH")){
                    curE.jLabel23.setEnabled(true);
                    curE.pressureTextField.setEnabled(true);
                    curE.pressureunitComboBox.setEnabled(true);
                    curE.jLabel46.setEnabled(true);
                    curE.transportComboBox.setEnabled(true);
                    curE.jLabel27.setEnabled(true);
                    curE.pistonmassTextField.setEnabled(true);
                    curE.pistonmassunitComboBox.setEnabled(true);
                }
                else {}
            }

            if (this.typeofensembleComboBox.getSelectedItem().toString().equals("GE")){
                this.gData.setTypeOfSimulation("MC");
                typeofsimulationComboBox.setSelectedItem(this.gData.getTypeOfSimulation());
                typeOfSimulationChanged();
            } 
            
            
            /* SVC */
        } else if (typeofsimulationComboBox.getSelectedItem().toString().equals("SVC")){
            jLabel10.setVisible(false);
            integratorComboBox.setVisible(false);
            jLabel11.setVisible(false);
            timestepTextField.setVisible(false);
//            timestepunitComboBox.setVisible(false);
            timestepunitlabel.setVisible(false);
            jLabel22.setVisible(false);
            mcorStepsTextField.setVisible(false);
            
            
            
            jLabel12.setVisible(false);
            acceptanceTextField.setVisible(false);
            optpressurelabel.setVisible(false);
            optpressurelabel.setEnabled(false);
            optPressureComboBox.setVisible(false);
            optPressureComboBox.setEnabled(false);
            commonequilabel.setVisible(false);
            commonequilabel.setEnabled(false);
            commonequiComboBox.setVisible(false);
            commonequiComboBox.setEnabled(false);
            if (this.gData.getMcorSteps()>0){
                this.acceptanceTextField.setEnabled(true);
                this.jLabel12.setEnabled(true);
            }
            
            jLabel23.setVisible(true);
            nOrientTextField.setVisible(true);
            jLabel24.setVisible(true);
            rStepsTextField.setVisible(true);
            jLabel25.setVisible(true);
            rMinRadiusTextField.setVisible(true);
            jLabel26.setVisible(true);
            rMaxRadiusTextField.setVisible(true);
            
            this.jLabel8.setEnabled(false);
            this.typeofensembleComboBox.setEnabled(false);
            this.jLabel13.setEnabled(false);
            this.nvtstepsTextField.setEnabled(false);
            this.jLabel14.setEnabled(false);
            this.nptstepsTextField.setEnabled(false);
            this.jLabel15.setEnabled(false);
            this.runstepsTextField.setEnabled(false);
            this.jLabel16.setEnabled(false);
            this.resultfreqTextField.setEnabled(false);
            this.jLabel17.setEnabled(false);
            this.errorfreqTextField.setEnabled(false);
            this.jLabel18.setEnabled(false);
            this.visualfreqTextField.setEnabled(false);
            this.jLabel21.setEnabled(false);
            this.cutofftypeComboBox.setEnabled(false);
            this.jLabel28.setEnabled(true);
            this.longRangeCorrComboBox.setEnabled(true);
            longRangeCorrChanged();
            this.jLabel33.setEnabled(true);
            this.intDegFreedComboBox.setEnabled(true);
            intDegFreedChanged();
      
            
            
            EnsemblePanel ECur = (EnsemblePanel) this.ensemblesTabbedPane.getSelectedComponent();
            
            java.util.List ComponentList = Arrays.asList(this.ensemblesTabbedPane.getComponents());
            Iterator CIt = ComponentList.iterator();
            
            while (CIt.hasNext()) {
                curE = (EnsemblePanel)CIt.next();
                curE.jLabel46.setEnabled(false);
            curE.transportComboBox.setEnabled(false);
                curE.jLabel23.setEnabled(false);
                curE.pressureTextField.setEnabled(false);
                curE.pressureunitComboBox.setEnabled(false);                
                curE.jLabel25.setEnabled(false);
                curE.densityTextField.setEnabled(false);
                curE.densityunitComboBox.setEnabled(false);
                curE.jLabel1.setEnabled(false);
                curE.vapDensityTextField.setEnabled(false);
                curE.vapdensityunitComboBox.setEnabled(false);
                //jLabel29.setEnabled(false);
                curE.jLabel27.setEnabled(false);
                curE.pistonmassTextField.setEnabled(false);
                curE.pistonmassunitComboBox.setEnabled(false);
/*                curE.jLabel8.setEnabled(false);
                curE.optPressureComboBox.setEnabled(false);
                if (curE.optPressureComboBox.getSelectedItem().toString().equals("No")){
                        curE.eData.setOptPressure("Yes");
                        curE.optPressureComboBox.setSelectedItem(curE.eData.getOptPressure());
                      //  curE.optPressureChanged();
                }*/
                curE.jLabel28.setEnabled(false);
                curE.particlesTextField.setEnabled(false);
                curE.jLabel24.setEnabled(false);
                curE.cutoffTextField.setEnabled(false);
                curE.jLabel36.setEnabled(false);
                curE.cutoffljTextField.setEnabled(false);
                curE.jLabel37.setEnabled(false);
                curE.cutoffddTextField.setEnabled(false);
                curE.jLabel38.setEnabled(false);
                curE.cutoffdqTextField.setEnabled(false);
                curE.jLabel39.setEnabled(false);
                curE.cutoffqqTextField.setEnabled(false);
                curE.jLabel40.setEnabled(false);
                curE.epsilonTextField.setEnabled(false);
                curE.calcMaxCutoff.setEnabled(false);
                
  /*ich*/       curE.jLabel40.setEnabled(false);
                    curE.corrLengthCFTextField.setEnabled(false);
                    curE.jLabel41.setEnabled(false);
                    curE.resFreqCFTextField.setEnabled(false);
                    curE.jLabel42.setEnabled(false);
                    curE.spanFunCFTextField.setEnabled(false);
                    curE.jLabel43.setEnabled(false);
                    curE.jLabel44.setEnabled(false);
                    curE.stepcorrfuntxt.setEnabled(false);
                    curE.jLabel45.setEnabled(false);
                    curE.viewFunCFTextField.setEnabled(false);              
                
                curE.updatePmButtons();
            }
        }
        
    }

    private void integratorChanged() {
        //JOptionPane.showMessageDialog(null, this.integratorComboBox.getSelectedItem().toString());
        this.gData.setIntegrator(this.integratorComboBox.getSelectedItem().toString());
    }

    
    private void cuttoffTypeChanged() {
        java.util.List ComponentList = Arrays.asList(this.ensemblesTabbedPane.getComponents());
        Iterator CIt = ComponentList.iterator();
        while (CIt.hasNext()) {
            EnsemblePanel curE = (EnsemblePanel)CIt.next();
            
            if (cutofftypeComboBox.getSelectedItem().toString().equals("Site")){
                curE.calcMaxCutoff.setEnabled(true);
                curE.jLabel24.setEnabled(true);
                curE.cutoffTextField.setEnabled(false);
                curE.jLabel36.setEnabled(true);
                curE.cutoffljTextField.setEnabled(true);
                curE.jLabel37.setEnabled(true);
                curE.cutoffddTextField.setEnabled(true);
                curE.jLabel38.setEnabled(true);
                curE.cutoffdqTextField.setEnabled(true);
                curE.jLabel39.setEnabled(true);
                curE.cutoffqqTextField.setEnabled(true);
                curE.cutoffljTextField.setText(String.valueOf(curE.eData.getCutOffLJ()));
                curE.cutoffdqTextField.setText(String.valueOf(curE.eData.getCutOffDQ()));
                curE.cutoffqqTextField.setText(String.valueOf(curE.eData.getCutOffQQ()));
                curE.cutoffddTextField.setText(String.valueOf(curE.eData.getCutOffDD()));
            } else {
                curE.calcMaxCutoff.setEnabled(true);
                curE.jLabel24.setEnabled(true);
                curE.cutoffTextField.setEnabled(true);
                curE.jLabel36.setEnabled(false);
                curE.cutoffljTextField.setEnabled(false);
                curE.jLabel37.setEnabled(false);
                curE.cutoffddTextField.setEnabled(false);
                curE.jLabel38.setEnabled(false);
                curE.cutoffdqTextField.setEnabled(false);
                curE.jLabel39.setEnabled(false);
                curE.cutoffqqTextField.setEnabled(false);
                curE.cutoffTextField.setText(String.valueOf(curE.eData.getCutOff()));
            }
        }
    }
    
    private void updatePanel() {
        
        systemofunitsComboBox.setSelectedItem(this.gData.getSystemOfUnits());
        systemOfUnitsChanged();
        lengthTextField.setText(String.valueOf(this.gData.getLengthUnit()));
        energyTextField.setText(String.valueOf(this.gData.getEnergyUnit()));
        massTextField.setText(String.valueOf(this.gData.getMassUnit()));
        typeofensembleComboBox.setSelectedItem(this.gData.getTypeOfEnsemble());
        typeOfEnsembleChanged();
        typeofsimulationComboBox.setSelectedItem(this.gData.getTypeOfSimulation());
        typeOfSimulationChanged();
        typeOfEnsembleChanged();
        integratorComboBox.setSelectedItem(this.gData.getIntegrator());
        integratorChanged();
        
        timestepTextField.setText(String.valueOf(this.gData.getTimeStep()));
//        timestepunitComboBox.setSelectedItem(this.gData.getTimeStepUnit());
        timestepunitlabel.setText(this.gData.getTimeStepUnit());
        mcorStepsTextField.setText(String.valueOf(this.gData.getMcorSteps()));
        if (this.gData.getMcorSteps()>0){
            this.acceptanceTextField.setEnabled(true);
            this.jLabel12.setEnabled(true);
        }
        acceptanceTextField.setText(String.valueOf(this.gData.getAcceptance()));
        optPressureComboBox.setSelectedItem(this.eData.getOptPressure());
        commonequiComboBox.setSelectedItem(this.eData.getcommonequi());
        nOrientTextField.setText(String.valueOf(this.gData.getNOrient()));
        rStepsTextField.setText(String.valueOf(this.gData.getRSteps()));
        rMinRadiusTextField.setText(String.valueOf(this.gData.getRMinRadius()));
        rMaxRadiusTextField.setText(String.valueOf(this.gData.getRMaxRadius()));
        
        nvtstepsTextField.setText(String.valueOf(this.gData.getNvtSteps()));
        nptstepsTextField.setText(String.valueOf(this.gData.getNptSteps()));
        mvttxtbox.setText(String.valueOf(this.gData.getmuevttSteps()));
        runstepsTextField.setText(String.valueOf(this.gData.getRunSteps()));
        resultfreqTextField.setText(String.valueOf(this.gData.getResultFreq()));
        ndbtxt.setText(String.valueOf(this.gData.getNumDinBins()));
        walltimetxt.setText(String.valueOf(this.gData.getWallTime()));
        errorfreqTextField.setText(String.valueOf(this.gData.getErrorFreq()));
        visualfreqTextField.setText(String.valueOf(this.gData.getVisualFreq()));
        cutofftypeComboBox.setSelectedItem(this.gData.getCutOffType());
        rdftxtbox.setText(String.valueOf(this.gData.getrdffreq()));
        numshelltxtbox.setText(String.valueOf(this.gData.getnumshell()));
        cuttoffTypeChanged();
        
        longRangeCorrComboBox.setSelectedItem(String.valueOf(this.gData.getLongRangeCorrType()));
        longRangeCorrChanged();
        kappaTextField.setText(String.valueOf(this.gData.getKappa()));
        nVecMaxTextField.setText(String.valueOf(this.gData.getNVecMax()));
        nsqMaxTextField.setText(String.valueOf(this.gData.getNsqMax()));
        nMaxTextField.setText(String.valueOf(this.gData.getNMax()));
        
        intDegFreedComboBox.setSelectedItem(String.valueOf(this.gData.getIntDegFreedType()));
        intDegFreedChanged();
        qshaker.setSelectedItem(String.valueOf(this.gData.getqshaker()));
        qshakerchanged();
        ljElComboBox.setSelectedItem(String.valueOf(this.gData.getLJ_El_14()));
        ljElchanged();
        
//        java.util.List ComponentList = Arrays.asList(this.ensemblesTabbedPane.getComponents());
//        Iterator CIt = ComponentList.iterator();
//            while (CIt.hasNext()) {
//                EnsemblePanel curE = (EnsemblePanel)CIt.next();
//      curE.transportComboBox.setSelectedItem(curE.eData.getTransport());
//      //JOptionPane.showMessageDialog(null, "update "+this.gData.getTransport());
//            }
    }

    /**
     *
     * @param curE
     */
    public void transportChanged(EnsemblePanel curE) {
              
                //JOptionPane.showMessageDialog(null, "in transportchanged : "+curE.transportComboBox.getSelectedItem().toString());
                if (curE.transportComboBox.getSelectedItem().toString().equals("Off"))
                {
                curE.jLabel41.setEnabled(false);
                curE.corrLengthCFTextField.setEnabled(false);
                curE.jLabel42.setEnabled(false);
                curE.resFreqCFTextField.setEnabled(false);
                curE.jLabel43.setEnabled(false);
                curE.spanFunCFTextField.setEnabled(false);
                curE.jLabel45.setEnabled(false);
                curE.jLabel44.setEnabled(false);
                curE.stepcorrfuntxt.setEnabled(false);
                curE.viewFunCFTextField.setEnabled(false);
            
        }
        else{
                curE.jLabel42.setEnabled(true);
                curE.corrLengthCFTextField.setEnabled(true);
                curE.jLabel43.setEnabled(true);
                curE.resFreqCFTextField.setEnabled(true);
                curE.jLabel44.setEnabled(true);
                curE.spanFunCFTextField.setEnabled(true);
                curE.jLabel45.setEnabled(true);
                curE.jLabel41.setEnabled(true);
                curE.stepcorrfuntxt.setEnabled(true);
                curE.viewFunCFTextField.setEnabled(true);
            
                }     
    }
//    public void transportChanged() {
//              
//            java.util.List ComponentList = Arrays.asList(this.ensemblesTabbedPane.getComponents());
//            Iterator CIt = ComponentList.iterator();
//            while (CIt.hasNext()) {
//                EnsemblePanel curE = (EnsemblePanel)CIt.next();
//                JOptionPane.showMessageDialog(null, "in transportchanged : "+curE.transportComboBox.getSelectedItem().toString());
//                if (curE.transportComboBox.getSelectedItem().toString().equals("Off"))
//                {
//                curE.jLabel41.setEnabled(false);
//                curE.corrLengthCFTextField.setEnabled(false);
//                curE.jLabel42.setEnabled(false);
//                curE.resFreqCFTextField.setEnabled(false);
//                curE.jLabel43.setEnabled(false);
//                curE.spanFunCFTextField.setEnabled(false);
//                curE.jLabel45.setEnabled(false);
//                curE.jLabel44.setEnabled(false);
//                curE.stepcorrfuntxt.setEnabled(false);
//                curE.viewFunCFTextField.setEnabled(false);
//            
//        }
//        else{
//                curE.jLabel42.setEnabled(true);
//                curE.corrLengthCFTextField.setEnabled(true);
//                curE.jLabel43.setEnabled(true);
//                curE.resFreqCFTextField.setEnabled(true);
//                curE.jLabel44.setEnabled(true);
//                curE.spanFunCFTextField.setEnabled(true);
//                curE.jLabel45.setEnabled(true);
//                curE.jLabel41.setEnabled(true);
//                curE.stepcorrfuntxt.setEnabled(true);
//                curE.viewFunCFTextField.setEnabled(true);
//            
//                }
//        }
//    }
    private void clearPanel() {
        
        initPanel();
        
    }
    
    private void addEnsemble() {
        
        String TabDesc = "";
        if (this.ensemblesTabbedPane.getComponentCount()>0){
            
            while (TabDesc.equals("")){
                TabDesc = (String)JOptionPane.showInputDialog("Description:","Ensemble " + String.valueOf(this.ensemblesTabbedPane.getComponentCount()+1));
            }
        } else {
            TabDesc = "Ensemble 1";
        }
        this.ensemblesTabbedPane.add(TabDesc, new EnsemblePanel(this));
        this.jLabel20.setText(String.valueOf(this.ensemblesTabbedPane.getComponentCount()));
        
    }
    
    private void removeEnsemble(){
        
        int ret = JOptionPane.showConfirmDialog(new JFrame(), "Are you sure, you want to remove this Ensemble ?","Delete?", 0);
        if (ret==0){
            this.ensemblesTabbedPane.remove(this.ensemblesTabbedPane.getSelectedComponent());
            this.jLabel20.setText(String.valueOf(this.ensemblesTabbedPane.getComponentCount()));
        } else {}
        
    }
    
    private void getResforEnsemble(EnsemblePanel curE, int e_num, String e_desc) {
        File pmdirDefault = null;
        pmdirDefault = this.loadDir;
        File resDefaultFile = null;
        
        /*Get Resultfilename by loaded par-File and Number of Ensemble*/
        String tempFile = this.loadFile.getName().substring(0,this.loadFile.getName().length()-4);
        tempFile = tempFile+"_"+e_num+".res";
        resDefaultFile = new File(tempFile);
        
        class pmFilter extends FileFilter {
            public String getDescription() {
                return "Only Result Files";
            }
            
            public boolean accept(File f) {
                if (f == null)
                    return false;
                if (f.isDirectory())
                    return true;
                return f.getName().toLowerCase().endsWith(".res");
            }
        }
        
        JFileChooser pmfileChooser = new JFileChooser();
        pmfileChooser.setDialogTitle("Open Resultfile for Ensemble:'"+ e_desc +"'");
        pmfileChooser.addChoosableFileFilter(new pmFilter());
        pmfileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
        pmfileChooser.setCurrentDirectory(pmdirDefault);
        pmfileChooser.setSelectedFile(resDefaultFile);
        int returnVal = pmfileChooser.showOpenDialog(this);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            
            try{
                readResFile(pmfileChooser.getSelectedFile(), curE);
                this.importRes=1;
                
            } catch (IOException e) {
                JOptionPane.showMessageDialog(new JFrame(), "An unknown error occurred");
            }
        }
    }
    
    private void readResFile(File file, EnsemblePanel curE ) throws IOException {
        FileInputStream input = new FileInputStream(file);
        BufferedReader in =    new BufferedReader(new InputStreamReader(input));
        String res = "";
        while (in.ready()) {
            String cur = in.readLine().replace("  "," ");
            cur = cur.replace("=","");
            if (cur.startsWith("#")){} else {
                res += cur + this.lineSep;
            }
        }
        in.close();
        while (res.contains("  ")){
            res = res.replace("  "," ");
        }
        res = res.substring(res.indexOf("Pressure"),res.indexOf("Cutoff")).trim();
        
        String dens [] = res.substring(res.indexOf("Density"),res.indexOf("Temperature")).split(lineSep);
        String enth [] = res.substring(res.indexOf("Enthalpy"),res.indexOf("Chemical")).split(lineSep);
        String beta [] = res.substring(res.indexOf("Isothermal"),res.indexOf("dH/dP")).split(lineSep);
        String dHdP [] = res.substring(res.indexOf("dH/dP"),res.indexOf("Isobaric")).split(lineSep);
        
        String dens_val [] = null;
        String enth_val [] = null;
        String beta_val [] = null;
        String dHdP_val [] = null;
        
        curE.eData.setGePressure0(curE.eData.getPressure());
        this.gData.setTypeOfEnsemble("GE");
        this.gData.setTypeOfSimulation("MC");
        
        if (this.gData.getSystemOfUnits().equals("SI")){
            
            dens_val = dens[1].substring(dens[1].indexOf(":")+2,dens[1].length()).split(" ");
            enth_val = enth[1].substring(enth[1].indexOf(":")+2,enth[1].length()).split(" ");
            beta_val = beta[1].substring(beta[1].indexOf(":")+2,beta[1].length()).split(" ");
            dHdP_val = dHdP[1].substring(dHdP[1].indexOf(":")+2,dHdP[1].length()).split(" ");
            
        } else { /*reduced*/
            dens_val = dens[0].substring(dens[0].indexOf(":")+2,dens[0].length()).split(" ");
            enth_val = enth[0].substring(enth[0].indexOf(":")+2,enth[0].length()).split(" ");
            beta_val = beta[0].substring(beta[0].indexOf(":")+2,beta[0].length()).split(" ");
            dHdP_val = dHdP[0].substring(dHdP[0].indexOf(":")+2,dHdP[0].length()).split(" ");
            
        }
        
        //Density
        curE.eData.setGeLiqDensity(this.validateInputDouble(dens_val[0]));
        curE.eData.setGeVarDensity(this.validateInputDouble(dens_val[1]));
        //Enthalpy
        curE.eData.setGeLiqEnthalp(this.validateInputDouble(enth_val[0]));
        curE.eData.setGeVarEnthalp(this.validateInputDouble(enth_val[1]));
        
        if (curE.eData.potModelList.size() > 0){
            for (int i=0;i<curE.eData.potModelList.size();i++){
                PotentialModel curPm = curE.eData.potModelList.get(i);
                String pmName = curPm.getFileName();
                String chemPot = res.substring(res.indexOf("Chemical potential of "+pmName),res.indexOf(lineSep,res.indexOf("Chemical potential of "+pmName)));
                String [] chemPot_val = chemPot.substring(chemPot.indexOf(":")+2,chemPot.length()).split(" ");
                //ChemPot
                curPm.setGeLiqChemPot(this.validateInputDouble(chemPot_val[0]));
                curPm.setGeVarChemPot(this.validateInputDouble(chemPot_val[1]));
                
                /*PartMolVol from Result-File*/
                if (curE.eData.potModelList.size() > 1){
                    String pmv = res.substring(res.indexOf("Partial"),res.indexOf("Isothermal"));
                    String [] part = pmv.substring(pmv.indexOf("Partial molar volume of "+pmName),pmv.indexOf(lineSep,pmv.indexOf("Partial molar volume of "+pmName)+61)).split(lineSep);
                    
                    String [] part_val = part[0].substring(part[0].indexOf(":")+2,part[0].length()).split(" ");
                    
                    curPm.setGeLiqPartMolVol(this.validateInputDouble(part_val[0]));
                    curPm.setGeVarPartMolVol(this.validateInputDouble(part_val[1]));
                    
                    /*PartMolVol 1/Density*/
                } else {
                    
                    dens_val = dens[0].substring(dens[0].indexOf(":")+2,dens[0].length()).split(" ");
                    
                    curPm.setGeLiqPartMolVol(1.0/this.validateInputDouble(dens_val[0]));
                    curPm.setGeVarPartMolVol(1.0/this.validateInputDouble(dens_val[1]));
                    
                }
                
                curPm.setGeMolarFrac(curPm.getMolarFrac());
            }
        }
        //BetaT
        curE.eData.setGeLiqBetaT(this.validateInputDouble(beta_val[0]));
        curE.eData.setGeVarBetaT(this.validateInputDouble(beta_val[1]));
        //dH/dP
        curE.eData.setGeLiqdHdP(this.validateInputDouble(dHdP_val[0]));
        curE.eData.setGeVardHdP(this.validateInputDouble(dHdP_val[1]));
        
        curE.eData.setDensity(curE.eData.getVapDensity());
        
        
        //update
        curE.updateEnsemblePanel();
        this.updatePanel();
        
        
        
    }
    
    private void readParFile(File file) throws IOException {
        try {
            FileInputStream input = new FileInputStream(file);
            BufferedReader in = new BufferedReader(new InputStreamReader(input));
            
            /* Reset Form and Data */
            int nEns = 0;
            this.gData = new General();
            this.ensemblesTabbedPane.removeAll();
            
            String i = "";
            String gInput = "";
            String[] dump;
            while (in.ready()) {
                //String cur = in.readLine().replace(" ",""); old shitty line!
                String cur = in.readLine();
                dump = cur.split("\\s+");
                if(dump.length >= 9){} else cur = cur.replace(" ","");
                
                if (cur.startsWith("#")){} else {
                    i += cur + this.lineSep;
                }
            }
            in.close();
            /* seperate gData from Ensembles */
            int ensembleIndex = i.indexOf("Temperature");
            gInput = i.substring(0, ensembleIndex-1).trim();
            
            /* seperate keys from values for each row */
            String row[] = gInput.split(this.lineSep);
            for (int j=0;j<row.length;j++){
                String[] data;
                if (row[j].contains("=")){
                    data = row[j].split("=");
                    String key = data[0].trim();
                    String val = data[1].trim();
                    readgData(key,val);
                    if (key.equals("NEnsembles"))
                        nEns = Integer.valueOf(data[1].trim()).intValue();
                }
            }
            updatePanel();
            
            /* generate Ensembles */
            String[] data = row[row.length-1].split("=");
            /*int nEns = Integer.valueOf(data[1]).intValue();*/
            
            String[] inputArray = new String[nEns];
            int ensembleStartIndex = 0;
            int ensembleEndIndex = 0;
            
            for ( int ens=0; ens<nEns ;ens++){
                ensembleStartIndex = i.indexOf("Temperature", ensembleStartIndex);
                ensembleStartIndex++;
                ensembleEndIndex = ensembleStartIndex;
                ensembleEndIndex = i.indexOf("Temperature", ensembleEndIndex);
                if (ensembleEndIndex==-1)
                    ensembleEndIndex = i.length()+1;
                inputArray[ens] = i.substring(ensembleStartIndex-1,ensembleEndIndex-1);
                ensembleEndIndex++;
                addEnsemble();
            }
            /* fill Ensemble with data */
            java.util.List ComponentList = Arrays.asList(this.ensemblesTabbedPane.getComponents());
            Iterator CIt = ComponentList.iterator();
            int e = 0;
            EnsemblePanel curE;
            while (CIt.hasNext()) {
               curE = (EnsemblePanel)CIt.next();
                //JOptionPane.showMessageDialog(null,inputArray[e]);
                readEnsemble(inputArray[e], curE);
                e++;
                /* update Panel includes Table and Buttons */
                curE.updateEnsemblePanel();
            }
            updatePanel();
            typeOfSimulationChanged();
        
//            this.timestepunitComboBox.setSelectedItem("reduced");
            //iscorrfunup = false;
            if(!isversionfound)
                versiontxt.setText("1.0");
            else
                isversionfound = false;
            
        } catch ( IOException e ) {
            
            JOptionPane.showMessageDialog(new JFrame(),"An error occured while reading "+file.getName());
        }
    }
    
    private void readgData(String key, String val){
        if(key.equals("ms2Version"))
        {
            versiontxt.setText(val);
            isversionfound = true;
        }
        else
        if (key.equals("WallTime"))
            this.gData.setWallTime(this.readInt(val));
        else
        if (key.equals("Units")){
            if (val.toUpperCase().equals("SI"))
                this.gData.setSystemOfUnits("SI");
            else if (val.toLowerCase().equals("reduced")) {
                this.gData.setTimeStepUnit("reduced");
                this.gData.setSystemOfUnits("Reduced");
            }
        } else if (key.equals("LengthUnit"))
            this.gData.setLengthUnit(this.readDouble(val));
        else if (key.equals("EnergyUnit"))
            this.gData.setEnergyUnit(this.readDouble(val));
        else if (key.equals("MassUnit"))
            this.gData.setMassUnit(this.readDouble(val));
        else if (key.equals("Simulation"))
            this.gData.setTypeOfSimulation(val.toUpperCase());
        else if(key.equals("IntDegFreed"))
            this.gData.setIntDegFreedType(val.toString());
        else if(key.equals("printIDF"))
            this.printidf.setSelectedItem(val.toString());
        else if(key.equals("Shake")){
            this.gData.setqshaker("On");
            this.tolerance.setText(val.toString());
        }
        else if(key.equals("IntraLJ_El"))
            this.intraLjElComboBox.setSelectedItem(val.toString());
        else if(key.equals("LJ_El_14"))
            this.gData.setLJ_El_14(val.toString());
        else if(key.equals("Scale"))
            this.scale.setText(val.toString());
       /* else if (key.equals("StepCorrfun"))
            this.eData.setstepFunCF(this.readInt(val));
        else if (key.equals("Corrlength"))
            this.eData.setcorrLengthCF(this.readInt(val));
        else if (key.equals("SpanCorrfun"))
            this.eData.setspanFunCF(this.readInt(val));
        else if (key.equals("ViewCorrfun"))
            this.eData.setviewFunCF(this.readInt(val));
        else if (key.equals("ResultFreqCF"))
            this.eData.setresFreqCF(this.readInt(val));*/
        else if (key.equals("Integrator")){
            if (val.equals("Gear"))
            { val = "Gear";}
            else if (val.equals("Leapfrog"))
            { val = "Leapfrog"; }
            this.gData.setIntegrator(val);
        } else if (key.equals("TimeStep"))
            this.gData.setTimeStep(this.readDouble(val));
        else if (key.equals("MCORSteps"))
            this.gData.setMcorSteps(this.readInt(val));
        else if (key.equals("Acceptance"))
            this.gData.setAcceptance(this.readDouble(val));
        else if (key.equals("NOrient"))
            this.gData.setNOrient(this.readInt(val));
        else if (key.equals("RSteps"))
            this.gData.setRSteps(this.readInt(val));
        else if (key.equals("RMinRadius"))
            this.gData.setRMinRadius(this.readDouble(val));
        else if (key.equals("RMaxRadius"))
            this.gData.setRMaxRadius(this.readDouble(val));
        else if (key.equals("Ensemble"))
            this.gData.setTypeOfEnsemble(val.toUpperCase());
        else if (key.equals("NVTSteps"))
            this.gData.setNvtSteps(this.readInt(val));
        else if (key.equals("NPTSteps"))
            this.gData.setNptSteps(this.readInt(val));
        else if (key.equals("mueVTSteps"))
            this.gData.setmuevttSteps(this.readInt(val));
        else if (key.equals("RunSteps"))
            this.gData.setRunSteps(this.readInt(val));
        else if (key.equals("ResultFreq"))
            this.gData.setResultFreq(this.readInt(val));
        else if (key.equals("ErrorsFreq"))
            this.gData.setErrorFreq(this.readInt(val));
        else if (key.equals("VisualFreq"))
            this.gData.setVisualFreq(this.readInt(val));
        else if (key.equals("RDFFreq"))
            this.gData.setrdffreq(this.readInt(val));
        else if (key.equals("NumShells"))
            this.gData.setnumshell(this.readInt(val));
        else if (key.equals("NumDinBins"))
            this.gData.setNumDinBins(this.readInt(val));
        else if (key.equals("CutoffMode")){
            if (val.toUpperCase().equals("COM"))
                val = "COM";
            if (val.toLowerCase().equals("site"))
                val = "Site";
            this.gData.setCutOffType(val);
        }
        else if (key.equals("LongRange")){
            if (val.toString().equals("Ewald"))
                val = "Ewald";
            if (val.toString().equals("ReactionField"))
                val = "ReactionField";
            if (val.toString().equals("Rodgers"))
                val = "Rodgers";
            this.gData.setLongRangeCorrType(val); 
        }
        else if (key.equals("Kappa"))
            this.gData.setKappa(this.readDouble(val));
        else if (key.equals("NsqMax"))
            this.gData.setNsqMax(this.readInt(val));
        else if (key.equals("NMax"))
            this.gData.setNMax(this.readInt(val));
        else if (key.equals("NVecMax"))
            this.gData.setNVecMax(this.readInt(val));
        /*else if (key.equals("IntDegFreed")){
            if (val.toLowerCase().equals("On"))
                val = "On";
            if (val.toLowerCase().equals("Off"))
                val = "Off";
            this.gData.setIntDegFreedType(val);
        }
        else if (key.equals("IntraLJ_El")){
            if (val.toLowerCase().equals("On"))
                val = "On";
            if (val.toLowerCase().equals("Off"))
                val = "Off";
            this.gData.setIntraLJ_Ej(val);
        } 
        else if (key.equals("LJ_El_14")){
            if (val.toLowerCase().equals("On"))
                val = "On";
            if (val.toLowerCase().equals("Off"))
                val = "Off";
            this.gData.setLJ_El_14(val);
        }*/
        else if (key.equals("NEnsembles")){}
        else if (key.equals("CorrfunMode")){
            if (val.toLowerCase().equals("yes"))
                val = "On";
            if (val.toLowerCase().equals("no"))
                val = "Off";
            this.eData.setTransport(val);
            //iscorrfunup = true;
           // JOptionPane.showMessageDialog(null, "in read : "+this.eData.getTransport());
        }
        
        else{
            /*File contains unknown parameters */
            JOptionPane.showMessageDialog(new JFrame()," Parameterfile contains unknown parameter: "+key);
        }
        
    }
    
    private void readEnsemble(String str, EnsemblePanel curE){
        String[] pmInputArray = str.split("NComponents=");
        //JOptionPane.showMessageDialog(null,pmInputArray[0]);
        String eInput = pmInputArray[0].trim();
        /* seperate keys from values for each row */
        String row[] = eInput.split(this.lineSep);
        for (int j=0;j<row.length;j++){
            String[] data;
            if (row[j].contains("=")){
                data = row[j].split("=");
                String key = data[0].trim();
                String val = data[1].trim();
                readEnsembleData(curE,key,val);
            }
        }

        int transportStart = pmInputArray[1].indexOf("CorrfunMode");
        if(transportStart <=0)
        {
            transportStart = pmInputArray[1].indexOf("StepsCorrfun");
        }
        if(transportStart <=0)
        {
            transportStart = pmInputArray[1].indexOf("Corrlength");
        }
        if(transportStart > 0)
        {
            int transportEnd = pmInputArray[1].indexOf("PotModel");
            if(transportEnd < 0)
            {
               transportEnd = pmInputArray[1].indexOf("Cutoff");
            }
            eInput = pmInputArray[1].substring(transportStart,transportEnd - 1);
            //JOptionPane.showMessageDialog(null,eInput);
            row = eInput.split(this.lineSep);
        for (int j=0;j<row.length;j++){
            String[] data;
            if (row[j].contains("=")){
                data = row[j].split("=");
                String key = data[0].trim();
                String val = data[1].trim();
                readEnsembleData(curE,key,val);
                
            }
        }
        }
        // extract the NHbond criteria 
        int hbondindexstart = pmInputArray[1].indexOf("NHBondCriteria");
        if(hbondindexstart > 0) {
            
            int hbondindexend = pmInputArray[1].indexOf("Cutoff");
            eInput = pmInputArray[1].substring(hbondindexstart,hbondindexend - 1);
            row = eInput.split(this.lineSep);
            String[] data;
            if (row[0].contains("=")){
                
                data = row[0].split("=");
                int hbondnum = Integer.parseInt(data[1].trim());
                curE.eData.numofcriteria = hbondnum;
                System.out.print(curE.eData.numofcriteria);
                curE.eData.hbondcriteria = new ArrayList<ArrayList<Object>>();
                for (int j=1;j<=hbondnum;j++) {
                    System.out.print(row[j]);
                    data = row[j].split("\\s");
                    System.out.print(Arrays.toString(data));
                    curE.eData.hbondcriteria.add(j - 1,new ArrayList<Object>());
                    for(int k = 0; k < data.length; k++) {
                        curE.eData.hbondcriteria.get(j - 1).add(data[k]);
                    }
                }
                
            }
            
        }
        int cutoffIndex = pmInputArray[1].indexOf("Cutoff");
        
        if (cutoffIndex < 1)
            cutoffIndex = pmInputArray[1].length();
        
        /* seperate PotentialModels */
        str = pmInputArray[1].substring(1,cutoffIndex-1).trim();

        readPotentialModel(str, curE);
        
        eInput = pmInputArray[1].substring(cutoffIndex-1,pmInputArray[1].length());
        row = eInput.split(this.lineSep);
        for (int j=0;j<row.length;j++){
            String[] data;
            if (row[j].contains("=")){
                data = row[j].split("=");
                String key = data[0].trim();
                String val = data[1].trim();
                readEnsembleData(curE,key,val);
            }
        }
    }
    
    private void readEnsembleData(EnsemblePanel curE, String key, String val) {
        //JOptionPane.showMessageDialog(null, key);
        //curE.transportComboBox.setSelectedItem("Off");
        if(key.equals("NHBondCriteria")){}
        else
        if (key.equals("Temperature"))
            curE.eData.setTemperature(this.readDouble(val));
        else if (key.equals("Hamiltonian"))
            curE.eData.setHamiltonian(this.readDouble(val));
        else if (key.equals("Enthalphy"))
            curE.eData.setEnthalphy(this.readDouble(val));
        else if (key.equals("Pressure"))
            curE.eData.setPressure(this.readDouble(val));
        /* GE parameter */
        else if (key.equals("Pressure0"))
            curE.eData.setGePressure0(this.readDouble(val));
        else if (key.equals("LiqDensity"))
            curE.eData.setGeLiqDensity(this.readDouble(val));
        else if (key.equals("VarDensity"))
            curE.eData.setGeVarDensity(this.readDouble(val));
        else if (key.equals("LiqEnthalpy"))
            curE.eData.setGeLiqEnthalp(this.readDouble(val));
        else if (key.equals("VarEnthalpy"))
            curE.eData.setGeVarEnthalp(this.readDouble(val));
        else if (key.equals("LiqBetaT"))
            curE.eData.setGeLiqBetaT(this.readDouble(val));
        else if (key.equals("VarBetaT"))
            curE.eData.setGeVarBetaT(this.readDouble(val));
        else if (key.equals("LiqdHdP"))
            curE.eData.setGeLiqdHdP(this.readDouble(val));
        else if (key.equals("VardHdP"))
            curE.eData.setGeVardHdP(this.readDouble(val));
        /* end of GE parameter */
        else if (key.equals("Density"))
            curE.eData.setDensity(this.readDouble(val));
        else if (key.equals("VapDensity"))
            curE.eData.setVapDensity(this.readDouble(val));
        else if (key.equals("PistonMass"))
            curE.eData.setPistonMass(this.readDouble(val));
        else if (key.equals("OptPressure")){
            if (val.toUpperCase().equals("Yes"))
                curE.eData.setOptPressure("Yes");
            else if (val.toLowerCase().equals("No"))
                curE.eData.setOptPressure("No");
        }
        else if (key.equals("NParticles"))
            curE.eData.setNPart(this.readInt(val));
        /* needed for old Par-Files -> NTest moved to PotentialModel */
        else if (key.equals("NTest"))
            curE.eData.setNTest(this.readInt(val));
        else if (key.equals("Cutoff"))
            curE.eData.setCutOff(this.readDouble(val));
        else if (key.equals("CutoffLJ"))
            curE.eData.setCutOffLJ(this.readDouble(val));
        else if (key.equals("CutoffDQ"))
            curE.eData.setCutOffDQ(this.readDouble(val));
        else if (key.equals("CutoffDD"))
            curE.eData.setCutOffDD(this.readDouble(val));
        else if (key.equals("CutoffQQ"))
            curE.eData.setCutOffQQ(this.readDouble(val));
        else if (key.equals("Epsilon"))
            curE.eData.setEps(this.readDouble(val));
        else if (key.equals("CorrfunMode")){
            if (val.toLowerCase().equals("yes"))
                val = "On";
            if (val.toLowerCase().equals("no"))
                val = "Off";
            
            this.eData.setTransport(val);
           // JOptionPane.showMessageDialog(null, "read transport property : "+this.eData.getTransport());
        }
        else if (key.equals("StepsCorrfun"))
        curE.eData.setstepFunCF(this.readInt(val));
        else if (key.equals("Corrlength"))
            curE.eData.setcorrLengthCF(this.readInt(val));
        else if (key.equals("SpanCorrfun"))
            curE.eData.setspanFunCF(this.readInt(val));
        else if (key.equals("ViewCorrfun"))
            curE.eData.setviewFunCF(this.readInt(val));
        else if (key.equals("ResultFreqCF"))
            curE.eData.setresFreqCF(this.readInt(val));
/*        else if (key.equals("FluctFreq"))
            curE.eData.setFluctFreq(this.readInt(val));
        else if (key.equals("NFullFluct"))
            curE.eData.setNFullFluct(this.readInt(val));
        else if (key.equals("MaxCounter"))
            curE.eData.setMaxCounter(this.readInt(val));*/
        else if (key.equals("")){
            
        } else{
            /* File contains unknown parameters */
            JOptionPane.showMessageDialog(new JFrame()," Parameterfile contains unknown parameter:"+key);
        }
    }
    
    private void readPotentialModel(String str, EnsemblePanel curE){
        
        String[] pmString = str.split("PotModel=");
        
        for (int j=1;j<pmString.length;j++){
            String [] row = pmString[j].split(this.lineSep);
            
            String pm = row[0];
            Double molfrac = 0.0;
            String perm = "No";
            Double therpartentahlpy = 0.0;
            String[] data = null;
            if (!this.gData.getTypeOfSimulation().equals("SVC")){
                data = row[1].split("=");
                if(data[0].trim().equals("MoleFract"))
                    molfrac = this.readDouble(data[1]);
                
                data = row[2].split("=");
                if(data[0].trim().equals("PartMolEnt"))
                   therpartentahlpy = Double.parseDouble(data[1].trim());
                
                data = row[3].split("=");
                if(data[0].trim().equals("Permeability"))
                    perm = data[1].trim();
            }
            

            String fn = this.loadFile.getName();
            String path = this.loadFile.getAbsolutePath();
            path = path.substring(0,(path.length()-fn.length()));
            File pmfile = new File(path+pm);
            
            PotentialModel curPm = new PotentialModel(pmfile,molfrac,perm);
            curPm.setpartmolarenthalpy(therpartentahlpy);
            curE.eData.potModelList.add(curPm);
            
            /* Copy NTest from old par-Files */
            if (curE.eData.getNTest()!= 0) {
                curPm.setChemPotMethod("Widom");
                curPm.setNTest(curE.eData.getNTest());
            }
            /* Create EtaXi List */
            if (curE.eData.potModelList.size() > 1){
                curE.addEtaXi(curPm);
                
            }
            /* Read ChemPotMethod */
            if (pmString[j].contains("ChemPotMethod")){
              
                String cpm = pmString[j];
                cpm =cpm.substring(cpm.indexOf("ChemPotMethod"),cpm.length());
                String [] cpmrow = cpm.split(this.lineSep);
                for (int k=0;k<cpmrow.length;k++){
                    if (cpmrow[k].contains("=")){
                        data = cpmrow[k].split("=");
                        String key = data[0];
                        String val = data[1];

                        readChemPotMethod(curPm,key,val,curE);
                    }
                }
                if (curPm.getChemPotMethod().equals("GradIns") && !curPm.getWeightFactorsType().toLowerCase().equals("auto")){
                    String wf = pmString[j];
                    int end = 0;
                    if (wf.toLowerCase().contains("eta"))
                        end = wf.toLowerCase().indexOf("eta");
                    else
                        end = wf.length();
                    wf = wf.substring(wf.indexOf("WeightFactors"),end);
                    curPm.setWeightFactors(wf.substring(wf.indexOf(lineSep)).trim());
                }
            }
            
            /* Read GE-Parameter if necessary */
            if (this.typeofensembleComboBox.getSelectedItem().equals("GE")){
                for (int k=1;k<row.length;k++){
                    if (row[k].contains("=")){
                        data = row[k].split("=");
                        String key = data[0];
                        String val = data[1];
                        readPotentialModelGE(curPm,key,val);
                    }
                }
            }
        }
        
        /* Read Eta and Xi from File */
        if (curE.eData.getEtaxi() != null){
            int etaXiIndex = str.toLowerCase().indexOf("eta");
            str = str.substring(etaXiIndex,str.length()).trim();
            String [] row = str.split(this.lineSep);
            
            EtaXi curEtaXi = curE.eData.getEtaxi();
            int k=0;
            while (curEtaXi != null){
                String[] data = null;
                /* Read Eta */
                if (row[k].contains("=")){
                data = row[k].split("=");
                curEtaXi.setEta(this.readDouble(data[1]));
                k++;
                }
                /* Read Xi */
                if (row[k].contains("=")){
                data = row[k].split("=");
                curEtaXi.setXi(this.readDouble(data[1]));
                k++;
                curEtaXi = curEtaXi.getNext();
                }
                else k++;
            }
        }
    }
    
    private void readPotentialModelGE(PotentialModel curPm, String key, String val){
        if (key.equals("MoleFract")){
            
        } else if (key.equals("LiqMoleFract"))
            curPm.setGeMolarFrac(this.readDouble(val));
        else if (key.equals("ChemPot"))
            curPm.setGeLiqChemPot(this.readDouble(val));
        else if (key.equals("VarChemPot"))
            curPm.setGeVarChemPot(this.readDouble(val));
        else if (key.equals("PartMolVol"))
            curPm.setGeLiqPartMolVol(this.readDouble(val));
        else if (key.equals("VarPartMolVol"))
            curPm.setGeVarPartMolVol(this.readDouble(val));
         else if (key.equals("eta") || key.equals("xi"))
                    {  }
        else if (key.equals("")){
            
        } else{
            /* File contains unknown parameters */
            JOptionPane.showMessageDialog(new JFrame()," Parameterfile contains unknown GE-Parameter: "+key);
        }
    }
    
    private void readChemPotMethod(PotentialModel curPm, String key, String val, EnsemblePanel curE){
        if (key.equals("MoleFract")){
            
        }
        else if (key.equals("ChemPotMethod")){
            if (val.toLowerCase().equals("gradins"))
                val = "GradIns";
            else if (val.toLowerCase().equals("widom"))
                val = "Widom";
            else
                val = val.trim();
             
            curPm.setChemPotMethod(val);
            
        }
        else if (key.equals("LambdaMin"))
            curPm.setlambdamin(Double.parseDouble(val));
        else if (key.equals("LambdaMax"))
            curPm.setlambdamax(Double.parseDouble(val));
        else if (key.equals("NBins"))
            curPm.setnbins(Integer.parseInt(val));
        else if (key.equals("LambdaStepMax"))
            curPm.setlabdastepmax(Double.parseDouble(val));
        else if (key.equals("LambdaExponent"))
            curPm.setlambdaexp(Integer.parseInt(val));
        else if (key.equals("NTest"))
            curPm.setNTest(this.readInt(val));
        else if (key.equals("WeightFactors")){
            if (val.toLowerCase().equals("optset"))
                val ="OptSet";
            else if (val.toLowerCase().equals("guess"))
                val = "Guess";
            else val = "Auto";
            curPm.setWeightFactorsType(val);
        } else if (key.equals("")){
        } else if (key.equals("eta")){
        } else if (key.equals("xi")){
        } else{
            /* File contains unknown parameters */
            JOptionPane.showMessageDialog(new JFrame()," Parameterfile contains unknown ChemPotMethod-Parameter: "+key);
        }
       
    }
    
    
    
    
    private double readDouble(String str){
        Double ret = 0.0;
        try {
            ret = Double.parseDouble(str);
        } catch(Exception e) {
            JOptionPane.showMessageDialog(new JFrame(),"Error: "+str+" is not a decimal number!");
        }
        //ret = Double.valueOf(str).doubleValue();
        return ret;
    }
    
    private int readInt(String str){
        int ret = 0;
        try {
            ret = Integer.parseInt(str);
        } catch(Exception e) {
            JOptionPane.showMessageDialog(new JFrame(),"Error: "+str+"is not a number!");
        }
        //ret = Integer.valueOf(str).intValue();
        return ret;
    }
    
    
   private void writeParFile(File file) throws IOException {
        
        try {
            PrintWriter o = new PrintWriter(
                    new BufferedWriter(
                    new FileWriter(file) ) );
            
            //java.text.DecimalFormat df = new java.text.DecimalFormat("0.0#####");
            
            o.println("#"+this.lineSep+"# MS2 Parameterfile created with ms2par"+this.lineSep+"#"+this.lineSep);
            o.println("#"+this.lineSep+"# General simulation parameters"+this.lineSep+"#");
            o.println("ms2Version     =       " + this.versiontxt.getText());
            o.println("WallTime     =       " + this.gData.getWallTime());
            o.println("Units       =       " + this.gData.getSystemOfUnits());
            if (!this.systemofunitsComboBox.getSelectedItem().toString().equals("Reduced")){
               o.println("LengthUnit  =       " + this.gData.getLengthUnit());
               o.println("EnergyUnit  =       " + this.gData.getEnergyUnit());
               o.println("MassUnit    =       " + this.gData.getMassUnit());}
            else{}
            o.println("Simulation  =       " + this.typeofsimulationComboBox.getSelectedItem().toString());
            o.println("IntDegFreed =       "+this.intDegFreedComboBox.getSelectedItem().toString());
            if(this.intDegFreedComboBox.getSelectedItem().toString().equals("On"))
            {
                o.println("printIDF    =       "+this.printidf.getSelectedItem().toString());
               // o.println("Shake       =       "+this.qshaker.getSelectedItem().toString());
                if(this.qshaker.getSelectedItem().toString().equals("On"))
                        o.println("Shake       =       "+this.tolerance.getText());
                   // o.println("Tolerance   =       "+this.tolerance.getText());
                o.println("IntraLJ_El  =       "+this.intraLjElComboBox.getSelectedItem().toString());
                o.println("LJ_El_14    =       "+this.ljElComboBox.getSelectedItem().toString());
                if(this.ljElComboBox.getSelectedItem().toString().equals("On"))
                    o.println("Scale       =       "+this.scale.getText());
            }
            if (this.typeofsimulationComboBox.getSelectedItem().equals("MD")){
                o.println("Integrator  =       " + this.gData.getIntegrator());
 //               if (this.timestepunitComboBox.getSelectedItem().toString().equals("femtosec"))
                if(this.systemofunitsComboBox.getSelectedItem().toString().equals("Reduced"))
                    o.println("TimeStep    =       " + this.gData.getReducedTimeStep(this.gData.getTimeStep()));
                else
                    o.println("TimeStep    =       " + this.gData.getTimeStep());
                
                if (this.gData.getMcorSteps()>0)
                    o.println("Acceptance    =       " + this.gData.getAcceptance());
            } else if (this.typeofsimulationComboBox.getSelectedItem().equals("MC"))
                o.println("Acceptance    =     " + this.gData.getAcceptance());
            else if (this.typeofsimulationComboBox.getSelectedItem().equals("SVC")){
                o.println("NOrient     =       " +this.gData.getNOrient());
                o.println("RSteps      =       " +this.gData.getRSteps());
                o.println("RMinRadius  =       " +this.gData.getRMinRadius());
                o.println("RMaxRadius  =       " +this.gData.getRMaxRadius());
                
            } else{}
            
            if (!this.typeofsimulationComboBox.getSelectedItem().equals("SVC")){
                
                o.println("Ensemble    =       " + this.typeofensembleComboBox.getSelectedItem().toString());
                if (this.typeofsimulationComboBox.getSelectedItem().equals("MD")){
                    o.println("MCORSteps   =       " + this.gData.getMcorSteps());
                }
                o.println("NVTSteps    =       " + this.gData.getNvtSteps());
                if (this.typeofensembleComboBox.getSelectedItem().equals("NPT"))
                    o.println("NPTSteps    =       " + this.gData.getNptSteps());
                if (this.typeofensembleComboBox.getSelectedItem().equals("NPH"))
                    o.println("NPTSteps    =       " + this.gData.getNptSteps());
                if (this.typeofensembleComboBox.getSelectedItem().equals("GE"))
                    o.println("mueVTSteps  =       " + this.gData.getNptSteps());
                o.println("RunSteps    =       " + this.gData.getRunSteps());
                o.println("ResultFreq  =       " + this.gData.getResultFreq());
                o.println("ErrorsFreq  =       " + this.gData.getErrorFreq());
                o.println("VisualFreq  =       " + this.gData.getVisualFreq());
                o.println("RDFFreq     =       " + this.gData.getrdffreq());
                o.println("CutoffMode  =       " + this.cutofftypeComboBox.getSelectedItem().toString());
            }
                o.println("NumDinBins  =       "+this.gData.getNumDinBins());
                // 
                o.println("LongRange   =       "+ this.longRangeCorrComboBox.getSelectedItem().toString());
                o.println("Kappa       =       "+ this.kappaTextField.getText());
                o.println("NsqMax      =       "+ this.nsqMaxTextField.getText());
                o.println("NMax        =       "+ this.nMaxTextField.getText());
                o.println("NVecMax     =       "+ this.nVecMaxTextField.getText());
                //
            o.println("NEnsembles  =       "+ this.ensemblesTabbedPane.getComponentCount());
            
           
            /* For each Ensemble */
            java.util.List ComponentList = Arrays.asList(this.ensemblesTabbedPane.getComponents());
            Iterator CIt = ComponentList.iterator();
            int e=1;
            while (CIt.hasNext()) {
                boolean cpmGI = false;
                EnsemblePanel curE = (EnsemblePanel)CIt.next();
                
                o.println(this.lineSep+"#"+this.lineSep+"# Ensemble "+ e + this.lineSep+"#");
                if (this.typeofensembleComboBox.getSelectedItem().equals("GE")){
                    o.println("\n#Fl�ssigkeit");
                }
                if (curE.temperatureunitComboBox.getSelectedItem().toString().equals("°C"))
                    o.println("Temperature = "+ curE.eData.CelsiusToKelvin(curE.eData.getTemperature()));
                else
                    o.println("Temperature = "+ curE.eData.getTemperature());
                
                if (!this.typeofsimulationComboBox.getSelectedItem().equals("SVC")){
                    
                    if (this.typeofensembleComboBox.getSelectedItem().toString().equals("NPT")) {
                        if (curE.pressureunitComboBox.getSelectedItem().toString().equals("bar"))
                            o.println("Pressure    = "+ curE.eData.barToMPa(curE.eData.getPressure()));
                        else
                            o.println("Pressure    = "+ curE.eData.getPressure());
                    }
                    if (this.typeofensembleComboBox.getSelectedItem().toString().equals("NPH")) {
                        o.println("Enthalphy    = "+ curE.eData.getEnthalphy());
                        if (curE.pressureunitComboBox.getSelectedItem().toString().equals("bar")) {
                            o.println("Pressure    = "+ curE.eData.barToMPa(curE.eData.getPressure()));
                        }
                        else {
                            o.println("Pressure    = "+ curE.eData.getPressure());
                        }
                            
                    }
                    if (this.typeofensembleComboBox.getSelectedItem().toString().equals("NVE")) {
                        o.println("Hamiltonian    = "+ curE.eData.getHamiltonian());
                    }
                    /* if GE write GE parameter */
                    if (this.typeofensembleComboBox.getSelectedItem().equals("GE")){
                        o.println("Pressure0   = "+curE.eData.getGePressure0());
                        o.println("LiqDensity  = "+curE.eData.getGeLiqDensity());
                        o.println("VarDensity  = "+curE.eData.getGeVarDensity());
                        o.println("LiqEnthalpy = "+curE.eData.getGeLiqEnthalp());
                        o.println("VarEnthalpy = "+curE.eData.getGeVarEnthalp());
                        o.println("LiqBetaT    = "+curE.eData.getGeLiqBetaT());
                        o.println("VarBetaT    = "+curE.eData.getGeVarBetaT());
                        o.println("LiqdHdP     = "+curE.eData.getGeLiqdHdP());
                        o.println("VardHdP     = "+curE.eData.getGeVardHdP());
                        
                        o.println("\n#Dampf");
                    }
                    
                    
                    
                    o.println("Density     = "+ curE.eData.getDensity());
                    
                    if (!this.typeofensembleComboBox.getSelectedItem().equals("GE")){
                        if (curE.eData.getVapDensity()!=0.0)
                            o.println("VapDensity  = "+curE.eData.getVapDensity());
                    }
                    
                    if (this.typeofensembleComboBox.getSelectedItem().toString().equals("NPT") && this.typeofsimulationComboBox.getSelectedItem().toString().equals("MD"))
                        o.println("PistonMass  = "+ curE.eData.getPistonMass());
                    if (this.typeofensembleComboBox.getSelectedItem().toString().equals("NPH") && this.typeofsimulationComboBox.getSelectedItem().toString().equals("MD"))
                        o.println("PistonMass  = "+ curE.eData.getPistonMass());

                    o.println("OptPressure = "+ curE.eData.getOptPressure());
                    
                    o.println("NParticles  = "+ curE.eData.getNPart());

/*                    if (transportComboBox.getSelectedItem().toString().equals("On")){
                        o.println("Corrlength  = "+ curE.eData.getcorrLengthCF());
                        o.println("SpanCorrfun = "+ curE.eData.getspanFunCF());
                        o.println("ViewCorrfun = "+ curE.eData.getviewFunCF());
                        o.println("ResultFreqCF= "+ curE.eData.getresFreqCF());
                    }*/

                /* Old Style NTest do not edit
                if (this.typeofensembleComboBox.getSelectedItem().toString().equals("GE")){} else
                    o.println("NTest       = "+ curE.eData.getNTest());
                 */
                    
                    /* For Each Potential Model */
                }
                //TODO PMLIST check
                
                if (curE.eData.potModelList.size() > 0){
                    
                    o.println("NComponents = "+curE.eData.potModelList.size());   

                    Iterator pmIt = curE.eData.potModelList.iterator();
                    while (pmIt.hasNext()) {
                        PotentialModel pm = (PotentialModel)pmIt.next();
                        
                        o.println("\nPotModel    = "+ pm.getFileName());
                        
                        if (!this.typeofsimulationComboBox.getSelectedItem().equals("SVC")){
                            o.println("MoleFract  = "+ pm.getMolarFrac());
                            o.println("PartMolEnt = "+ pm.getpartmolarenthalpy());
                            o.println("Permeability  = "+ pm.getpermeablity());
                            if (!this.typeofensembleComboBox.getSelectedItem().equals("GE")){
                                o.println("ChemPotMethod= "+ pm.getChemPotMethod());
                                
                                if (pm.getChemPotMethod().equals("Widom")){
                                    
                                    o.println("NTest        = "+ pm.getNTest());
                                } else if (pm.getChemPotMethod().equals("GradIns")){
                                    cpmGI = true;
                                    o.println("WeightFactors="+ pm.getWeightFactorsType());
                                    if (!pm.getWeightFactorsType().toLowerCase().equals("auto"))
                                        o.println(pm.getWeightFactors());
                                }else if(pm.getChemPotMethod().equals("ThermoInt")){
                                        o.println("LambdaMin    = "+ pm.getlambdamin());
                                        o.println("LambdaMax    = "+ pm.getlambdamax());
                                        o.println("NBins        = "+ pm.getnbins());
                                        o.println("LambdaStepMax= "+ pm.getlambdastepmax());
                                        o.println("LambdaExponent= "+ pm.getlambdaexp());
                                }
                            }
                            
                            /* if GE write GE parameter */
                            if (this.typeofensembleComboBox.getSelectedItem().equals("GE")){
                                o.println("LiqMoleFract = "+pm.getGeMolarFrac());
                                o.println("ChemPot       = "+pm.getGeLiqChemPot());
                                o.println("VarChemPot    = "+pm.getGeVarChemPot());
                                
                                //if (curE.eData.potModelList.size() > 1){
                                o.println("PartMolVol    = "+pm.getGeLiqPartMolVol());
                                o.println("VarPartMolVol = "+pm.getGeVarPartMolVol());
                                //}
                            }
                        }
                    }
                } else {
                    o.println("NComponents = 0");
                }
                
                if (curE.transportComboBox.getSelectedItem().toString().equals("On")){
                        o.println("CorrfunMode  = yes");
                        o.println("StepsCorrfun = "+ curE.eData.getstepFunCF());
                        o.println("Corrlength  = "+ curE.eData.getcorrLengthCF());
                        o.println("SpanCorrfun = "+ curE.eData.getspanFunCF());
                        o.println("ViewCorrfun = "+ curE.eData.getviewFunCF());
                        o.println("ResultFreqCF= "+ curE.eData.getresFreqCF());
                    } 

                /* Eta & Xi Parameter */
                EtaXi curEtaXi = curE.eData.getEtaxi();
                int count = 0;
                o.println("");
                if (curEtaXi !=null){
                    o.println("");
                    while (curEtaXi != null){
                        o.println("#"+curE.eData.Commentsdata.get(count).toString());
                        o.println("eta         = "+curEtaXi.getEta());
                        o.println("xi          = "+curEtaXi.getXi());
                         o.println("");
                        curEtaXi = curEtaXi.getNext();
                        count++;
                    }
                }
                o.println("");
                if(curE.eData.hbondcriteria != null)
                {
                    o.println("NHBondCriteria = "+curE.eData.numofcriteria);
                    for(int i=0; i<curE.eData.hbondcriteria.size();i++)
                    {
                        for(int j=0;j<9;j++)
                        {
                            o.print(curE.eData.hbondcriteria.get(i).get(j).toString()+" ");
                        }
                        o.print("\n");
                    }
                }
                o.println("");
                if (!this.typeofsimulationComboBox.getSelectedItem().equals("SVC")){
                    /* Cutoff */
                    if (this.cutofftypeComboBox.getSelectedItem().toString().equals("COM"))
                        o.println("Cutoff      = "+ curE.eData.getCutOff());
                    else if (this.cutofftypeComboBox.getSelectedItem().toString().equals("Site")){
                        o.println("CutoffLJ    = "+ curE.eData.getCutOffLJ());
                        o.println("CutoffDD    = "+ curE.eData.getCutOffDD());
                        o.println("CutoffDQ    = "+ curE.eData.getCutOffDQ());
                        o.println("CutoffQQ    = "+ curE.eData.getCutOffQQ());
                    } else {}
                    
                    o.println("Epsilon     = "+ curE.eData.getEps());
                    
                }
                o.flush();
                /* next Ensemble */
                e++;
            }
            
            o.flush();
            o.close();
        } catch (Exception e) {
            System.out.println("error123 " + e);
        }
        
    }
    
    private void copyPmFiles(String path, File saveDir) {
        this.allPmFileOption = false;
        this.pmFileOverwrite = false;
        /* For each Ensemble */
        java.util.List ComponentList = Arrays.asList(this.ensemblesTabbedPane.getComponents());
        Iterator CIt = ComponentList.iterator();
        while (CIt.hasNext()) {
            EnsemblePanel curE = (EnsemblePanel)CIt.next();
            
            /* For Each Potential Model */
            
            if (curE.eData.potModelList.size() > 0){
                
                String[] children = saveDir.list();
                if (children == null) {
                    // Either dir does not exist or is not a directory
                } else {
                    FilenameFilter pmFilter = new FilenameFilter() {
                        public boolean accept(File dir, String name) {
                            return name.endsWith(".pm");
                        }
                    };
                    children = saveDir.list(pmFilter);
                }
                
                Iterator pmIt = curE.eData.potModelList.iterator();
                
                while (pmIt.hasNext()) {
                    PotentialModel curPm = (PotentialModel)pmIt.next();
                    
                    String pmFileName = curPm.getFileName();
                    boolean pmexists = false;
                    for (int i=0; i<children.length; i++) {
                        String filename = children[i];
                        if(filename.equals(pmFileName)){
                            pmexists = true;
                        }
                    }
                    if (pmexists){
                        if (!allPmFileOption){
                            pmFileOverwrite = false;
                            
                            PmCopyDialog dialog = new PmCopyDialog(this,true,curPm);
                            dialog.setVisible(true);
                        }
                        if (pmFileOverwrite) {
                            copyPmFileData(curPm, path+curPm.getFileName());
                        }
                        
                    } else {
                        String pmCopyPath = path+pmFileName;
                        copyPmFileData(curPm, pmCopyPath);
                        
                    }
                }
                /* necessary if filename changed */
                curE.updatePmTable();
            }
        }
    }
    
    
    private void copyPmFileData(PotentialModel curPm, String pmCopyPath){
        File file = new File(pmCopyPath);
        
        try {
            PrintWriter pmo = new PrintWriter(
                    new BufferedWriter(
                    new FileWriter(file) ) );
            pmo.print(curPm.getPmFileData());
            pmo.close();
        } catch(Exception e){
            JOptionPane.showMessageDialog(new JFrame(),pmCopyPath+" could not be copied!");
        }
        
    }
    
    
    private void checkMolFrac() {
        double eps = 1.E-3;
        /* For each Ensemble */
        java.util.List ComponentList = Arrays.asList(this.ensemblesTabbedPane.getComponents());
        Iterator CIt = ComponentList.iterator();
        while (CIt.hasNext()) {
            EnsemblePanel curE = (EnsemblePanel)CIt.next();
            Iterator pmIt = curE.eData.potModelList.iterator();
            double molfracsum = 0.0;
            while (pmIt.hasNext()) {
                PotentialModel pm = (PotentialModel)pmIt.next();
                molfracsum += pm.getMolarFrac();
            }
            if (molfracsum <= 1.0-eps || molfracsum >= 1.0+eps ){
                int ret = JOptionPane.showConfirmDialog(new JFrame(), "The sum of the mole fractions is not \"1\"!\nWould you like to do an automatic correction?","",0);
                if (ret==0) {
                    for (int i=0;i<curE.eData.potModelList.size();i++ ){
                        PotentialModel pm = curE.eData.potModelList.get(i);
                        pm.setMolarFrac(pm.getMolarFrac()/molfracsum);
                    }
                }
            }
        }
    }
    private void calcNVecMax(){

        if (this.gData.getNsqMax()==0)
            JOptionPane.showMessageDialog(new JFrame()," NsqMax = 0 ");
        else if (this.gData.getNMax()==0)
            JOptionPane.showMessageDialog(new JFrame()," NMax = 0 ");
        else {
            int i,j,k;
            int counter;
            double sum;
            sum = 0.0;
            counter = 0;
            for (i=0;i<=this.gData.getNMax();i++)
               {for (j=-this.gData.getNMax();j<=this.gData.getNMax();j++)
                  {for (k=-this.gData.getNMax();k<=this.gData.getNMax();k++)
                   {sum = i*i + j*j + k*k;
                    if ((sum <= this.gData.getNsqMax()) && (sum != 0.0))
                         {counter = counter + 1;}
                   }
                  }
               }
            this.gData.setNVecMax(counter);
            this.nVecMaxTextField.setText(String.valueOf(this.gData.getNVecMax()));
        }

    }
    
    private void showConsoleMessage() {
        System.out.println("********************************************");
        System.out.println("*             Version 2.1.1                *\n" +
                           "*                Authors                   *\n" +
                           "*   Arnaud Diffo Kaze / Stephan Deublein   *");
        System.out.println("*     This tool is under development,      *\n" +
                           "*      please report bugs and errors       *\n" +
                           "*    with short description to Stephan     *");
        System.out.println("********************************************");
    }
    
    /**
     *
     */
    public boolean allPmFileOption = false;

    /**
     *
     */
    public boolean pmFileOverwrite = false;

    /**
     *
     */
    public boolean isversionfound=false;

    /**
     *
     */
    public General gData;

    /**
     *
     */
    public Ensemble eData;
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JDialog EnsembleAddDialog;
    private javax.swing.JTextField acceptanceTextField;
    private javax.swing.JButton clearButton;
    public javax.swing.JComboBox commonequiComboBox;
    public javax.swing.JLabel commonequilabel;
    public javax.swing.JComboBox cutofftypeComboBox;
    private javax.swing.JTextField energyTextField;
    public static javax.swing.JTabbedPane ensemblesTabbedPane;
    private javax.swing.JButton ensemblesaddButton;
    private javax.swing.JButton ensemblesremoveButton;
    private javax.swing.JTextField errorfreqTextField;
    private javax.swing.JLabel filenameLabel;
    private javax.swing.JPanel generalPanel;
    private javax.swing.JButton helpButton;
    private javax.swing.JComboBox intDegFreedComboBox;
    private javax.swing.JComboBox integratorComboBox;
    private javax.swing.JComboBox intraLjElComboBox;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel10;
    private javax.swing.JLabel jLabel11;
    private javax.swing.JLabel jLabel12;
    private javax.swing.JLabel jLabel13;
    private javax.swing.JLabel jLabel14;
    private javax.swing.JLabel jLabel15;
    private javax.swing.JLabel jLabel16;
    private javax.swing.JLabel jLabel17;
    private javax.swing.JLabel jLabel18;
    private javax.swing.JLabel jLabel19;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel20;
    private javax.swing.JLabel jLabel21;
    private javax.swing.JLabel jLabel22;
    private javax.swing.JLabel jLabel23;
    private javax.swing.JLabel jLabel24;
    private javax.swing.JLabel jLabel25;
    private javax.swing.JLabel jLabel26;
    private javax.swing.JLabel jLabel27;
    private javax.swing.JLabel jLabel28;
    private javax.swing.JLabel jLabel29;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel30;
    private javax.swing.JLabel jLabel31;
    private javax.swing.JLabel jLabel32;
    private javax.swing.JLabel jLabel33;
    private javax.swing.JLabel jLabel34;
    private javax.swing.JLabel jLabel35;
    private javax.swing.JLabel jLabel36;
    private javax.swing.JLabel jLabel37;
    private javax.swing.JLabel jLabel38;
    private javax.swing.JLabel jLabel39;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel40;
    private javax.swing.JLabel jLabel41;
    private javax.swing.JLabel jLabel42;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel6;
    private javax.swing.JLabel jLabel7;
    private javax.swing.JLabel jLabel8;
    private javax.swing.JLabel jLabel9;
    private javax.swing.JTextField kappaTextField;
    private javax.swing.JScrollPane leftjscroll;
    private javax.swing.JTextField lengthTextField;
    private javax.swing.JComboBox ljElComboBox;
    private javax.swing.JButton loadButton;
    private javax.swing.JComboBox longRangeCorrComboBox;
    private javax.swing.JPanel mainjpanel;
    private javax.swing.JScrollPane mainjscroll;
    private javax.swing.JTextField massTextField;
    private javax.swing.JTextField mcorStepsTextField;
    private javax.swing.JLabel mvtlabel;
    private javax.swing.JTextField mvttxtbox;
    private javax.swing.JTextField nMaxTextField;
    private javax.swing.JTextField nOrientTextField;
    private javax.swing.JTextField nVecMaxTextField;
    private javax.swing.JTextField ndbtxt;
    private javax.swing.JTextField nptstepsTextField;
    private javax.swing.JTextField nsqMaxTextField;
    private javax.swing.JLabel numshelllabel;
    private javax.swing.JTextField numshelltxtbox;
    private javax.swing.JTextField nvtstepsTextField;
    public javax.swing.JComboBox optPressureComboBox;
    public javax.swing.JLabel optpressurelabel;
    private javax.swing.JComboBox printidf;
    private javax.swing.JComboBox qshaker;
    private javax.swing.JTextField rMaxRadiusTextField;
    private javax.swing.JTextField rMinRadiusTextField;
    private javax.swing.JTextField rStepsTextField;
    private javax.swing.JLabel rdflabel;
    private javax.swing.JTextField rdftxtbox;
    private javax.swing.JTextField resultfreqTextField;
    private javax.swing.JScrollPane rightjscroll;
    private javax.swing.JTextField runstepsTextField;
    private javax.swing.JButton saveButton;
    private javax.swing.JTextField scale;
    private javax.swing.JPanel simPanel;
    public javax.swing.JComboBox systemofunitsComboBox;
    private javax.swing.JTextField timestepTextField;
    private javax.swing.JLabel timestepunitlabel;
    private javax.swing.JTextField tolerance;
    private javax.swing.JPanel topjpanel;
    private javax.swing.JScrollPane topjscroll;
    public javax.swing.JComboBox typeofensembleComboBox;
    public javax.swing.JComboBox typeofsimulationComboBox;
    private javax.swing.JTextField versiontxt;
    private javax.swing.JButton viewparbtn;
    private javax.swing.JButton viewresbtn;
    private javax.swing.JTextField visualfreqTextField;
    private javax.swing.JTextField walltimetxt;
    // End of variables declaration//GEN-END:variables
    
}
