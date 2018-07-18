/*
 * Ms2chart.java
 *
 * Created on October 24, 2006, 1:41 PM
 *
 * Copyright (C) 2006, Anupam Srivastava
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */

package ms_newcharts;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Component;
import java.awt.Cursor;
import java.awt.Desktop;
import java.awt.GraphicsDevice;
import java.awt.GraphicsEnvironment;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.StringTokenizer;
import javax.imageio.ImageIO;
import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * The Chart Program for MS2 output files.
 * Released under LGPL. For complete license, see copying.txt. Alternatively visit: http://www.gnu.org/licenses/lgpl.html
 * <P>
 * NOTE:<BR>
 * 1. The panel which holds the chart has to be <B>BorderLayout</B>. In <I>NetBeans 5.0</I>, right-click on the panel and select <B>BorderLayout</B>.<BR>
 * 2. It is important to note that RUN files are treated first than RAV file. This is to simplify the algorithm.<BR>
 * </P>
 * @author Anupam Srivastava, Fri, 24 Nov 2006 12:36:19 +0100
 * @version 2.3beta
 */
public class Ms2chart extends javax.swing.JFrame {
    /**
     * Returns the index of the largest integer in the given array of type <CODE>Integer</CODE>.
     * <P>For example,</P>
     * <P><BLOCKQUOTE><CODE>
     * int[] a1 = { 0, 1, 2, 100, 4, 50, 22 };<BR>
     * int i = findIndexOfLargestInt(a1); // i = 2<BR>
     * System.out.println(a1[i]); // ai[i] = 100<BR>
     * </CODE></BLOCKQUOTE></P>
     * @param arrayOfInts Array of the set of integers.
     * @return The index of the largest integer in the input array <B>arrayOfInts</B>.
     */
  
    public int findIndexOfLargestInt(int[] arrayOfInts) {
        int[] tempArray = new int[arrayOfInts.length];
        System.arraycopy(arrayOfInts, 0, tempArray, 0, arrayOfInts.length);
        Arrays.sort(tempArray);
        int ArrLength = tempArray.length;
        if(tempArray.length > 3)
        {
            ArrLength = ArrLength - 1;
        }
        
        return binarySearch(arrayOfInts, tempArray[ArrLength - 1]);
        
    }
       
    /**
     * Performs the standard binary search using two comparisons per level.
     * @param arrayOfInts
     * @param intToBeFound
     * @return The index where the item is found, or -1.
     */
    public static int binarySearch(int[] arrayOfInts, int intToBeFound) {
        int temphigh = arrayOfInts.length;
        if(arrayOfInts.length > 3)
        {
            temphigh = temphigh - 1;
        }
        int low = 0;
        int high = temphigh - 1;
        int mid;
        //JOptionPane.showMessageDialog(null,intToBeFound);
        while(low <= high) {
           // JOptionPane.showMessageDialog(null,low +" | "+ high +" =:in loop");
            mid = (low + high) / 2;
           // JOptionPane.showMessageDialog(null,mid + "=:mid var");
           // JOptionPane.showMessageDialog(null,arrayOfInts[mid] + "=:array of ints");
            if(arrayOfInts[mid] > intToBeFound) {
                
                low = mid + 1;
                //JOptionPane.showMessageDialog(null,low +"=low :in mid");
               // JOptionPane.showMessageDialog(null,arrayOfInts[mid] + " :mid var");
            } else if(arrayOfInts[mid] < intToBeFound) {
                high = mid - 1;
               // JOptionPane.showMessageDialog(null,high +"=high :in mid");
            } else {
                //JOptionPane.showMessageDialog(null,mid +"=return mid");
                return mid;
            }
        }
        //JOptionPane.showMessageDialog(null,":out of loop");
        return -1;
    }
    
    /**
     * Does initialization.
     */
    public Ms2chart() {
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
        UserDirectory = System.getProperty("user.dir");
        selectedDirectory = startingDirectory;
        comboCaseList.setModel(new javax.swing.DefaultComboBoxModel(getDirList()));
        if (isHavingCases) {
            prepareData();
        } else {
            cbRUN.setEnabled(false);
            cbRAV.setEnabled(false);
            cbNVTE.setEnabled(false);
            cbEQUI.setEnabled(false);
            cbPROD.setEnabled(false);
            tfXMax.setEnabled(false);
            tfXMin.setEnabled(false);
            tfYMax.setEnabled(false);
            tfYMin.setEnabled(false);
            tfError.setEnabled(false);
            buttonLOG.setEnabled(false);
            buttonPAR.setEnabled(false);
            buttonRES.setEnabled(false);
            buttonDrawChart.setEnabled(false);
            buttonViewData.setEnabled(false);
            buttonPNGsave.setEnabled(false);
            buttonXAxisChange.setEnabled(false);
            buttonYAxisChange.setEnabled(false);
            buttonClear.setEnabled(false);
            comboXAxis.setModel(new javax.swing.DefaultComboBoxModel(new String[] {""}));
            listYAxises.setModel(new javax.swing.AbstractListModel() {
                String[] strings = new String[] {""};
                public int getSize() { return strings.length; }
                public Object getElementAt(int i) { return strings[i]; }
            });
            setLabelStatusMessage("Browse for a directory...");
            panelChartContainer.removeAll();
            panelChartContainer.setVisible(false);
            panelChartContainer.setVisible(true);
        }
        setwindowsize();
        readHelpxlsxfile();
        settooltipforlabels();
        addShortCutsForButtonActions();
    }
    private void addShortCutsForButtonActions() {

    //Shortcut for load button
    Action loadbuttonAction = new AbstractAction("Load") {
        @Override
        public void actionPerformed(ActionEvent evt) {
           buttonBrowseActionPerformed(evt);
        }
    };
 
    this.buttonBrowse.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(
        KeyStroke.getKeyStroke(KeyEvent.VK_A, KeyEvent.CTRL_DOWN_MASK), "Load");
    this.buttonBrowse.getActionMap().put("Load", loadbuttonAction);
    
   //shortcut for help button
    Action helpbuttonAction = new AbstractAction("Help") {

        @Override
        public void actionPerformed(ActionEvent evt) {
            HelpbtnActionPerformed(evt);
        }
    };

    this.Helpbtn.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(
            KeyStroke.getKeyStroke(KeyEvent.VK_F1, 0), "Help");
    this.Helpbtn.getActionMap().put("Help", helpbuttonAction);

}
    /**
     * Get all the label and buttons from the UI
     * Assign tool tip from the populated hashMap (helpmap)
     * Call addlinebreaks function to make tool tip in paragraph form
     * use html tags inside tool tip to formate paragrah form.
    */
    private void settooltipforlabels()
    {
        for(Component cm: main_panel.getComponents())
        {
                if(cm instanceof JPanel)
                {
                    JPanel jp = (JPanel) cm;
                    for(Component cmp: jp.getComponents())
                    {
                        if(cmp instanceof JLabel)
                        {
                            JLabel lbl = (JLabel) cmp;
                            if(helpMap.containsKey(lbl.getText().toLowerCase().trim()))
                            {
                            String tooltip = addLinebreaks(helpMap.get(lbl.getText().toLowerCase().trim()),30);
                            lbl.setToolTipText("<html><body><p>"+tooltip+"</p></body></html>");
                            }
                        }
                        
                        if(cmp instanceof JButton)
                        {
                            JButton btn = (JButton) cmp;
                            if(helpMap.containsKey(btn.getText().toLowerCase().trim()))
                            {
                            String tooltip = addLinebreaks(helpMap.get(btn.getText().toLowerCase().trim()),30);
                            btn.setToolTipText("<html><body><p>"+tooltip+"</p></body></html>");
                            }
                        }
                        else
                        if(cmp instanceof JCheckBox)
                        {
                            JCheckBox chk = (JCheckBox) cmp;
                            if(helpMap.containsKey(chk.getText().toLowerCase().trim()))
                            {
                            String tooltip = addLinebreaks(helpMap.get(chk.getText().toLowerCase().trim()),30);
                            chk.setToolTipText("<html><body><p>"+tooltip+"</p></body></html>");
                            }
                        }
                        else
                        if(cmp instanceof JRadioButton)
                        {
                            JRadioButton rbtn = (JRadioButton) cmp;
                            if(helpMap.containsKey(rbtn.getText().toLowerCase().trim()))
                            {
                            String tooltip = addLinebreaks(helpMap.get(rbtn.getText().toLowerCase().trim()),30);
                            rbtn.setToolTipText("<html><body><p>"+tooltip+"</p></body></html>");
                            }
                        }
                    }
                }
        }
    }

    /**
     * To add line breaks after every specified characters length
     * Input: String which needed to be break in line
     * maxLineLength: maximum length specified to break the line
     * 
     * @param input Help string
     * @param maxLineLength Length of the line
     * 
     * @return Formatted String
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
      * Read help file from the help folder
      * Find the "Name" & "Meaning" coloumn inside the help xlsx file (we need only these two)
      * populate the hashmap (helpMap) to save all the names and meanings
      *     Key: Text from the "Name" column
      *     Value: Text from the "Meaning"" column.
      * 
     */
     private void readHelpxlsxfile()
    {
        try 
        {
        // Get the workbook instance for XLS file
        String filePath = UserDirectory+""+File.separator+"help"+File.separator+"ms2chart_help_format.xlsx";
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
     
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        buttonGroupFileType = new javax.swing.ButtonGroup();
        mainscroll = new javax.swing.JScrollPane();
        main_panel = new javax.swing.JPanel();
        panelChartContainer = new javax.swing.JPanel();
        panelStatusBar = new javax.swing.JPanel();
        jLabel1 = new javax.swing.JLabel();
        labelStatus = new javax.swing.JLabel();
        buttonViewData = new javax.swing.JButton();
        Helpbtn = new javax.swing.JButton();
        panelControl = new javax.swing.JPanel();
        jLabel3 = new javax.swing.JLabel();
        buttonBrowse = new javax.swing.JButton();
        comboCaseList = new javax.swing.JComboBox();
        cbRUN = new javax.swing.JCheckBox();
        cbRAV = new javax.swing.JCheckBox();
        jSeparator1 = new javax.swing.JSeparator();
        cbNVTE = new javax.swing.JCheckBox();
        cbEQUI = new javax.swing.JCheckBox();
        cbPROD = new javax.swing.JCheckBox();
        jLabel4 = new javax.swing.JLabel();
        comboXAxis = new javax.swing.JComboBox();
        jLabel5 = new javax.swing.JLabel();
        jScrollPane1 = new javax.swing.JScrollPane();
        listYAxises = new javax.swing.JList();
        jLabel6 = new javax.swing.JLabel();
        tfXMax = new javax.swing.JTextField();
        tfXMin = new javax.swing.JTextField();
        buttonXAxisChange = new javax.swing.JButton();
        jLabel7 = new javax.swing.JLabel();
        tfYMax = new javax.swing.JTextField();
        tfYMin = new javax.swing.JTextField();
        buttonYAxisChange = new javax.swing.JButton();
        jLabel8 = new javax.swing.JLabel();
        tfError = new javax.swing.JTextField();
        buttonDrawChart = new javax.swing.JButton();
        buttonPNGsave = new javax.swing.JButton();
        cbAutoDraw = new javax.swing.JCheckBox();
        cbAskDirectory = new javax.swing.JCheckBox();
        jLabel9 = new javax.swing.JLabel();
        buttonLOG = new javax.swing.JButton();
        buttonPAR = new javax.swing.JButton();
        buttonRES = new javax.swing.JButton();
        jSeparator2 = new javax.swing.JSeparator();
        buttonRestart = new javax.swing.JButton();
        buttonClear = new javax.swing.JButton();
        rbRUN_RAV = new javax.swing.JRadioButton();
        rbRTR = new javax.swing.JRadioButton();
        jSeparator3 = new javax.swing.JSeparator();
        jLabel10 = new javax.swing.JLabel();
        vimtovmd_btn = new javax.swing.JButton();
        buttonQuit1 = new javax.swing.JButton();
        rbTHI = new javax.swing.JRadioButton();
        rbDCP = new javax.swing.JRadioButton();
        jLabel11 = new javax.swing.JLabel();
        jLabel12 = new javax.swing.JLabel();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setTitle("Ms2 Chart");
        setName("ms2chart_frame"); // NOI18N

        panelChartContainer.setLayout(new java.awt.BorderLayout());

        jLabel1.setText("Status:");

        buttonViewData.setText("View data");
        buttonViewData.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                buttonViewDataActionPerformed(evt);
            }
        });

        Helpbtn.setText("Help");
        Helpbtn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                HelpbtnActionPerformed(evt);
            }
        });

        org.jdesktop.layout.GroupLayout panelStatusBarLayout = new org.jdesktop.layout.GroupLayout(panelStatusBar);
        panelStatusBar.setLayout(panelStatusBarLayout);
        panelStatusBarLayout.setHorizontalGroup(
            panelStatusBarLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(panelStatusBarLayout.createSequentialGroup()
                .addContainerGap()
                .add(jLabel1)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(labelStatus)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .add(buttonViewData)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(Helpbtn, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 97, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );
        panelStatusBarLayout.setVerticalGroup(
            panelStatusBarLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(panelStatusBarLayout.createSequentialGroup()
                .addContainerGap(org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .add(panelStatusBarLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(org.jdesktop.layout.GroupLayout.TRAILING, panelStatusBarLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                        .add(jLabel1)
                        .add(labelStatus))
                    .add(org.jdesktop.layout.GroupLayout.TRAILING, panelStatusBarLayout.createSequentialGroup()
                        .add(panelStatusBarLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                            .add(buttonViewData)
                            .add(Helpbtn, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addContainerGap())))
        );

        panelControl.setPreferredSize(new java.awt.Dimension(260, 539));

        jLabel3.setText("Cases:");

        buttonBrowse.setText("Browse");
        buttonBrowse.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                buttonBrowseActionPerformed(evt);
            }
        });

        comboCaseList.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                comboCaseListActionPerformed(evt);
            }
        });

        cbRUN.setText(".run");
        cbRUN.setBorder(javax.swing.BorderFactory.createEmptyBorder(0, 0, 0, 0));
        cbRUN.setMargin(new java.awt.Insets(0, 0, 0, 0));
        cbRUN.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                cbRUNItemStateChanged(evt);
            }
        });

        cbRAV.setSelected(true);
        cbRAV.setText(".rav");
        cbRAV.setBorder(javax.swing.BorderFactory.createEmptyBorder(0, 0, 0, 0));
        cbRAV.setMargin(new java.awt.Insets(0, 0, 0, 0));
        cbRAV.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                cbRAVItemStateChanged(evt);
            }
        });

        cbNVTE.setText("NVT-Equilibriation");
        cbNVTE.setBorder(javax.swing.BorderFactory.createEmptyBorder(0, 0, 0, 0));
        cbNVTE.setMargin(new java.awt.Insets(0, 0, 0, 0));
        cbNVTE.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                cbNVTEItemStateChanged(evt);
            }
        });

        cbEQUI.setText("NPT-Equilibration");
        cbEQUI.setBorder(javax.swing.BorderFactory.createEmptyBorder(0, 0, 0, 0));
        cbEQUI.setMargin(new java.awt.Insets(0, 0, 0, 0));
        cbEQUI.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                cbEQUIItemStateChanged(evt);
            }
        });

        cbPROD.setSelected(true);
        cbPROD.setText("Production");
        cbPROD.setBorder(javax.swing.BorderFactory.createEmptyBorder(0, 0, 0, 0));
        cbPROD.setMargin(new java.awt.Insets(0, 0, 0, 0));
        cbPROD.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                cbPRODItemStateChanged(evt);
            }
        });

        jLabel4.setText("X-Axis");

        comboXAxis.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                comboXAxisActionPerformed(evt);
            }
        });

        jLabel5.setText("Y-Axis");

        listYAxises.addListSelectionListener(new javax.swing.event.ListSelectionListener() {
            public void valueChanged(javax.swing.event.ListSelectionEvent evt) {
                listYAxisesValueChanged(evt);
            }
        });
        jScrollPane1.setViewportView(listYAxises);

        jLabel6.setText("X_max");

        tfXMax.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                tfXMaxActionPerformed(evt);
            }
        });

        tfXMin.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                tfXMinActionPerformed(evt);
            }
        });

        buttonXAxisChange.setText("OK");
        buttonXAxisChange.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                buttonXAxisChangeActionPerformed(evt);
            }
        });

        jLabel7.setText("Y_max");

        tfYMax.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                tfYMaxActionPerformed(evt);
            }
        });

        tfYMin.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                tfYMinActionPerformed(evt);
            }
        });

        buttonYAxisChange.setText("OK");
        buttonYAxisChange.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                buttonYAxisChangeActionPerformed(evt);
            }
        });

        jLabel8.setText("Statistical error");

        tfError.setEditable(false);

        buttonDrawChart.setText("Draw Chart");
        buttonDrawChart.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                buttonDrawChartActionPerformed(evt);
            }
        });

        buttonPNGsave.setText("Save Image");
        buttonPNGsave.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                buttonPNGsaveActionPerformed(evt);
            }
        });

        cbAutoDraw.setSelected(true);
        cbAutoDraw.setText("Auto draw");
        cbAutoDraw.setBorder(javax.swing.BorderFactory.createEmptyBorder(0, 0, 0, 0));
        cbAutoDraw.setMargin(new java.awt.Insets(0, 0, 0, 0));
        cbAutoDraw.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                cbAutoDrawItemStateChanged(evt);
            }
        });

        cbAskDirectory.setSelected(true);
        cbAskDirectory.setText("Ask directory");
        cbAskDirectory.setBorder(javax.swing.BorderFactory.createEmptyBorder(0, 0, 0, 0));
        cbAskDirectory.setMargin(new java.awt.Insets(0, 0, 0, 0));
        cbAskDirectory.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                cbAskDirectoryItemStateChanged(evt);
            }
        });

        jLabel9.setText("View:");

        buttonLOG.setText("log");
        buttonLOG.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                buttonLOGActionPerformed(evt);
            }
        });

        buttonPAR.setText("par");
        buttonPAR.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                buttonPARActionPerformed(evt);
            }
        });

        buttonRES.setText("res");
        buttonRES.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                buttonRESActionPerformed(evt);
            }
        });

        buttonRestart.setText("Restart");
        buttonRestart.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                buttonRestartActionPerformed(evt);
            }
        });

        buttonClear.setText("Clear");
        buttonClear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                buttonClearActionPerformed(evt);
            }
        });

        buttonGroupFileType.add(rbRUN_RAV);
        rbRUN_RAV.setSelected(true);
        rbRUN_RAV.setText("run/rav");
        rbRUN_RAV.setBorder(javax.swing.BorderFactory.createEmptyBorder(0, 0, 0, 0));
        rbRUN_RAV.setMargin(new java.awt.Insets(0, 0, 0, 0));
        rbRUN_RAV.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                rbRUN_RAVActionPerformed(evt);
            }
        });

        buttonGroupFileType.add(rbRTR);
        rbRTR.setText("rtr");
        rbRTR.setBorder(javax.swing.BorderFactory.createEmptyBorder(0, 0, 0, 0));
        rbRTR.setMargin(new java.awt.Insets(0, 0, 0, 0));
        rbRTR.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                rbRTRActionPerformed(evt);
            }
        });

        jLabel10.setText("Conversion to VMD:");

        vimtovmd_btn.setText("Convert");
        vimtovmd_btn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                vimtovmd_btnActionPerformed(evt);
            }
        });

        buttonQuit1.setText("Quit");
        buttonQuit1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                buttonQuit1ActionPerformed(evt);
            }
        });

        buttonGroupFileType.add(rbTHI);
        rbTHI.setText("thi");
        rbTHI.setBorder(javax.swing.BorderFactory.createEmptyBorder(0, 0, 0, 0));
        rbTHI.setMargin(new java.awt.Insets(0, 0, 0, 0));
        rbTHI.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                rbTHIActionPerformed(evt);
            }
        });

        buttonGroupFileType.add(rbDCP);
        rbDCP.setText("dcp");
        rbDCP.setBorder(javax.swing.BorderFactory.createEmptyBorder(0, 0, 0, 0));
        rbDCP.setMargin(new java.awt.Insets(0, 0, 0, 0));
        rbDCP.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                rbDCPActionPerformed(evt);
            }
        });

        jLabel11.setText("X_min");

        jLabel12.setText("Y_min");

        org.jdesktop.layout.GroupLayout panelControlLayout = new org.jdesktop.layout.GroupLayout(panelControl);
        panelControl.setLayout(panelControlLayout);
        panelControlLayout.setHorizontalGroup(
            panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(org.jdesktop.layout.GroupLayout.TRAILING, panelControlLayout.createSequentialGroup()
                .addContainerGap()
                .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, jSeparator1)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, panelControlLayout.createSequentialGroup()
                        .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING, false)
                            .add(panelControlLayout.createSequentialGroup()
                                .add(jLabel7)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                .add(tfYMax, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                            .add(org.jdesktop.layout.GroupLayout.TRAILING, panelControlLayout.createSequentialGroup()
                                .add(jLabel6)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .add(tfXMax, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 43, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)))
                        .add(6, 6, 6)
                        .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING, false)
                            .add(panelControlLayout.createSequentialGroup()
                                .add(jLabel12)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .add(tfYMin, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                            .add(panelControlLayout.createSequentialGroup()
                                .add(jLabel11)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                .add(tfXMin, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 64, Short.MAX_VALUE)))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                            .add(buttonXAxisChange, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .add(buttonYAxisChange, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                    .add(org.jdesktop.layout.GroupLayout.LEADING, panelControlLayout.createSequentialGroup()
                        .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING, false)
                            .add(jLabel5, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .add(jLabel4, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                            .add(comboXAxis, 0, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .add(jScrollPane1, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 0, Short.MAX_VALUE)))
                    .add(org.jdesktop.layout.GroupLayout.LEADING, comboCaseList, 0, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, panelControlLayout.createSequentialGroup()
                        .add(jLabel3, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 88, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .add(buttonBrowse))
                    .add(org.jdesktop.layout.GroupLayout.LEADING, panelControlLayout.createSequentialGroup()
                        .add(buttonQuit1, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 63, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING)
                            .add(panelControlLayout.createSequentialGroup()
                                .add(0, 0, Short.MAX_VALUE)
                                .add(buttonLOG)
                                .add(18, 18, 18)
                                .add(buttonPAR, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 63, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                                .add(18, 18, 18)
                                .add(buttonRES, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 64, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                            .add(panelControlLayout.createSequentialGroup()
                                .add(buttonClear)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .add(buttonRestart, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 81, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))))
                    .add(org.jdesktop.layout.GroupLayout.LEADING, panelControlLayout.createSequentialGroup()
                        .add(12, 12, 12)
                        .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                            .add(cbRUN)
                            .add(cbRAV))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING, false)
                            .add(org.jdesktop.layout.GroupLayout.TRAILING, rbRTR, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .add(org.jdesktop.layout.GroupLayout.TRAILING, rbRUN_RAV, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .add(rbTHI, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .add(rbDCP, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                    .add(org.jdesktop.layout.GroupLayout.LEADING, jSeparator2)
                    .add(panelControlLayout.createSequentialGroup()
                        .add(jLabel10)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .add(vimtovmd_btn, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 87, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                    .add(org.jdesktop.layout.GroupLayout.LEADING, jSeparator3)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, panelControlLayout.createSequentialGroup()
                        .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING)
                            .add(org.jdesktop.layout.GroupLayout.LEADING, jLabel9)
                            .add(org.jdesktop.layout.GroupLayout.LEADING, panelControlLayout.createSequentialGroup()
                                .add(12, 12, 12)
                                .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                                    .add(cbPROD)
                                    .add(cbEQUI)
                                    .add(cbNVTE)))
                            .add(org.jdesktop.layout.GroupLayout.LEADING, panelControlLayout.createSequentialGroup()
                                .add(buttonPNGsave)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                .add(cbAskDirectory))
                            .add(org.jdesktop.layout.GroupLayout.LEADING, panelControlLayout.createSequentialGroup()
                                .add(buttonDrawChart)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                .add(cbAutoDraw)))
                        .add(0, 0, Short.MAX_VALUE))
                    .add(org.jdesktop.layout.GroupLayout.LEADING, panelControlLayout.createSequentialGroup()
                        .add(jLabel8)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(tfError)))
                .add(45, 45, 45))
        );

        panelControlLayout.linkSize(new java.awt.Component[] {tfXMax, tfXMin, tfYMax, tfYMin}, org.jdesktop.layout.GroupLayout.HORIZONTAL);

        panelControlLayout.linkSize(new java.awt.Component[] {buttonDrawChart, buttonPNGsave}, org.jdesktop.layout.GroupLayout.HORIZONTAL);

        panelControlLayout.setVerticalGroup(
            panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(panelControlLayout.createSequentialGroup()
                .addContainerGap()
                .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel3)
                    .add(buttonBrowse))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(comboCaseList, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(cbRUN)
                    .add(rbRUN_RAV))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(cbRAV)
                    .add(rbRTR))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(rbTHI)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(rbDCP)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(jSeparator1, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 18, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(cbNVTE)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(cbEQUI)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(cbPROD)
                .add(18, 18, 18)
                .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(comboXAxis, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(jLabel4))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(jScrollPane1, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 102, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(jLabel5))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.UNRELATED)
                .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(buttonXAxisChange)
                    .add(tfXMin, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(tfXMax, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(jLabel6)
                    .add(jLabel11))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(buttonYAxisChange)
                    .add(tfYMin, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(tfYMax, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(jLabel7)
                    .add(jLabel12))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.UNRELATED)
                .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel8)
                    .add(tfError, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .add(18, 18, 18)
                .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(buttonDrawChart)
                    .add(cbAutoDraw))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(buttonPNGsave)
                    .add(cbAskDirectory))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.UNRELATED)
                .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel9)
                    .add(buttonRES)
                    .add(buttonPAR)
                    .add(buttonLOG))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(jSeparator2, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 11, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel10)
                    .add(vimtovmd_btn))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(jSeparator3, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 12, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(panelControlLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(buttonQuit1, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 35, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(buttonClear, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 35, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(buttonRestart, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 35, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(151, Short.MAX_VALUE))
        );

        panelControlLayout.linkSize(new java.awt.Component[] {tfXMax, tfXMin, tfYMax, tfYMin}, org.jdesktop.layout.GroupLayout.VERTICAL);

        panelControlLayout.linkSize(new java.awt.Component[] {buttonXAxisChange, buttonYAxisChange}, org.jdesktop.layout.GroupLayout.VERTICAL);

        panelControlLayout.linkSize(new java.awt.Component[] {buttonDrawChart, buttonPNGsave}, org.jdesktop.layout.GroupLayout.VERTICAL);

        org.jdesktop.layout.GroupLayout main_panelLayout = new org.jdesktop.layout.GroupLayout(main_panel);
        main_panel.setLayout(main_panelLayout);
        main_panelLayout.setHorizontalGroup(
            main_panelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(main_panelLayout.createSequentialGroup()
                .addContainerGap()
                .add(main_panelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(panelStatusBar, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(main_panelLayout.createSequentialGroup()
                        .add(panelControl, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 321, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                        .add(18, 18, 18)
                        .add(panelChartContainer, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 989, Short.MAX_VALUE)))
                .add(25, 25, 25))
        );
        main_panelLayout.setVerticalGroup(
            main_panelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(main_panelLayout.createSequentialGroup()
                .add(5, 5, 5)
                .add(main_panelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(panelChartContainer, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(panelControl, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 1001, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(panelStatusBar, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );

        mainscroll.setViewportView(main_panel);

        org.jdesktop.layout.GroupLayout layout = new org.jdesktop.layout.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(mainscroll)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(mainscroll)
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents
    
    /**
     * 
     * Set the initial size of the window according to the resolution of the screen.
     * 
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
        
        if(screenwidth < Pwidth && screenheight < Pheight)
        {
         parentframe.setSize(screenwidth - 10, screenheight - 50);
         
        }
        else
        if(screenwidth < Pwidth)
            {
             parentframe.setSize(screenwidth - 10, Pheight);
            }
            else
            if(screenheight < Pheight)
            {
               parentframe.setSize(Pwidth, screenheight - 50);
   
            }
        mainscroll.getVerticalScrollBar().setUnitIncrement(16);
    }
    
    /**
     * RTR radio button event handler
     * Enable and disable UI controls depending upon condition
     * 
     * @see Ms2chart#prepareData()
     * 
     * @param evt 
     */
    private void rbRTRActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_rbRTRActionPerformed
        panelChartContainer.removeAll();
        setLabelStatusMessage("Listing .rtr files");
        comboCaseList.setModel(new javax.swing.DefaultComboBoxModel(getDirList()));
        if (isHavingCases) {
            prepareData();
        } else {
            cbRUN.setEnabled(false);
            cbRAV.setEnabled(false);
            cbNVTE.setEnabled(false);
            cbEQUI.setEnabled(false);
            cbPROD.setEnabled(false);
            tfXMax.setEnabled(false);
            tfXMin.setEnabled(false);
            tfYMax.setEnabled(false);
            tfYMin.setEnabled(false);
            tfError.setEnabled(false);
            buttonLOG.setEnabled(false);
            buttonPAR.setEnabled(false);
            buttonRES.setEnabled(false);
            buttonDrawChart.setEnabled(false);
            buttonViewData.setEnabled(false);
            buttonPNGsave.setEnabled(false);
            buttonXAxisChange.setEnabled(false);
            buttonYAxisChange.setEnabled(false);
            buttonClear.setEnabled(false);
            comboXAxis.setModel(new javax.swing.DefaultComboBoxModel(new String[] {""}));
            listYAxises.setModel(new javax.swing.AbstractListModel() {
                String[] strings = new String[] {""};
                public int getSize() { return strings.length; }
                public Object getElementAt(int i) { return strings[i]; }
            });
            setLabelStatusMessage("No file found. Browse again...");
            panelChartContainer.removeAll();
            panelChartContainer.setVisible(false);
            panelChartContainer.setVisible(true);
        }
    }//GEN-LAST:event_rbRTRActionPerformed
    
    /**
     * RUN radio button event handler
     * Enable and disable UI controls depending upon condition.
     * 
     * @see Ms2chart#prepareData()
     * 
     * @param evt 
     */
    private void rbRUN_RAVActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_rbRUN_RAVActionPerformed
        panelChartContainer.removeAll();
        setLabelStatusMessage("Listing .run and .rav files");
        comboCaseList.setModel(new javax.swing.DefaultComboBoxModel(getDirList()));
        if (isHavingCases) {
            prepareData();
        } else {
            cbRUN.setEnabled(false);
            cbRAV.setEnabled(false);
            cbNVTE.setEnabled(false);
            cbEQUI.setEnabled(false);
            cbPROD.setEnabled(false);
            tfXMax.setEnabled(false);
            tfXMin.setEnabled(false);
            tfYMax.setEnabled(false);
            tfYMin.setEnabled(false);
            tfError.setEnabled(false);
            buttonLOG.setEnabled(false);
            buttonPAR.setEnabled(false);
            buttonRES.setEnabled(false);
            buttonDrawChart.setEnabled(false);
            buttonViewData.setEnabled(false);
            buttonPNGsave.setEnabled(false);
            buttonXAxisChange.setEnabled(false);
            buttonYAxisChange.setEnabled(false);
            buttonClear.setEnabled(false);
            comboXAxis.setModel(new javax.swing.DefaultComboBoxModel(new String[] {""}));
            listYAxises.setModel(new javax.swing.AbstractListModel() {
                String[] strings = new String[] {""};
                public int getSize() { return strings.length; }
                public Object getElementAt(int i) { return strings[i]; }
            });
            setLabelStatusMessage("No file found. Browse again...");
            panelChartContainer.removeAll();
            panelChartContainer.setVisible(false);
            panelChartContainer.setVisible(true);
        }
    }//GEN-LAST:event_rbRUN_RAVActionPerformed
    
    /**
     * @see Ms2chart#buttonYAxisChangeActionPerformed(java.awt.event.ActionEvent)
     * 
     * @param evt 
     */
    private void tfYMinActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_tfYMinActionPerformed
        buttonYAxisChangeActionPerformed(evt);
    }//GEN-LAST:event_tfYMinActionPerformed
    
    /**
     * @see Ms2chart#buttonYAxisChangeActionPerformed(java.awt.event.ActionEvent)
     * 
     * @param evt 
     */
    private void tfYMaxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_tfYMaxActionPerformed
        buttonYAxisChangeActionPerformed(evt);
    }//GEN-LAST:event_tfYMaxActionPerformed
    
    /**
     * @see Ms2chart#buttonXAxisChangeActionPerformed(java.awt.event.ActionEvent)
     * 
     * @param evt 
     */
    private void tfXMinActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_tfXMinActionPerformed
        buttonXAxisChangeActionPerformed(evt);
    }//GEN-LAST:event_tfXMinActionPerformed
    
    /**
     * @see Ms2chart#buttonXAxisChangeActionPerformed(java.awt.event.ActionEvent)
     * 
     * @param evt 
     */
    private void tfXMaxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_tfXMaxActionPerformed
        buttonXAxisChangeActionPerformed(evt);
    }//GEN-LAST:event_tfXMaxActionPerformed
    
    /**
     * Checkbox production event handler.
     * 
     * @see Ms2chart#plot()
     * 
     * @param evt 
     */
    private void cbPRODItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_cbPRODItemStateChanged
        if((columnInRunFileForXAxis != -1 || columnInRavFileForXAxis != -1) && (columnsInRunFileForYAxis.length != 0 || columnsInRavFileForYAxis.length != 0) && (cbRUN.isEnabled() && cbRUN.isSelected() || cbRAV.isEnabled() && cbRAV.isSelected()) && (cbNVTE.isEnabled() && cbNVTE.isSelected() || cbEQUI.isEnabled() && cbEQUI.isSelected() || cbPROD.isEnabled() && cbPROD.isSelected())) {
            buttonDrawChart.setEnabled(true);
            if(cbAutoDraw.isSelected())
                plot();
        } else {
            buttonDrawChart.setEnabled(false);
        }
    }//GEN-LAST:event_cbPRODItemStateChanged
    
    /**
     * Checkbox NPT - Equilibrium event handler.
     * 
     * @see Ms2chart#plot()
     * 
     * @param evt 
     */
    private void cbEQUIItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_cbEQUIItemStateChanged
        if((columnInRunFileForXAxis != -1 || columnInRavFileForXAxis != -1) && (columnsInRunFileForYAxis.length != 0 || columnsInRavFileForYAxis.length != 0) && (cbRUN.isEnabled() && cbRUN.isSelected() || cbRAV.isEnabled() && cbRAV.isSelected()) && (cbNVTE.isEnabled() && cbNVTE.isSelected() || cbEQUI.isEnabled() && cbEQUI.isSelected() || cbPROD.isEnabled() && cbPROD.isSelected())) {
            buttonDrawChart.setEnabled(true);
            if(cbAutoDraw.isSelected())
                plot();
        } else {
            buttonDrawChart.setEnabled(false);
        }
    }//GEN-LAST:event_cbEQUIItemStateChanged
    
    /**
     * Checkbox NVT - Equilibrium event handler.
     * 
     * @see Ms2chart#plot()
     * 
     * @param evt 
     */
    private void cbNVTEItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_cbNVTEItemStateChanged
        if((columnInRunFileForXAxis != -1 || columnInRavFileForXAxis != -1) && (columnsInRunFileForYAxis.length != 0 || columnsInRavFileForYAxis.length != 0) && (cbRUN.isEnabled() && cbRUN.isSelected() || cbRAV.isEnabled() && cbRAV.isSelected()) && (cbNVTE.isEnabled() && cbNVTE.isSelected() || cbEQUI.isEnabled() && cbEQUI.isSelected() || cbPROD.isEnabled() && cbPROD.isSelected())) {
            buttonDrawChart.setEnabled(true);
            if(cbAutoDraw.isSelected())
                plot();
        } else {
            buttonDrawChart.setEnabled(false);
        }
    }//GEN-LAST:event_cbNVTEItemStateChanged
    
    /**
     * Case dropdown event handler.
     * 
     * @see Ms2chart#prepareData()
     * 
     * @param evt 
     */
    private void comboCaseListActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_comboCaseListActionPerformed
        if (isHavingCases) {
            cbEQUI.setSelected(false);
            cbEQUI.setEnabled(false);
            prepareData();
        }
    }//GEN-LAST:event_comboCaseListActionPerformed
    
    /**
     * x-Axis drop down event handler
     * Checks different user input condition
     * set the title for x-axis
     * and calls plot function accordingly
     * 
     * @see Ms2chart#plot()
     * @see Ms2chart#settitleXAxis(java.lang.String)
     * 
     * @param evt 
     */
    String titleXAxis = null;
    private void comboXAxisActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_comboXAxisActionPerformed
        String selectedValue = (String) comboXAxis.getSelectedItem();
        columnInRtrFileForXAxis = -1;
        columnInRunFileForXAxis = -1;
        columnInRavFileForXAxis = -1;
        String title = (String) selectedValue;
        
        if (rbRTR.isSelected() || rbTHI.isSelected() || rbDCP.isSelected()) {/// i added rbTHI check with rtr
            for (int j = 0; j < arrayOfRtrFile[0].length; ++j) {
                if (title.equals(arrayOfRtrFile[0][j])) {
                    columnInRtrFileForXAxis = j; 
                    settitleXAxis(title);
                }
            }
        }
        if (cbRUN.isEnabled()) {
            for (int j = 0; j < arrayOfRunFile[0].length; ++j) {
                if (title.equals(arrayOfRunFile[0][j])) {
                    columnInRunFileForXAxis = j;
                    settitleXAxis(title);
                }
            }
        }
        if (cbRAV.isEnabled()) {
            for (int j = 0; j < arrayOfRavFile[0].length; ++j) {
                if (title.equals(arrayOfRavFile[0][j])) {
                    columnInRavFileForXAxis = j;
                    settitleXAxis(title);
                    
                }
            }
        }
        if (columnInRunFileForXAxis == -1) {
            cbRUN.setEnabled(false);
        }
        if (columnInRavFileForXAxis == -1) {
            cbRAV.setEnabled(false);
        }
        if (rbRTR.isSelected() || rbTHI.isSelected() || rbDCP.isSelected()) { /////// i added rbTHI check with rtr
            if (columnInRtrFileForXAxis != -1
                    && columnsInRtrFileForYAxis.length != 0) {
                buttonDrawChart.setEnabled(true);
                if (cbAutoDraw.isSelected()) {
                    plot();
                        if(rbTHI.isSelected() && selectedValue.trim().equals("LAMBDA"))
                           {
                                domainAxis.setUpperBound(1.0);
                                domainAxis.setLowerBound(0.0);
                                tfXMax.setText(Math.round(domainAxis.getRange().getUpperBound() * 100) / 100.0 + "");
                                tfXMin.setText(Math.round(domainAxis.getRange().getLowerBound() * 100) / 100.0 + "");
                           }
                }
            } else {
                buttonDrawChart.setEnabled(false);
            }
        } else {
            if ((columnInRunFileForXAxis != -1
                    || columnInRavFileForXAxis != -1)
                    && (columnsInRunFileForYAxis.length != 0
                    || columnsInRavFileForYAxis.length != 0)
                    && ((cbRUN.isEnabled()
                    && cbRUN.isSelected())
                    || (cbRAV.isEnabled()
                    && cbRAV.isSelected()))
                    && ((cbNVTE.isEnabled()
                    && cbNVTE.isSelected())
                    || (cbEQUI.isEnabled()
                    && cbEQUI.isSelected())
                    || (cbPROD.isEnabled()
                    && cbPROD.isSelected()))) {
                buttonDrawChart.setEnabled(true);
                if (cbAutoDraw.isSelected()) {
                    plot();
                }
            } else {
                buttonDrawChart.setEnabled(false);
            }
        }
    }//GEN-LAST:event_comboXAxisActionPerformed
    
    /**
     * Set the proper title along the x-axis.
     * 
     * @see Ms2chart#comboXAxisActionPerformed(java.awt.event.ActionEvent) 
     * 
     * @param Otitle Title keywords
     */
    void settitleXAxis(String Otitle)
    {
        if(Otitle.equals("PROC"))
            titleXAxis = "Number of Processes";
        else
        if(Otitle.equals("TIME[ps]"))
            titleXAxis = "Timestep [ps]";
        else
        if(Otitle.equals("NR"))
            titleXAxis = "Number of steps";
        else
        if(Otitle.equals("PRESS"))
            titleXAxis = "Reduced pressure";
         else
        if(Otitle.equals("DENSITY"))
            titleXAxis = "Reduced density";
        else
        if(Otitle.equals("EPOT"))
            titleXAxis = "Reduced potential energy";
        else
        if(Otitle.equals("ENTLP"))
            titleXAxis = "Reduced enthalpy";
         else
        if(Otitle.equals("TEMP"))
            titleXAxis = "Reduced temperature";
        else
        if(Otitle.equals("MUE"))
            titleXAxis = "Chemical potential";
        else
        if(Otitle.equals("FRACT"))
            titleXAxis = "Mole fraction";
         else
        if(Otitle.equals("DISP"))
            titleXAxis = "Displacement";
        else
        if(Otitle.equals("VW"))
            titleXAxis = "Partial Molar Volume";
         else
        if(Otitle.equals("D_i"))
            titleXAxis = "ACF Self-diff. coeff.";
        else
        if(Otitle.equals("VS"))
            titleXAxis = "ACF Shear-Viscosity";
        else
        if(Otitle.equals("VB"))
            titleXAxis = "ACF Bulk-Viscosity";
         else
        if(Otitle.equals("CO"))
            titleXAxis = "ACF Thermal Conductivity";
        else
        if(Otitle.equals("EP_Intra"))
            titleXAxis = "Reduced intramolecular potential energy";
        else
        if(Otitle.equals("EP_Bonds"))
            titleXAxis = "Reduced bond energy";
         else
        if(Otitle.equals("EP_Angles"))
            titleXAxis = "Reduced angle energy";
        else
        if(Otitle.equals("EP_Dihed"))
            titleXAxis = "Reduced torsion energy";
        else
        if(Otitle.equals("EP_14_15"))
            titleXAxis = "Reduced 1-4/1-5 energy";
         else
        if(Otitle.equals("EP_Inter"))
            titleXAxis = "Reduced intermolecular potential energy";
        else
        if(Otitle.equals("Vir_Intra"))
            titleXAxis = "Virial intramolecular energy";
        else
        if(Otitle.equals("Vir_Inter"))
            titleXAxis = "Virial intermolecular energy";
        else
            titleXAxis = Otitle;
    }
    
    /**
     * Y-Axis control event handler
     * Checks different user input condition
     * set the title for x-axis
     * Prepare the arrays for calculations
     * and calls plot function accordingly
     * 
     * @see Ms2chart#plot()
     * @see Ms2chart#settitleXAxis(java.lang.String)
     * 
     * @param evt 
     */
    private void listYAxisesValueChanged(javax.swing.event.ListSelectionEvent evt) {//GEN-FIRST:event_listYAxisesValueChanged
        if (evt.getValueIsAdjusting() == false) {
            Object[] selectedValues = listYAxises.getSelectedValues();
            columnsInRtrFileForYAxis = new int[0];
            columnsInRunFileForYAxis = new int[0];
            columnsInRavFileForYAxis = new int[0];
            for (int i = 0; i < selectedValues.length; ++i) {
                String titleSelected = (String) selectedValues[i];
                if (rbRTR.isSelected() || rbTHI.isSelected() || rbDCP.isSelected()) { ///// I added rbTHI check with rtr
                    for (int j = 0; j < arrayOfRtrFile[0].length; ++j) {
                        if (titleSelected.equals(arrayOfRtrFile[0][j])) {
                            int[] tempArray = new int[columnsInRtrFileForYAxis.length + 1];
                            System.arraycopy(columnsInRtrFileForYAxis, 0, tempArray, 0, columnsInRtrFileForYAxis.length);
                            tempArray[columnsInRtrFileForYAxis.length] = j;
                            columnsInRtrFileForYAxis = tempArray;
                        }
                    }
                }
                if (cbRUN.isEnabled()) {
                    for (int j = 0; j < arrayOfRunFile[0].length; ++j) {
                        if (titleSelected.equals(arrayOfRunFile[0][j])) {
                            int[] tempArray = new int[columnsInRunFileForYAxis.length + 1];
                            System.arraycopy(columnsInRunFileForYAxis, 0, tempArray, 0, columnsInRunFileForYAxis.length);
                            tempArray[columnsInRunFileForYAxis.length] = j;
                            columnsInRunFileForYAxis = tempArray;
                        }
                    }
                }
                if (cbRAV.isEnabled()) {
                    for (int j = 0; j < arrayOfRavFile[0].length; ++j) {
                        if (titleSelected.equals(arrayOfRavFile[0][j])) {
                            int[] tempArray = new int[columnsInRavFileForYAxis.length + 1];
                            System.arraycopy(columnsInRavFileForYAxis, 0, tempArray, 0, columnsInRavFileForYAxis.length);
                            tempArray[columnsInRavFileForYAxis.length] = j;
                            columnsInRavFileForYAxis = tempArray;
                        }
                    }
                }
            }
        }
        if (columnInRunFileForXAxis == -1) {
            cbRUN.setEnabled(false);
        }
        if (columnInRavFileForXAxis == -1) {
            cbRAV.setEnabled(false);
        }
        if (rbRTR.isSelected() || rbTHI.isSelected() || rbDCP.isSelected()) { //// I added rbTHI check with rtr
            if (columnInRtrFileForXAxis != -1
                    && columnsInRtrFileForYAxis.length != 0) {
                buttonDrawChart.setEnabled(true);
                if (cbAutoDraw.isSelected()) {
                    plot();
                }
            } else {
                buttonDrawChart.setEnabled(false);
            }
        } else {
            if ((columnInRunFileForXAxis != -1
                    || columnInRavFileForXAxis != -1)
                    && (columnsInRunFileForYAxis.length != 0
                    || columnsInRavFileForYAxis.length != 0)
                    && ((cbRUN.isEnabled()
                    && cbRUN.isSelected())
                    || (cbRAV.isEnabled()
                    && cbRAV.isSelected()))
                    && ((cbNVTE.isEnabled()
                    && cbNVTE.isSelected())
                    || (cbEQUI.isEnabled()
                    && cbEQUI.isSelected())
                    || (cbPROD.isEnabled()
                    && cbPROD.isSelected()))) {
                buttonDrawChart.setEnabled(true);
                if (cbAutoDraw.isSelected()) {
                    plot();
                }
            } else {
                buttonDrawChart.setEnabled(false);
            }
        }
    }//GEN-LAST:event_listYAxisesValueChanged
    
    /**
     * Prepare data to show in the Table Frame.
     * @see ShowFileFrame
     * @param evt 
     */
    private void buttonViewDataActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_buttonViewDataActionPerformed
        XYDataset xydataset = jPanelOfChart.getChart().getXYPlot().getDataset();
        
        int[] arrayOfItemCounts = new int[xydataset.getSeriesCount()];
        for (int i = 0; i < xydataset.getSeriesCount(); ++i) {
            arrayOfItemCounts[i] = xydataset.getItemCount(i);
          
        }
        String[][] data = new String[arrayOfItemCounts[findIndexOfLargestInt(arrayOfItemCounts)]][xydataset.getSeriesCount() + 1];
       //JOptionPane.showMessageDialog(null,data.length);
        for (int series = 0; series < data[data.length - 1].length; ++series) {
            for (int item = 0; item < data.length; ++item) {
                if (series == 0) {
                   data[item][series] = xydataset.getXValue(findIndexOfLargestInt(arrayOfItemCounts) , item) + "";
                    // JOptionPane.showMessageDialog(null,xydataset.getXValue(findIndexOfLargestInt(arrayOfItemCounts) , item) + "");
                } else {
                    if (xydataset.getItemCount(series - 1) > item) {
                        if(acfIsNeeded)
                        {
                        data[item][series] = arrayOfRtrFile[arrayOfRtrFileStartingIndices[arrayOfRtrFileStartingIndices.length - 1] + 1 + item][columnsInRtrFileForYAxis[0]] + "";
                        
                        }
                        else
                        {
                        data[item][series] = xydataset.getYValue(series - 1, item) + "";
                        
                        }
                        
                        //JOptionPane.showMessageDialog(null,xydataset.getYValue(series - 1, item) + "");
                    } else {
                        data[item][series] = null;
                    }
                }
            }
        }
        String[] yAxisTitle = jPanelOfChart.getChart().getXYPlot().getRangeAxis().getLabel().toString().split(", ");
        String[] title = new String[yAxisTitle.length + 1];
        System.arraycopy(yAxisTitle, 0, title, 1, yAxisTitle.length);
        title[0] = titleXAxis;
        ShowTableFrame frame;
        frame = new ShowTableFrame(data, title);
        frame.setVisible(true);
    }//GEN-LAST:event_buttonViewDataActionPerformed
    
    private void cbAutoDrawItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_cbAutoDrawItemStateChanged
        setLabelStatusMessage("Toggles between use of Draw Button");
    }//GEN-LAST:event_cbAutoDrawItemStateChanged
    
    private void cbAskDirectoryItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_cbAskDirectoryItemStateChanged
        setLabelStatusMessage("Press Save Image to save in PNG format");
    }//GEN-LAST:event_cbAskDirectoryItemStateChanged
    
    /**
     * Prepare and open help file.
     * 
     * @param evt 
     */
    private void HelpbtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_HelpbtnActionPerformed
         
        
        try{
                String filePath = UserDirectory+""+File.separator+"help"+File.separator+"ms2chart_help_format.xlsx";
                Desktop.getDesktop().open(new File(filePath));
    }catch (Exception e){
      JOptionPane.showMessageDialog(null,"Exception caught ="+e.getMessage());
    }
    }//GEN-LAST:event_HelpbtnActionPerformed
    
    /**
     * Reset all the UI controls
     * 
     * @param evt 
     */
    private void buttonClearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_buttonClearActionPerformed
        buttonViewData.setEnabled(false);
        buttonPNGsave.setEnabled(false);
        buttonXAxisChange.setEnabled(false);
        buttonYAxisChange.setEnabled(false);
        buttonClear.setEnabled(false);
        panelChartContainer.removeAll();
        panelChartContainer.setVisible(false);
        panelChartContainer.setVisible(true);
    }//GEN-LAST:event_buttonClearActionPerformed
    
    /**
     * Reset global variable and restart the program
     * 
     * @param evt 
     */
    private void buttonRestartActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_buttonRestartActionPerformed
        setVisible(false);
        dispose();
        columnInRunFileForXAxis = -1;
        columnInRavFileForXAxis = -1;
        columnsInRunFileForYAxis = new int[0];
        columnsInRavFileForYAxis = new int[0];
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new Ms2chart().setVisible(true);
            }
        });
    }//GEN-LAST:event_buttonRestartActionPerformed
    
    /**
     * Trim the name of the file.
     * 
     * @return Trimmed name of the file
     */
    private String TrimCasetitle()
    {
        String TrimTitle = caseTitle.substring(0, caseTitle.length() - 2);
        return TrimTitle;
    }
    
    /**
     * Open the log file.
     * 
     * @param evt 
     */
    private void buttonLOGActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_buttonLOGActionPerformed
        ShowFileFrame frame;
        if (isHavingCases) {
            if(rbRTR.isSelected())
            {frame = new ShowFileFrame(selectedDirectory.getAbsolutePath() + java.io.File.separator + TrimCasetitle() + ".log");}
            else
            {frame = new ShowFileFrame(selectedDirectory.getAbsolutePath() + java.io.File.separator + caseTitle + ".log");}
        } else {
            frame = new ShowFileFrame();
        }
        frame.setVisible(true);
    }//GEN-LAST:event_buttonLOGActionPerformed
    
    /**
     * Prepare and Open PAR file.
     * 
     * @param evt 
     */
    private void buttonPARActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_buttonPARActionPerformed
        ShowFileFrame frame;
        if (isHavingCases) {
            if(rbRTR.isSelected())
            {frame = new ShowFileFrame(selectedDirectory.getAbsolutePath() + java.io.File.separator + TrimCasetitle() + ".par");}
            else
            {frame = new ShowFileFrame(selectedDirectory.getAbsolutePath() + java.io.File.separator + caseTitle + ".par");}
        } else {
            frame = new ShowFileFrame();
        }
        frame.setVisible(true);
    }//GEN-LAST:event_buttonPARActionPerformed
    
    /**
     * Prepare and Open Res file.
     * 
     * @param evt 
     */
    private void buttonRESActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_buttonRESActionPerformed
        ShowFileFrame frame;
        if (isHavingCases) {
            if (rbRTR.isSelected()) {
                frame = new ShowFileFrame(selectedDirectory.getAbsolutePath() + java.io.File.separator + caseTitle + ".res");
            } else if (rbRUN_RAV.isSelected()) {
                frame = new ShowFileFrame(selectedDirectory.getAbsolutePath() + java.io.File.separator + caseTitle + "_1.res");
            } else {
                frame = new ShowFileFrame();
            }
        } else {
            frame = new ShowFileFrame();
        }
        frame.setVisible(true);
    }//GEN-LAST:event_buttonRESActionPerformed
    
    /**
     * Save the png image of Graph chart in the selected user directory.
     * 
     * @param evt 
     */
    private void buttonPNGsaveActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_buttonPNGsaveActionPerformed
        try {
            JFreeChart chart = jPanelOfChart.getChart();
            Object[] yTitles = listYAxises.getSelectedValues();
            String yName = "";
            for (int i = 0; i < yTitles.length; ++i) {
                yName += ((String) yTitles[i]).replaceAll("\\s", " ") + "_";
            }
            String xName = ((String) comboXAxis.getSelectedItem()).replaceAll("\\s", " ");
            File file = new File(selectedDirectory.getAbsolutePath() + java.io.File.separator + caseTitle + "_" + yName + "_" + xName + ".png");
            if (cbAskDirectory.isSelected()) {
                JFileChooser fileChooser = new JFileChooser();
                fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
                if (fileChooser.showSaveDialog(this) == JFileChooser.APPROVE_OPTION) {
                    file = new File(fileChooser.getSelectedFile().getAbsolutePath() + java.io.File.separator + caseTitle + "_" + yName + "_" + xName + ".png");
                }
            }
            try {
                if (file.createNewFile()) {
                    ImageIO.write(chart.createBufferedImage(panelChartContainer.getWidth(), panelChartContainer.getHeight()), "png", file);
                } else {
                    setLabelStatusMessage("File already exists");
                }
            } catch (Exception e) {
                setLabelStatusMessage("ERROR!");
            }
        } catch (Exception e) {
            setLabelStatusMessage("Error reading chart!");
        }
    }//GEN-LAST:event_buttonPNGsaveActionPerformed
    
    private void buttonDrawChartActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_buttonDrawChartActionPerformed
        plot();
    }//GEN-LAST:event_buttonDrawChartActionPerformed
    
    /**
     * Set minimum and maximum value for the y-axis
     * 
     * @param evt 
     */
    private void buttonYAxisChangeActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_buttonYAxisChangeActionPerformed
        if(this.readDouble(tfYMax.getText()) && this.readDouble(tfYMin.getText())) {
            rangeAxis.setUpperBound(Double.parseDouble(tfYMax.getText()));
            rangeAxis.setLowerBound(Double.parseDouble(tfYMin.getText()));
            tfYMax.setText(Math.round(rangeAxis.getRange().getUpperBound() * 100) / 100.0 + "");
            tfYMin.setText(Math.round(rangeAxis.getRange().getLowerBound() * 100) / 100.0 + "");
        }
        else{
            tfYMax.setText(Math.round(rangeAxis.getRange().getUpperBound() * 100) / 100.0 + "");
            tfYMin.setText(Math.round(rangeAxis.getRange().getLowerBound() * 100) / 100.0 + "");
        }
    }//GEN-LAST:event_buttonYAxisChangeActionPerformed
    
    /**
     * Set minimum and maximum value for the x-axis
     * 
     * @param evt 
     */
    private void buttonXAxisChangeActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_buttonXAxisChangeActionPerformed
        if(this.readDouble(tfXMax.getText()) && this.readDouble(tfXMin.getText())) {
            domainAxis.setUpperBound(Double.parseDouble(tfXMax.getText()));
            domainAxis.setLowerBound(Double.parseDouble(tfXMin.getText()));
            tfXMax.setText(Math.round(domainAxis.getRange().getUpperBound() * 100) / 100.0 + "");
            tfXMin.setText(Math.round(domainAxis.getRange().getLowerBound() * 100) / 100.0 + "");
        }
        else{
            tfXMax.setText(Math.round(domainAxis.getRange().getUpperBound() * 100) / 100.0 + "");
            tfXMin.setText(Math.round(domainAxis.getRange().getLowerBound() * 100) / 100.0 + "");
        }
    }//GEN-LAST:event_buttonXAxisChangeActionPerformed
    
    private boolean readDouble(String str){
        double ret = 0.0;
        try {
            ret = Double.parseDouble(str);
        } catch(Exception e) {
            JOptionPane.showMessageDialog(new JFrame(),"Error: "+str+" is not a decimal number!");
            return false;
        }
        
        return true;
    }
    
    /**
     * Show open dialog,
     * Set user Directory,
     * Populate combo box,
     * Enable and disable UI accordingly
     * 
     * @param evt 
     */
    private void buttonBrowseActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_buttonBrowseActionPerformed
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
        fileChooser.setCurrentDirectory(selectedDirectory);
        FileFilter filter_par = new FileNameExtensionFilter("Complete Simulation","res","vim","log","rav","rtr","run","rst","par","txt","thi","dcp");
        fileChooser.setFileFilter(filter_par);
        FileFilter filter_vim = new FileNameExtensionFilter("Single vim file","vim");
        fileChooser.addChoosableFileFilter(filter_vim);
        fileChooser.setAcceptAllFileFilterUsed(false);
        
        if (fileChooser.showOpenDialog(this) == JFileChooser.APPROVE_OPTION) {
            if(fileChooser.getSelectedFile().isFile())
            {
                if(fileChooser.getSelectedFile().getName().endsWith(".vim"))
                {
                    String[] VimCaseTitle = {fileChooser.getSelectedFile().getName().replaceFirst("[_][^_]+$", "")};
                    comboCaseList.setModel(new javax.swing.DefaultComboBoxModel(VimCaseTitle));
                    selectedDirectory = fileChooser.getSelectedFile().getParentFile();
                }
                else
                    JOptionPane.showMessageDialog(null,"File cannot be selected. Please select directory"); 
                
            }
            else
            {
            selectedDirectory = fileChooser.getSelectedFile();
            setLabelStatusMessage("In " + selectedDirectory.getAbsolutePath());
            comboCaseList.setModel(new javax.swing.DefaultComboBoxModel(getDirList()));
            }
            if (isHavingCases) {
                cbEQUI.setSelected(false);
                cbEQUI.setEnabled(false);
                prepareData();
            } else {
                cbRUN.setEnabled(false);
                cbRAV.setEnabled(false);
                cbNVTE.setEnabled(false);
                cbEQUI.setEnabled(false);
                cbPROD.setEnabled(false);
                tfXMax.setEnabled(false);
                tfXMin.setEnabled(false);
                tfYMax.setEnabled(false);
                tfYMin.setEnabled(false);
                tfError.setEnabled(false);
                buttonLOG.setEnabled(false);
                buttonPAR.setEnabled(false);
                buttonRES.setEnabled(false);
                buttonDrawChart.setEnabled(false);
                buttonViewData.setEnabled(false);
                buttonPNGsave.setEnabled(false);
                buttonXAxisChange.setEnabled(false);
                buttonYAxisChange.setEnabled(false);
                buttonClear.setEnabled(false);
                comboXAxis.setModel(new javax.swing.DefaultComboBoxModel(new String[] {""}));
                listYAxises.setModel(new javax.swing.AbstractListModel() {
                    String[] strings = new String[] {""};
                    public int getSize() { return strings.length; }
                    public Object getElementAt(int i) { return strings[i]; }
                });
                setLabelStatusMessage("No .rav or .run file found. Browse again...");
                panelChartContainer.removeAll();
                panelChartContainer.setVisible(false);
                panelChartContainer.setVisible(true);
            }
        /// outer else bracket}
        }
    }//GEN-LAST:event_buttonBrowseActionPerformed
    
    /**
     * RAV value changed event handler
     * Checks condition for RAV.
     * 
     * @see Ms2chart#plot()
     * 
    */
    private void cbRAVItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_cbRAVItemStateChanged
        if ((columnInRunFileForXAxis != -1
                || columnInRavFileForXAxis != -1)
                && (columnsInRunFileForYAxis.length != 0
                || columnsInRavFileForYAxis.length != 0)
                && ((cbRUN.isEnabled()
                && cbRUN.isSelected())
                || (cbRAV.isEnabled()
                && cbRAV.isSelected()))
                && ((cbNVTE.isEnabled()
                && cbNVTE.isSelected())
                || (cbEQUI.isEnabled()
                && cbEQUI.isSelected())
                || (cbPROD.isEnabled()
                && cbPROD.isSelected()))) {
            buttonDrawChart.setEnabled(true);
            if (cbAutoDraw.isSelected()) {
                plot();
            }
        } else {
            buttonDrawChart.setEnabled(false);
        }
    }//GEN-LAST:event_cbRAVItemStateChanged
    
    /**
     * RUN value changed event handler
     * Checks condition for RAV.
     * 
     * @see Ms2chart#plot()
     * 
    */
    private void cbRUNItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_cbRUNItemStateChanged
        if ((columnInRunFileForXAxis != -1
                || columnInRavFileForXAxis != -1)
                && (columnsInRunFileForYAxis.length != 0
                || columnsInRavFileForYAxis.length != 0)
                && ((cbRUN.isEnabled()
                && cbRUN.isSelected())
                || (cbRAV.isEnabled()
                && cbRAV.isSelected()))
                && ((cbNVTE.isEnabled()
                && cbNVTE.isSelected())
                || (cbEQUI.isEnabled()
                && cbEQUI.isSelected())
                || (cbPROD.isEnabled()
                && cbPROD.isSelected()))) {
            buttonDrawChart.setEnabled(true);
            if (cbAutoDraw.isSelected()) {
                plot();
            }
        } else {
            buttonDrawChart.setEnabled(false);
        }
    }//GEN-LAST:event_cbRUNItemStateChanged
    
    /**
     * Prepare vim file for read operation
     * 
     * @see Ms2chart#Readvimfile(java.lang.String, java.lang.String)
     * 
     * @param evt 
     */
    private void vimtovmd_btnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_vimtovmd_btnActionPerformed
        String casetitle = (String) comboCaseList.getSelectedItem();
        String file = casetitle + "_1.vim";
        boolean check = new File(selectedDirectory.getAbsolutePath(),file).exists();
        if(!check)
        {
            setLabelStatusMessage("Require "+file+" not found in the specified directory");
            JOptionPane.showMessageDialog(null,"Required "+file+" not found in the specified directory");
        }
        else
        {
            this.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
            this.setEnabled(false);
            if(Readvimfile(file,casetitle))
            {
            JOptionPane.showMessageDialog(null, "VMD file have been created successfuly");
            }
            this.setEnabled(true);
            this.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
            
        }
    }//GEN-LAST:event_vimtovmd_btnActionPerformed

    private void buttonQuit1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_buttonQuit1ActionPerformed
        // TODO add your handling code here:
        System.exit(0);
        
    }//GEN-LAST:event_buttonQuit1ActionPerformed
    /**
     * THI radio button event handler
     * 
     * @param evt 
     */
    private void rbTHIActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_rbTHIActionPerformed
        panelChartContainer.removeAll();
        setLabelStatusMessage("Listing .thi files");
        comboCaseList.setModel(new javax.swing.DefaultComboBoxModel(getDirList()));
        if (isHavingCases) {
            prepareData();
        } else {
            cbRUN.setEnabled(false);
            cbRAV.setEnabled(false);
            cbNVTE.setEnabled(false);
            cbEQUI.setEnabled(false);
            cbPROD.setEnabled(false);
            tfXMax.setEnabled(false);
            tfXMin.setEnabled(false);
            tfYMax.setEnabled(false);
            tfYMin.setEnabled(false);
            tfError.setEnabled(false);
            buttonLOG.setEnabled(false);
            buttonPAR.setEnabled(false);
            buttonRES.setEnabled(false);
            buttonDrawChart.setEnabled(false);
            buttonViewData.setEnabled(false);
            buttonPNGsave.setEnabled(false);
            buttonXAxisChange.setEnabled(false);
            buttonYAxisChange.setEnabled(false);
            buttonClear.setEnabled(false);
            comboXAxis.setModel(new javax.swing.DefaultComboBoxModel(new String[] {""}));
            listYAxises.setModel(new javax.swing.AbstractListModel() {
                String[] strings = new String[] {""};
                public int getSize() { return strings.length; }
                public Object getElementAt(int i) { return strings[i]; }
            });
            setLabelStatusMessage("No file found. Browse again...");
            panelChartContainer.removeAll();
            panelChartContainer.setVisible(false);
            panelChartContainer.setVisible(true);
        }
    }//GEN-LAST:event_rbTHIActionPerformed
    
    /**
     * DCP radio button event handler
     * 
     * @param evt 
     */
    private void rbDCPActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_rbDCPActionPerformed
       panelChartContainer.removeAll();
        setLabelStatusMessage("Listing .dcp files");
        comboCaseList.setModel(new javax.swing.DefaultComboBoxModel(getDirList()));
        if (isHavingCases) {
            prepareData();
        } else {
            cbRUN.setEnabled(false);
            cbRAV.setEnabled(false);
            cbNVTE.setEnabled(false);
            cbEQUI.setEnabled(false);
            cbPROD.setEnabled(false);
            tfXMax.setEnabled(false);
            tfXMin.setEnabled(false);
            tfYMax.setEnabled(false);
            tfYMin.setEnabled(false);
            tfError.setEnabled(false);
            buttonLOG.setEnabled(false);
            buttonPAR.setEnabled(false);
            buttonRES.setEnabled(false);
            buttonDrawChart.setEnabled(false);
            buttonViewData.setEnabled(false);
            buttonPNGsave.setEnabled(false);
            buttonXAxisChange.setEnabled(false);
            buttonYAxisChange.setEnabled(false);
            buttonClear.setEnabled(false);
            comboXAxis.setModel(new javax.swing.DefaultComboBoxModel(new String[] {""}));
            listYAxises.setModel(new javax.swing.AbstractListModel() {
                String[] strings = new String[] {""};
                public int getSize() { return strings.length; }
                public Object getElementAt(int i) { return strings[i]; }
            });
            setLabelStatusMessage("No file found. Browse again...");
            panelChartContainer.removeAll();
            panelChartContainer.setVisible(false);
            panelChartContainer.setVisible(true);
        }
    }//GEN-LAST:event_rbDCPActionPerformed
    
    /**
     * Read vim file frame vise and write frame data in a seperate files (.xyz)
     * in a specific folder.
     * 
     * @see Ms2chart#calcTowrite(java.util.ArrayList, double, int, int)
     * @see Ms2chart#writexyzfiles(java.lang.StringBuffer, java.lang.String)
     * 
     * @param fileName Name of the file
     * @param casetitle Case which is opened by the user
     * 
     * @return True - False
     */
    private boolean Readvimfile(String fileName, String casetitle)
    {
        BufferedReader inputStream;
        String line;
        int rows = 0,columns = 0, frnum = 0, frmrows = 0;
        double vol = 0;
        ArrayList<String> frmdata = new ArrayList<String>();
        StringBuffer linedata = new StringBuffer();
        StringBuffer completedata = new StringBuffer();
        
        File inputFile = new File(selectedDirectory.getAbsolutePath() + java.io.File.separator + fileName);
        try {
            inputStream = new BufferedReader(new FileReader(inputFile));
            //String dirname = makedir(casetitle);
                while (inputStream.ready()) {
                    line = inputStream.readLine();
                    String[] linearray = line.split("\\s+");
                    
                    if(linearray[0].equals("~"))
                    {
                        if(columns == 0)
                        {
                        columns = linearray.length - 1;
                        }
                        
                        rows++;
                        linedata.append(line + " ");
                    }
                    else
                        if(linearray[0].equals("!"))
                        {
                            //framedata.append(line + " ");
                            //JOptionPane.showMessageDialog(null, frnum+" "+line);
                            frmdata.add(frmrows,line);
                            frmrows++;
                        }
                    else
                       if(linearray[0].equals("#"))
                       {
                           if(frnum == 0)
                           {
                               calcTypes(linedata,rows);
                               vol =  Double.parseDouble(linearray[1]);
                               frnum++;
                           }
                           else
                           {
                              completedata.append(calcTowrite(frmdata,vol,frmrows,frnum));
                              //writexyzfiles(dirname,calcTowrite(frmdata,vol,frmrows),frnum,casetitle);
                              frmdata.clear();
                              vol =  Double.parseDouble(linearray[1]);
                              frnum++;
                              frmrows = 0;
                           }                               
                       }
                }
                   // else
                   // if(linearray[0].equals("##") || line.equals("##"))// EOF ....
                   // {
                              completedata.append(calcTowrite(frmdata,vol,frmrows,frnum));
                              frmdata.clear();
                   //     break;
                   // }
            //}
                
                writexyzfiles(completedata,casetitle);
        }
        catch(Exception e)
        {
            JOptionPane.showMessageDialog(null, "An error occured while processing... Please try again");
            setLabelStatusMessage("An error occured while processing... Please try again");
            e.printStackTrace();
            return false;
            //setLabelStatusMessage("error while processing "+fileName+" file. Please try again" + e.toString());
        }
        
        return true;
    }
    
    int[] nsites;
    double[][] Xsort,Ysort,Zsort,Colsort;
    /**
     * Prepare arrays of Type, X , Y, Z, Sigma and Columns from the data
     * in the start of the .vim file
     * 
     * @see Ms2chart#Readvimfile(java.lang.String, java.lang.String)
     * 
     * @param lineData Data from start of the frame
     * @param rows Number of rows
     */
    private void calcTypes(StringBuffer lineData,int rows)
    {
        StringTokenizer st = new StringTokenizer(lineData.toString());
        String Type[] = new String[rows];
        String X[] = new String[rows];
        String Y[] = new String[rows];
        String Z[] = new String[rows];
        String Sigma[] = new String[rows];
        String Col[] = new String[rows];
        for (int i = 0; i < rows; ++i) {
                if (st.hasMoreTokens()) {
                    st.nextToken();
                    Type[i] = st.nextToken();
                    st.nextToken();
                    X[i] = st.nextToken();
                    Y[i] = st.nextToken();
                    Z[i] = st.nextToken();
                    Sigma[i] = st.nextToken();
                    Col[i] = st.nextToken();
                } else {
                    break;
                }
        }
        getnumofsites(Type,X,Y,Z,Col);
    }
    
    /**
     * Calculates the first character for every line in .xyz file
     * from the strings of predefined characters.
     * 
     * @see Ms2chart#Readvimfile(java.lang.String, java.lang.String)
     * 
     * @param frmdata Frame data
     * @param vol     Volume of a frame. (In the start of the frame)
     * @param frmrows Number of rows
     * @param frmnum  Number of the current frame
     * 
     * 
     * @return StringBuffer (Data prepared to be written)
     */
    private StringBuffer calcTowrite(ArrayList<String> frmdata, double vol, int frmrows, int frmnum)
    {
        double[] writepos = {0,0,0};
        StringBuffer wData = new StringBuffer();
        wData.append(String.format("%d\n VIEW    1%03d\n",frmrows,frmnum));
        String str="HZOBSPPPCNNCOOH";
                                for(String put:frmdata)
                                {
                                   String[] getd = put.split("\\s+");
                                   double a = Math.round(((Double.parseDouble(getd[2])/1000) * vol) * 10000.0)/10000.0;
                                   double b = Math.round(((Double.parseDouble(getd[3])/1000) * vol) * 10000.0)/10000.0;
                                   double c = Math.round(((Double.parseDouble(getd[4])/1000) * vol) * 10000.0)/10000.0;
                                   double[] pos = {a,b,c};
                                   double[][] M = calnorm(getd);
                                   int type = Integer.parseInt(getd[1])-1;
                                   for(int k=0;k<nsites[type]; k++)
                                   {
                                    double[] sitelocal = {Xsort[type][k],Ysort[type][k],Zsort[type][k]};
                                    for(int i=0;i<sitelocal.length;i++)
                                    {
                                        for(int j=0;j<sitelocal.length;j++)
                                        {
                                            writepos[i] += M[i][j] * sitelocal[j];
                                        }
                                    }
                                    for(int l=0; l<pos.length;l++)
                                    {
                                        writepos[l] += pos[l];
                                    }
                                    wData.append(str.charAt((int)Colsort[type][k])+"    "+writepos[0]+"     "+writepos[1]+"     "
                                            +writepos[2]+"\n");
                                    writepos[0] = writepos[1] = writepos[2] = 0;
                                   }
                                }
             
        return wData;
    }
    
    /**
     * Get the number of sites into the global arrays and use it later.
     * 
     * @see Ms2chart#calcTypes(java.lang.StringBuffer, int)
     * 
     * @param type
     * @param X
     * @param Y
     * @param Z
     * @param Col 
     */
    private void getnumofsites(String[] type, String[] X, String[] Y,String[] Z,String[] Col)
    {
        int[] gettypes = calntypes(type);
        nsites = new int[gettypes[0]];
        Xsort = new double[type.length][gettypes[1]];
        Ysort = new double[type.length][gettypes[1]];
        Zsort = new double[type.length][gettypes[1]];
        Colsort = new double[type.length][gettypes[1]];
        int temp = 1,count = -1,index = 0;
        //JOptionPane.showMessageDialog(null, "columns: "+gettypes[1]+" rows: "+type.length+" diffT: "+gettypes[0]);
        for(int i=0;i<type.length;i++)
        {
            int num = Integer.parseInt(type[i]);
            if(num==temp)
            {
                nsites[index] = nsites[index] + 1;
                count = count + 1;
            }
            else
            {
                temp = temp + 1;
                index++;
                nsites[index] = nsites[index] + 1;
                count = 0;
            }
            //JOptionPane.showMessageDialog(null, "temp: "+temp+" count: "+count+" i: "+i);
            Xsort[index][count] = Double.parseDouble(X[i]);
            Ysort[index][count] = Double.parseDouble(Y[i]);
            Zsort[index][count] = Double.parseDouble(Z[i]);
            Colsort[index][count] = Double.parseDouble(Col[i]);
        }
    }
    
    /**
     * Calculate the N number of types from the Types array.
     * 
     * @see Ms2chart#getnumofsites(java.lang.String[], java.lang.String[], java.lang.String[], java.lang.String[], java.lang.String[])
     * 
     * @param type Array of types
     * 
     * @return Array of int
     */
    private int[] calntypes(String[] type)
    {
        int[] Ntypes = {1,1};
        int prenum = 1,precount = 0;
        for(int i=0;i<type.length;i++)
        {
            int num = Integer.parseInt(type[i]);
            if(num > Ntypes[0])
            {
                Ntypes[0] = Ntypes[0] + 1;
                prenum = num; 
                precount = 1;
            }
            else
            if(prenum == num)
            {
                    precount = precount + 1;
            }
            
            
            if(Ntypes[1]<precount)
             {
                    Ntypes[1] = precount;
                    
             } 
        }
        
        return Ntypes;
    }
    
    /**
     * Some calculations for the norm.
     * 
     * @see Ms2chart#calcTowrite(java.util.ArrayList, double, int, int)
     * 
     * @param numb Array of line data
     * 
     * @return Array of norm
     */
    private double[][] calnorm(String[] numb)
    {
        double a = Double.parseDouble(numb[5]);
        double b = Double.parseDouble(numb[6]);
        double c = Double.parseDouble(numb[7]);
        double d = Double.parseDouble(numb[8]);
       
        double norm = Math.sqrt((a*a)+(b*b)+(c*c)+(d*d));
        
        a = a/norm;
        b = b/norm;
        c = c/norm;
        d = d/norm;
        
        double[][] M = new double[3][3];
        
        M[0][0] = 1-(2*((c*c)+(d*d))); M[0][1] = 2*((b*c)-(a*d)); M[0][2] = 2*((b*d)+(a*c));
        M[1][0] = 2*((b*c)+(a*d)); M[1][1] = 1-(2*((d*d)+(b*b))); M[1][2] = 2*((c*d)-(a*b));
        M[2][0] = 2*((b*d)-(a*c)); M[2][1] = 2*((c*d)+(a*b)); M[2][2] = 1-(2*((b*b)+(c*c)));
        
        return M;
    }
    
    /**
     * Write .xyz file into the specified directory
     * 
     * @see Ms2chart#Readvimfile(java.lang.String, java.lang.String) 
     * 
     * @param Wdata Data to be written
     * @param Casetitle Name of the file
     */
    private void writexyzfiles(StringBuffer Wdata,String Casetitle)
    {
        try {
          File file = new File(String.format(selectedDirectory.getAbsolutePath() + java.io.File.separator + Casetitle +".xyz"));
          BufferedWriter output = new BufferedWriter(new FileWriter(file));
          
          output.write(Wdata.toString());
          
          output.close();
        } catch (Exception e ) {
           e.printStackTrace();
        }
    }
    /**
     * main function for Ms2chart. Starts the program.
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        startingDirectory = new File(new File("").getAbsolutePath());
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new Ms2chart().setVisible(true);
            }
        });
    }
    
    /**
     * Setting the appropriate filter and getting the list of specified files
     * 
     * @return Array of specified list of files
     */
    private String[] getDirList() {
        if (rbRTR.isSelected()) {
            String[] dirListForRtrFile;
            FilenameFilter rtrFilter = new FilenameFilter() {
                public boolean accept(File dir, String name) {
                    return name.endsWith(".rtr");
                }
            };
            dirListForRtrFile = selectedDirectory.list(rtrFilter);
            if (dirListForRtrFile.length > 0) {
                for (int i = 0; i < dirListForRtrFile.length; ++i) {
                    dirListForRtrFile[i] = dirListForRtrFile[i].substring(0, dirListForRtrFile[i].length() - 4);
                }
            }
            Arrays.sort(dirListForRtrFile);
            if (dirListForRtrFile.length > 0) {
                isHavingCases = true;
                return dirListForRtrFile;
            } else {
                isHavingCases = false;
                return new String[] {""};
            }
        }else
            if (rbTHI.isSelected()) {
            String[] dirListForRtrFile;
            FilenameFilter rtrFilter = new FilenameFilter() {
                public boolean accept(File dir, String name) {
                    return name.endsWith(".thi");
                }
            };
            dirListForRtrFile = selectedDirectory.list(rtrFilter);
            if (dirListForRtrFile.length > 0) {
                for (int i = 0; i < dirListForRtrFile.length; ++i) {
                    dirListForRtrFile[i] = dirListForRtrFile[i].substring(0, dirListForRtrFile[i].length() - 6);
                }
            }
            Arrays.sort(dirListForRtrFile);
            if (dirListForRtrFile.length > 0) {
                isHavingCases = true;
                return dirListForRtrFile;
            } else {
                isHavingCases = false;
                return new String[] {""};
            }
        } else
            if (rbDCP.isSelected()) {
            String[] dirListForRtrFile;
            FilenameFilter rtrFilter = new FilenameFilter() {
                public boolean accept(File dir, String name) {
                    return name.endsWith(".dcp");
                }
            };
            dirListForRtrFile = selectedDirectory.list(rtrFilter);
            if (dirListForRtrFile.length > 0) {
                for (int i = 0; i < dirListForRtrFile.length; ++i) {
                    dirListForRtrFile[i] = dirListForRtrFile[i].substring(0, dirListForRtrFile[i].length() - 6);
                }
            }
            Arrays.sort(dirListForRtrFile);
            if (dirListForRtrFile.length > 0) {
                isHavingCases = true;
                return dirListForRtrFile;
            } else {
                isHavingCases = false;
                return new String[] {""};
            }
        }
        else {
            String[] dirListForRunFile, dirListForRavFile, tempDirList;
            FilenameFilter runFilter = new FilenameFilter() {
                public boolean accept(File dir, String name) {
                    return name.endsWith(".run");
                }
            };
            FilenameFilter ravFilter = new FilenameFilter() {
                public boolean accept(File dir, String name) {
                    return name.endsWith(".rav");
                }
            };
            dirListForRunFile = selectedDirectory.list(ravFilter);
            dirListForRavFile = selectedDirectory.list(runFilter);
            if (dirListForRunFile.length > 0) {
                for (int i = 0; i < dirListForRunFile.length; ++i) {
                    dirListForRunFile[i] = dirListForRunFile[i].substring(0, dirListForRunFile[i].length() - 6);
                }
            }
            if (dirListForRavFile.length > 0) {
                for (int i = 0; i < dirListForRavFile.length; ++i) {
                    dirListForRavFile[i] = dirListForRavFile[i].substring(0, dirListForRavFile[i].length() - 6);
                }
            }
            tempDirList = new String[dirListForRunFile.length + dirListForRavFile.length];
            System.arraycopy(dirListForRunFile, 0, tempDirList, 0, dirListForRunFile.length);
            System.arraycopy(dirListForRavFile, 0, tempDirList, dirListForRunFile.length, dirListForRavFile.length);
            Arrays.sort(tempDirList);
            if (tempDirList.length > 0) {
                isHavingCases = true;
                String[] dirList = new String[1];
                dirList[0] = tempDirList[0];
                for (int i = 1; i < tempDirList.length; ++i) {
                    if (!tempDirList[i].equals(tempDirList[i - 1])) {
                        String[] tempArray = new String[dirList.length + 1];
                        System.arraycopy(dirList, 0, tempArray, 0, dirList.length);
                        tempArray[dirList.length] = tempDirList[i];
                        dirList = tempArray;
                    }
                }
                return dirList;
            } else {
                isHavingCases = false;
                return new String[] {""};
            }
        }
    }
    
    /**
     * Open file according to the specified case and enable/disable specified controls
     * 
     * @see Ms2chart#readFileIn2DStringArray(java.lang.String)
     * @see Ms2chart#getArrayOfDatasetLengths(java.lang.String[][], int[])
     * @see Ms2chart#getArrayOfStartingIndices(java.lang.String[][])
     */
    private void prepareData() {
        FileInputStream fstream;
        caseTitle = (String) comboCaseList.getSelectedItem();
        if (rbRTR.isSelected()) {
            try {
                fstream = new FileInputStream(selectedDirectory.getAbsolutePath() + java.io.File.separator + caseTitle + ".rtr");
                fstream.close();
                arrayOfRtrFile = readFileIn2DStringArray(caseTitle + ".rtr");
                arrayOfRtrFileStartingIndices = getArrayOfStartingIndices(arrayOfRtrFile);
                arrayOfRtrFileDatasetLengths = getArrayOfDatasetLengths(arrayOfRtrFile, arrayOfRtrFileStartingIndices);
                cbRUN.setEnabled(false);
                cbRAV.setEnabled(false);
                cbNVTE.setEnabled(false);
                cbEQUI.setEnabled(false);
                cbPROD.setEnabled(false);
                
            } catch (Exception e) {
                setLabelStatusMessage("Error selecting RTR files");
            }
        } else if(rbTHI.isSelected()){
            try {
                fstream = new FileInputStream(selectedDirectory.getAbsolutePath() + java.io.File.separator + caseTitle + "_1.thi");
                fstream.close();
                arrayOfRtrFile = readFileIn2DStringArray(caseTitle + "_1.thi");
                arrayOfRtrFileStartingIndices = getArrayOfStartingIndices(arrayOfRtrFile);
                arrayOfRtrFileDatasetLengths = getArrayOfDatasetLengths(arrayOfRtrFile, arrayOfRtrFileStartingIndices);
                cbRUN.setEnabled(false);
                cbRAV.setEnabled(false);
                cbNVTE.setEnabled(false);
                cbEQUI.setEnabled(false);
                cbPROD.setEnabled(false);
                
            } catch (Exception e) {
                setLabelStatusMessage("Error selecting thi files");
            }
        }
        else if(rbDCP.isSelected()){
            try {
                fstream = new FileInputStream(selectedDirectory.getAbsolutePath() + java.io.File.separator + caseTitle + "_1.dcp");
                fstream.close();
                arrayOfRtrFile = readFileIn2DStringArray(caseTitle + "_1.dcp");
               arrayOfRtrFileStartingIndices = getArrayOfStartingIndices(arrayOfRtrFile);
                arrayOfRtrFileDatasetLengths = getArrayOfDatasetLengths(arrayOfRtrFile, arrayOfRtrFileStartingIndices);
                cbRUN.setEnabled(false);
                cbRAV.setEnabled(false);
                cbNVTE.setEnabled(false);
                cbEQUI.setEnabled(false);
                cbPROD.setEnabled(false);
                
            } catch (Exception e) {
                setLabelStatusMessage("Error selecting dcp files");
            }
        }
        else{
            cbPROD.setEnabled(true);
            try {
                fstream = new FileInputStream(selectedDirectory.getAbsolutePath() + java.io.File.separator + caseTitle + "_1.run");
                fstream.close();
                arrayOfRunFile = readFileIn2DStringArray(caseTitle + "_1.run");
                arrayOfRunFileStartingIndices = getArrayOfStartingIndices(arrayOfRunFile);
                arrayOfRunFileDatasetLengths = getArrayOfDatasetLengths(arrayOfRunFile, arrayOfRunFileStartingIndices);
                if (arrayOfRunFileDatasetLengths.length > 2) {
                    cbEQUI.setEnabled(true);
                }
                if (arrayOfRunFileDatasetLengths.length > 1) {
                    cbNVTE.setEnabled(true);
                }
                cbRUN.setEnabled(true);
            } catch (Exception e) {
                cbRUN.setEnabled(false);
            }
            try {
                fstream = new FileInputStream(selectedDirectory.getAbsolutePath() + java.io.File.separator + caseTitle + "_1.rav");
                fstream.close();
                arrayOfRavFile = readFileIn2DStringArray(caseTitle + "_1.rav");
                setRangeIndexForErrors();
                arrayOfRavFileStartingIndices = getArrayOfStartingIndices(arrayOfRavFile);
                arrayOfRavFileDatasetLengths = getArrayOfDatasetLengths(arrayOfRavFile, arrayOfRavFileStartingIndices);
                if (arrayOfRavFileDatasetLengths.length > 2) {
                    cbEQUI.setEnabled(true);
                }
                if (arrayOfRavFileDatasetLengths.length > 1) {
                    cbNVTE.setEnabled(true);
                }
                cbRAV.setEnabled(true);
            } catch (Exception e) {
                cbRAV.setEnabled(false);
            }
        }
        try { ///////////////////////////////// I added rbTHI check here with rtr
            if(rbRTR.isSelected() || rbTHI.isSelected() || rbDCP.isSelected()){fstream = new FileInputStream(selectedDirectory.getAbsolutePath() + java.io.File.separator + TrimCasetitle() + ".log");}
            else{fstream = new FileInputStream(selectedDirectory.getAbsolutePath() + java.io.File.separator + caseTitle + ".log");}
            fstream.close();
            buttonLOG.setEnabled(true);
        } catch (Exception e) {
            buttonLOG.setEnabled(false);
        }
        try { ///////////////////////////////// I added rbTHI check here with rtr
            if(rbRTR.isSelected() || rbTHI.isSelected() || rbDCP.isSelected()){fstream = new FileInputStream(selectedDirectory.getAbsolutePath() + java.io.File.separator + TrimCasetitle() + ".par");}
            else{fstream = new FileInputStream(selectedDirectory.getAbsolutePath() + java.io.File.separator + caseTitle + ".par");}
            fstream.close();
            buttonPAR.setEnabled(true);
        } catch (Exception e) {
            buttonPAR.setEnabled(false);
        }
        try {
            if (rbRTR.isSelected() || rbTHI.isSelected() || rbDCP.isSelected()) { ///////////////////////////////// I added rbTHI check here with rtr
                fstream = new FileInputStream(selectedDirectory.getAbsolutePath() + java.io.File.separator + caseTitle + ".res");
                fstream.close();
                buttonRES.setEnabled(true);
            } else if (rbRUN_RAV.isSelected()) {
                fstream = new FileInputStream(selectedDirectory.getAbsolutePath() + java.io.File.separator + caseTitle + "_1.res");
                fstream.close();
                buttonRES.setEnabled(true);
            } else {
                buttonRES.setEnabled(false);
            }
        } catch (Exception e) {
            buttonRES.setEnabled(false);
        }
        comboXAxis.setModel(new javax.swing.DefaultComboBoxModel(getListOfXAxises()));
        comboXAxis.setSelectedIndex(0);
        int[] temp = listYAxises.getSelectedIndices();
        listYAxises.setModel(new javax.swing.AbstractListModel() {
            String[] strings = getListOfYAxises();
            public int getSize() { return strings.length; }
            public Object getElementAt(int i) { return strings[i]; }
        });
        listYAxises.setSelectedIndices(temp);
    }
    
    /**
     * Read the specified file data in two dimensional string
     * according to the different cases.
     * 
     * @see Ms2chart#prepareData()
     * 
     * @param fileName Name of the selected file
     * 
     * @return Two Dimensional array of string
     */
    private String[][] readFileIn2DStringArray(String fileName) {
        int numberOfColumns = 0, numberOfRows = 0;
        BufferedReader inputStream;
        String line;
        ArrayList<String> raw = new ArrayList<String>();
        StringBuffer fileData = new StringBuffer();
        File inputFile = new File(selectedDirectory.getAbsolutePath() + java.io.File.separator + fileName);
        try {
            inputStream = new BufferedReader(new FileReader(inputFile));
            if (rbRTR.isSelected()) {
                while ((line = inputStream.readLine()) != null) {
                    if (line.length() >= 6) {
                        ++numberOfRows;
                        String[] para = line.split("\\s+");
                        // ListIterator listiterator = raw.listIterator();
                        if (numberOfRows == 1) {
                           // numberOfColumns = line.length() / 10 + 1;
                        
                        int i=0,pos=0;
                        while(i<para.length)
                        {
                            if(para[i].matches("[0-9]"))
                            {
                                //raw[i-1] = raw[i-1]+""+para[i];
                               raw.set(pos-1,raw.get(pos-1)+" "+para[i]);
                               i++;
                            }
                            else
                                if(para[i].equals("Int"))
                                {
                                   // raw[i] = para[i]+""+para[i+1];
                                    raw.add(para[i]+" "+para[i+1]);
                                    i+=2;pos++;
                                }
                            else
                                {
                                    raw.add(para[i]);
                                    i++;pos++;
                                }
                        }
                        //JOptionPane.showMessageDialog(null,"hehe");
                      for(String put:raw)
                        {
                            //JOptionPane.showMessageDialog(null,"Line : " +put);
                            fileData.append(put+" ");
                        }
                        }
                        else
                        {
                            for(int i=0; i<para.length; i++)
                            {
                                fileData.append(para[i]+" ");
                            }
                        }
                        
                        if(numberOfRows==2)
                        {
                            numberOfColumns = para.length -1;
                        }
                    }
                }
            } else {
                while ((line = inputStream.readLine()) != null) {
                    if (line.length() >= 6 && !line.contains("Component") && !line.contains("currentlambda") && 
                            !line.trim().startsWith("#")) {
                        ++numberOfRows;
                       // System.out.println(">>> "+numberOfRows+" : ");
                       String[] para=line.split("\\s+");
                       int length = para.length;
                       
                      
                       if(numberOfRows == 2)
                        numberOfColumns = length - 1;
                       
                       for(int i=0; i<length; i++)
                            {
                               // System.out.print(para[i]+"  ");
                               /* if(numberOfRows==2)
                                {
                                    numberOfColumns = length - 1;
                                }*/
                                if(para[i].equals("-Infinity") || para[i].equals("Infinity") || para[i].equals("*******"))
                                {
                                    para[i] = "0.00000";  
                                }
                          
                                fileData.append(para[i]+" ");
                            }
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        
        setLabelStatusMessage("R: " + numberOfRows + " C: " + numberOfColumns + " " + fileName);
        StringTokenizer st = new StringTokenizer(fileData.toString());
        String[][] data = new String[numberOfRows][numberOfColumns];
        for (int i = 0; i < numberOfRows; ++i) {
            for (int j = 0; j < numberOfColumns; ++j) {
                if (st.hasMoreTokens()) {
                    data[i][j] = st.nextToken();
                } else {
                    break;
                }
              if (data[i][j].equals("MUE")
                || data[i][j].equals("FRACT")
                || data[i][j].equals("VW")
                || data[i][j].equals("HM")
                || data[i][j].equals("Int")
                || data[i][j].equals("IntD_i")
                || data[i][j].endsWith("D_i")
                || data[i][j].equals("Number")){
                    data[i][j] = data[i][j] + " " + st.nextToken();
                    if (data[i][j].equals("Int D")) {
                        data[i][j] = data[i][j] + " " + st.nextToken();
                    }
                }
            }
        }
        return data;
    }
    
    /**
     *  Check and save the indexes of the splited blocks of data contained in a file
     *  
     * @see Ms2chart#prepareData()
     * 
     * @param arrayOfFile Array of file data
     * 
     * @return Array of indexes
     */
    private int[] getArrayOfStartingIndices(String[][] arrayOfFile) {
        int[] arrayOfFileStartingIndices = new int[0];
        for (int i = 0; i < arrayOfFile.length; ++i) {
            if (arrayOfFile[i][0].equals(arrayOfFile[0][0])) {
                System.out.println(i+" 0 = "+arrayOfFile[i][0]+" = "+arrayOfFile[0][0]);
                int[] tempArray = new int[arrayOfFileStartingIndices.length + 1];
                System.arraycopy(arrayOfFileStartingIndices, 0, tempArray, 0, arrayOfFileStartingIndices.length);
                tempArray[arrayOfFileStartingIndices.length] = i;
                arrayOfFileStartingIndices = tempArray;
            }
        }
        
        return arrayOfFileStartingIndices;
    }
    
    /**
     * Check and save the length of each data block contained in a file.
     * 
     * @see Ms2chart#prepareData()
     * 
     * @param arrayOfFile Array of file data
     * @param arrayOfFileStartingIndices Array of indexes
     * 
     * @return 
     */
    private int[] getArrayOfDatasetLengths(String[][] arrayOfFile, int[] arrayOfFileStartingIndices) {
        int[] arrayOfFileDatasetLengths = new int[arrayOfFileStartingIndices.length];
        for (int i = 0; i < arrayOfFileDatasetLengths.length - 1; ++i) {
            arrayOfFileDatasetLengths[i] = arrayOfFileStartingIndices[i + 1] - arrayOfFileStartingIndices[i] - 1;
        }
        arrayOfFileDatasetLengths[arrayOfFileDatasetLengths.length - 1] = arrayOfFile.length - arrayOfFileStartingIndices[arrayOfFileStartingIndices.length - 1] - 1;

        return arrayOfFileDatasetLengths;
    }
    
    /**
     * Retrieve all the titles for x-axis from the array of file data
     * 
     * @see Ms2chart#prepareData()
     * 
     * @return Array of titles
     */
    private String[] getListOfXAxises() {
        String[] arrayOfRunFileXAxises = null;
        String[] arrayOfRavFileXAxises = null;
        String[] arrayOfRtrFileXAxises = null;
        String[] listOfXAxises, tempXlist;
        if (cbRUN.isEnabled()) {
            arrayOfRunFileXAxises = new String[arrayOfRunFile[0].length];
            for (int i = 0; i < arrayOfRunFileXAxises.length; ++i) {
                arrayOfRunFileXAxises[i] = arrayOfRunFile[0][i];
            }
        }
        if (cbRAV.isEnabled()) {
            arrayOfRavFileXAxises = new String[arrayOfRavFile[0].length];
            for (int i = 0; i < arrayOfRavFileXAxises.length; ++i) {
                arrayOfRavFileXAxises[i] = arrayOfRavFile[0][i];
            }
        }
        if (cbRUN.isEnabled() && cbRAV.isEnabled()) {
            boolean isPresent;
            for (int i = 0; i < arrayOfRavFileXAxises.length; ++i) {
                isPresent = false;
                for (int j = 0; j < arrayOfRunFileXAxises.length; ++j) {
                    if (arrayOfRavFileXAxises[i].equals(arrayOfRunFileXAxises[j])) {
                        isPresent = true;
                        break;
                    }
                }
                if (isPresent == false) {
                    String[] tempArray = new String[arrayOfRunFileXAxises.length + 1];
                    System.arraycopy(arrayOfRunFileXAxises, 0, tempArray, 0, arrayOfRunFileXAxises.length);
                    tempArray[arrayOfRunFileXAxises.length] = arrayOfRavFileXAxises[i];
                    arrayOfRunFileXAxises = tempArray;
                }
            }
            listOfXAxises = arrayOfRunFileXAxises;
        } else if (cbRUN.isEnabled()) {
            listOfXAxises = arrayOfRunFileXAxises;
        } else if (cbRAV.isEnabled()) {
            listOfXAxises = arrayOfRavFileXAxises;
        } else if (rbRTR.isSelected() || rbTHI.isSelected() || rbDCP.isSelected()) { ///////////////////////////////// I added rbTHI check here with rtr
            arrayOfRtrFileXAxises = new String[arrayOfRtrFile[0].length];
            for (int i = 0; i < arrayOfRtrFileXAxises.length; ++i) {
                arrayOfRtrFileXAxises[i] = arrayOfRtrFile[0][i];
            }
            listOfXAxises = arrayOfRtrFileXAxises;
        }
        else {
            listOfXAxises = new String[] {""};
        }
        return listOfXAxises;
    }
    
    /**
     * Retrieve all the titles for y-axis from the array of file data
     * 
     * @see Ms2chart#prepareData()
     * 
     * @return Array of titles
     */
    private String[] getListOfYAxises() {
        String[] arrayOfRunFileYAxises = null;
        String[] arrayOfRavFileYAxises = null;
        String[] arrayOfRtrFileYAxises = null;
        String[] listOfYAxises, tempYList;
        if (cbRUN.isEnabled()) {
            arrayOfRunFileYAxises = new String[arrayOfRunFile[0].length];
            for (int i = 0; i < arrayOfRunFileYAxises.length; ++i) {
                arrayOfRunFileYAxises[i] = arrayOfRunFile[0][i];
            }
        }
        if (cbRAV.isEnabled()) {
            arrayOfRavFileYAxises = new String[arrayOfRavFile[0].length];
            for (int i = 0; i < arrayOfRavFileYAxises.length; ++i) {
                arrayOfRavFileYAxises[i] = arrayOfRavFile[0][i];
            }
        }
        if (cbRUN.isEnabled() && cbRAV.isEnabled()) {
            boolean isPresent;
            for (int i = 0; i < arrayOfRavFileYAxises.length; ++i) {
                isPresent = false;
                for (int j = 0; j < arrayOfRunFileYAxises.length; ++j) {
                    if (arrayOfRavFileYAxises[i].equals(arrayOfRunFileYAxises[j])) {
                        isPresent = true;
                        break;
                    }
                }
                if (isPresent == false) {
                    String[] tempArray = new String[arrayOfRunFileYAxises.length + 1];
                    System.arraycopy(arrayOfRunFileYAxises, 0, tempArray, 0, arrayOfRunFileYAxises.length);
                    tempArray[arrayOfRunFileYAxises.length] = arrayOfRavFileYAxises[i];
                    arrayOfRunFileYAxises = tempArray;
                }
            }
            listOfYAxises = arrayOfRunFileYAxises;
        } else if (cbRUN.isEnabled()) {
            listOfYAxises = arrayOfRunFileYAxises;
        } else if (cbRAV.isEnabled()) {
            listOfYAxises = arrayOfRavFileYAxises;
        } else if (rbRTR.isSelected() || rbTHI.isSelected()) { ///////////////////////////////// I added rbTHI check here with rtr
            arrayOfRtrFileYAxises = new String[arrayOfRtrFile[0].length];
            for (int i = 1; i < arrayOfRtrFileYAxises.length; ++i) {
                arrayOfRtrFileYAxises[i] = arrayOfRtrFile[0][i];
            }
            listOfYAxises = arrayOfRtrFileYAxises;
        }
        else if (rbDCP.isSelected()) { ///////////////////////////////// I added rbdcp check here with rtr
            arrayOfRtrFileYAxises = new String[arrayOfRtrFile[0].length];
            for (int i = 0; i < arrayOfRtrFileYAxises.length; ++i) {
                arrayOfRtrFileYAxises[i] = arrayOfRtrFile[0][i];
            }
            listOfYAxises = arrayOfRtrFileYAxises;
        }
        else {
            listOfYAxises = new String[] {""};
        }
        return listOfYAxises;
    }
    
    /**
     * Setting the range for the errors
     */
    private void setRangeIndexForErrors()
    {
        for(int i=0;i<arrayOfRavFile[0].length;i++)
        {
            if(arrayOfRavFile[0][i].equals("NR"))
            {
                RangeIndex = i;
                
                break;
            }
        }
    }
    
    /**
     * Prepare the series collection with the array of file data
     * also, the arrays for the error bars according to the specified cases
     * and pass it to the JFreechart liberary function to actually plot the data
     * Change the titles into more readable title for y-axis.
     */
    private void plot() {
        Object[] selectedValues =  listYAxises.getSelectedValues();
        String titleForYAxis = "", error = "", acf = "";
        String[] arrayOfErrorsPLUS = new String [selectedValues.length];
        String[] arrayOfErrorsMINUS = new String [selectedValues.length];
        acfIsNeeded = false;
        if (selectedValues.length > 1) {
            for (int i = 0; i < selectedValues.length - 1; ++i) {
                String titleSelected = (String) selectedValues[i];
                titleForYAxis += titleSelected + ", ";
                if (buttonRES.isEnabled()) {
                    arrayOfErrorsPLUS [i] = "+" + getError(titleSelected);
                    arrayOfErrorsMINUS [i] = "-" + getError(titleSelected);
                    error += "+" + arrayOfErrorsMINUS [i] + ", ";
                    
                }
            }
        }
        titleForYAxis += (String) selectedValues[selectedValues.length - 1];
        if ((rbRTR.isSelected() || rbTHI.isSelected() || rbDCP.isSelected()) && buttonRES.isEnabled()) { ///// i added rbTHI check with rtr
            if (!titleForYAxis.startsWith("Int") && // Bad programming here. Needs better check!
                    (titleForYAxis.indexOf("D_") != -1 ||
                    titleForYAxis.indexOf("VS") != -1 ||
                    titleForYAxis.indexOf("VB") != -1) ||
                    titleForYAxis.indexOf("EC") != -1){
                acf = getAverageRTR();
                acfIsNeeded = true;
            }
        }
        if(titleXAxis.equals("Number of steps"))
        {
        titleForYAxis = titleForYAxis.replaceAll("NR", "Number of Steps");
        }
        else{
        titleForYAxis = titleForYAxis.replaceAll("NR", "Number of runs");}
        
        titleForYAxis = titleForYAxis.replaceAll("PROC", "Number of Processes");
        titleForYAxis = titleForYAxis.replaceAll("PRESS", "Reduced pressure");
        titleForYAxis = titleForYAxis.replaceAll("DENSITY", "Reduced density");
        titleForYAxis = titleForYAxis.replaceAll("TEMP", "Reduced temperature");
        titleForYAxis = titleForYAxis.replaceAll("EPOT", "Reduced potential energy");
        titleForYAxis = titleForYAxis.replaceAll("ENTLP", "Reduced enthalpy");
        titleForYAxis = titleForYAxis.replaceAll("MUE", "Chemical potential");
        titleForYAxis = titleForYAxis.replaceAll("FRACT", "Mole fraction");
        titleForYAxis = titleForYAxis.replaceAll("DISP", "Displacement");
        titleForYAxis = titleForYAxis.replaceAll("VW", "Partial Molar Volume");
        titleForYAxis = titleForYAxis.replaceAll("HM", "Partial Molar Enthalpy");
        titleForYAxis = titleForYAxis.replaceAll("IntD_i", "Self-diff. coeff.");
        titleForYAxis = titleForYAxis.replaceAll("Int VS", "Shear-Viscosity");
        titleForYAxis = titleForYAxis.replaceAll("Int VB", "Bulk-Viscosity");
        titleForYAxis = titleForYAxis.replaceAll("Int CO", "Thermal Conductivity");
        titleForYAxis = titleForYAxis.replaceAll("D_i", "ACF Self-diff. coeff.");
        titleForYAxis = titleForYAxis.replaceAll("Int C", "Thermal Conductivity");
        titleForYAxis = titleForYAxis.replaceAll("VS", "ACF Shear-Viscosity");
        titleForYAxis = titleForYAxis.replaceAll("VB", "ACF Bulk-Viscosity");
        titleForYAxis = titleForYAxis.replaceAll("CO", "ACF Thermal Conductivity");
        titleForYAxis = titleForYAxis.replaceAll("EP_Intra", "Reduced intramolecular potential energy");
        titleForYAxis = titleForYAxis.replaceAll("EP_Bonds", "Reduced bond energy");
        titleForYAxis = titleForYAxis.replaceAll("EP_Angles", "Reduced angle energy");
        titleForYAxis = titleForYAxis.replaceAll("EP_Dihed", "Reduced torsion energy");
        titleForYAxis = titleForYAxis.replaceAll("EP_14_15", "Reduced 1-4/1-5 energy");
        titleForYAxis = titleForYAxis.replaceAll("EP_Inter", "Reduced intermolecular potential energy");
        titleForYAxis = titleForYAxis.replaceAll("Vir_Intra", "Virial intramolecular energy");
        titleForYAxis = titleForYAxis.replaceAll("Vir_Inter", "Virial intermolecular energy");
        titleForYAxis = titleForYAxis.replaceAll("Int_", " ");
        titleForYAxis = titleForYAxis.replaceAll("Int", " ");
        XYSeriesCollection xyseriescollection = new XYSeriesCollection();
        if (buttonRES.isEnabled()) {
            arrayOfErrorsPLUS[selectedValues.length - 1] = "+" + getError((String) selectedValues[selectedValues.length - 1]);
            arrayOfErrorsMINUS[selectedValues.length - 1] = "-" + getError((String) selectedValues[selectedValues.length - 1]);
            error += "+" + arrayOfErrorsMINUS[selectedValues.length - 1];
            //JOptionPane.showMessageDialog(null, getError((String) selectedValues[selectedValues.length - 1]));
            tfError.setEnabled(true);
            tfError.setText(error);
        }
        if (rbRTR.isSelected()|| rbTHI.isSelected() || rbDCP.isSelected()) {///// i added rbTHI check with rtr
            oneByESeries = new XYSeries("1/e");
            averageSeries = new XYSeries("Time span");
            ZeroLine = new XYSeries("Zero");
            xyseriesForRtr = new XYSeries[columnsInRtrFileForYAxis.length][arrayOfRtrFileDatasetLengths.length];
////////..................
            errorseriesForRtr = new XYSeries[columnsInRtrFileForYAxis.length][2];
            if (acfIsNeeded) {
                for (int j = 0; j < arrayOfRtrFileDatasetLengths[arrayOfRtrFileDatasetLengths.length - 1]; ++j) {
                    try {
                        oneByESeries.add(
                                Double.parseDouble(
                                arrayOfRtrFile[arrayOfRtrFileStartingIndices[arrayOfRtrFileStartingIndices.length - 1] + 1 + j][0]),
                                0.367879441);
                        averageSeries.add(Double.parseDouble(acf), 0);// made change from -0.5 to 0
                        averageSeries.add(Double.parseDouble(acf), 1);
                        ZeroLine.add(Double.parseDouble(
                                arrayOfRtrFile[arrayOfRtrFileStartingIndices[arrayOfRtrFileStartingIndices.length - 1] + 1 + j][0]),0);
                    } catch (Exception e) {
                        System.out.println("Incorrect data point at: " + j + " 1/e");
                        oneByESeries.add(
                                Double.parseDouble(
                                arrayOfRtrFile[arrayOfRtrFileStartingIndices[arrayOfRtrFileStartingIndices.length - 1] + 1 + j][0]),
                                null);
                        averageSeries.add(null, 0);// made change from -0.5 to 0
                        averageSeries.add(null, 1);
                    }
                }
            }
            if (acfIsNeeded) {
                xyseriescollection.addSeries(oneByESeries);
                xyseriescollection.addSeries(ZeroLine);
                xyseriescollection.addSeries(averageSeries);
            }
            for (int i = 0; i < columnsInRtrFileForYAxis.length; ++i) {
                if(rbDCP.isSelected())
                    xyseriesForRtr[i][xyseriesForRtr[i].length - 1] = new XYSeries("DCP file");
                else
                if(rbTHI.isSelected())
                    xyseriesForRtr[i][xyseriesForRtr[i].length - 1] = new XYSeries("THI file");
                else
                xyseriesForRtr[i][xyseriesForRtr[i].length - 1] = new XYSeries("RTR file");
                ///................                
                errorseriesForRtr[i][0] = new XYSeries("(+) Error");
                errorseriesForRtr[i][1] = new XYSeries("(-) Error");
               ////////............
                for (int j = 0; j < arrayOfRtrFileDatasetLengths[arrayOfRtrFileDatasetLengths.length - 1]; ++j) {
                    try {
                        xyseriesForRtr[i][xyseriesForRtr[i].length - 1].add(
                                Double.parseDouble(
                                arrayOfRtrFile[arrayOfRtrFileStartingIndices[arrayOfRtrFileStartingIndices.length - 1] + 1 + j][columnInRtrFileForXAxis]),
                                Double.parseDouble(
                                arrayOfRtrFile[arrayOfRtrFileStartingIndices[arrayOfRtrFileStartingIndices.length - 1] + 1 + j][columnsInRtrFileForYAxis[i]]));
                    } catch (Exception e) {
                        System.out.println("Incorrect data point at: " + j + " for column " + i);
                        xyseriesForRtr[i][xyseriesForRtr[i].length - 1].add(
                                Double.parseDouble(
                                arrayOfRtrFile[arrayOfRtrFileStartingIndices[arrayOfRtrFileStartingIndices.length - 1] + 1 + j][columnInRtrFileForXAxis]),
                                null);
                    }
                }
                /////////////............
                if (buttonRES.isEnabled() && !arrayOfErrorsPLUS[i].startsWith("+ ")) {
                    
                    for (int j = 0; j < arrayOfRtrFileDatasetLengths[arrayOfRtrFileDatasetLengths.length - 1]; ++j) {
                        try {
                            errorseriesForRtr[i][0].add(
                                    Double.parseDouble(
                                    arrayOfRtrFile[arrayOfRtrFileStartingIndices[arrayOfRtrFileStartingIndices.length - 1] + 1 + j][0]),
                                    Double.parseDouble(
                                    arrayOfRtrFile[arrayOfRtrFile.length - 1][columnsInRtrFileForYAxis[i]]) + Double.parseDouble(arrayOfErrorsPLUS[i]));
                            errorseriesForRtr[i][1].add(
                                    Double.parseDouble(
                                    arrayOfRtrFile[arrayOfRtrFileStartingIndices[arrayOfRtrFileStartingIndices.length - 1] + 1 + j][0]),
                                    Double.parseDouble(
                                    arrayOfRtrFile[arrayOfRtrFile.length - 1][columnsInRtrFileForYAxis[i]]) + Double.parseDouble(arrayOfErrorsMINUS[i]));

                        } catch (Exception e) {
                            System.out.println("Incorrect data point (for error) at: " + j + " for column " + i);
                           errorseriesForRtr[i][0].add(
                                    Double.parseDouble(
                                    arrayOfRtrFile[arrayOfRtrFileStartingIndices[arrayOfRtrFileStartingIndices.length - 1] + 1 + j][0]),
                                    null);
                            errorseriesForRtr[i][1].add(
                                    Double.parseDouble(
                                    arrayOfRtrFile[arrayOfRtrFileStartingIndices[arrayOfRtrFileStartingIndices.length - 1] + 1 + j][0]),
                                    null);
                        }
                    }
                }
                ////////////.......
                xyseriescollection.addSeries(xyseriesForRtr[i][xyseriesForRtr[i].length - 1]);
                //////////...............
               if (buttonRES.isEnabled() && !arrayOfErrorsPLUS[i].startsWith("+ ")) {
                
                    xyseriescollection.addSeries(errorseriesForRtr[i][0]);
                    xyseriescollection.addSeries(errorseriesForRtr[i][1]);
                }
                ////............
            }
        } else {
            if (cbRUN.isEnabled() && cbRUN.isSelected()) {
                xyseriesForRun = new XYSeries[columnsInRunFileForYAxis.length][arrayOfRunFileDatasetLengths.length];
                int offset;
                for (int i = 0; i < columnsInRunFileForYAxis.length; ++i) {
                    if (cbPROD.isEnabled() && cbPROD.isSelected()) {
                        offset = 1;
                        xyseriesForRun[i][xyseriesForRun[i].length - offset] = new XYSeries("Production (run)");
                        for (int j = 0; j < arrayOfRunFileDatasetLengths[arrayOfRunFileDatasetLengths.length - offset]; ++j) {
                            try {
                                xyseriesForRun[i][xyseriesForRun[i].length - offset].add(
                                        Double.parseDouble(
                                        arrayOfRunFile[arrayOfRunFileStartingIndices[arrayOfRunFileStartingIndices.length - offset] + 1 + j][columnInRunFileForXAxis]),
                                        Double.parseDouble(
                                        arrayOfRunFile[arrayOfRunFileStartingIndices[arrayOfRunFileStartingIndices.length - offset] + 1 + j][columnsInRunFileForYAxis[i]]));
                            } catch (Exception e) {
                                System.out.println("Incorrect data point at: " + j + " for column " + i);
                                xyseriesForRun[i][xyseriesForRun[i].length - offset].add(
                                        Double.parseDouble(
                                        arrayOfRunFile[arrayOfRunFileStartingIndices[arrayOfRunFileStartingIndices.length - offset] + 1 + j][columnInRunFileForXAxis]),
                                        null);
                            }
                        }
                        xyseriescollection.addSeries(xyseriesForRun[i][xyseriesForRun[i].length - offset]);
                    }
                    if (cbEQUI.isEnabled() && cbEQUI.isSelected()) {
                        offset = 2;
                        xyseriesForRun[i][xyseriesForRun[i].length - offset] = new XYSeries("Equilibration (run)");
                        for (int j = 0; j < arrayOfRunFileDatasetLengths[arrayOfRunFileDatasetLengths.length - offset]; ++j) {
                            try {
                                xyseriesForRun[i][xyseriesForRun[i].length - offset].add(
                                        Double.parseDouble(
                                        arrayOfRunFile[arrayOfRunFileStartingIndices[arrayOfRunFileStartingIndices.length - offset] + 1 + j][columnInRunFileForXAxis]),
                                        Double.parseDouble(
                                        arrayOfRunFile[arrayOfRunFileStartingIndices[arrayOfRunFileStartingIndices.length - offset] + 1 + j][columnsInRunFileForYAxis[i]]));
                            } catch (Exception e) {
                                System.out.println("Incorrect data point at: " + j + " for column " + i);
                                xyseriesForRun[i][xyseriesForRun[i].length - offset].add(
                                        Double.parseDouble(
                                        arrayOfRunFile[arrayOfRunFileStartingIndices[arrayOfRunFileStartingIndices.length - offset] + 1 + j][columnInRunFileForXAxis]),
                                        null);
                            }
                        }
                        xyseriescollection.addSeries(xyseriesForRun[i][xyseriesForRun[i].length - offset]);
                    }
                    if (cbNVTE.isEnabled() && cbNVTE.isSelected()) {
                        if (cbEQUI.isEnabled()) {
                            offset = 3;
                        } else {
                            offset = 2;
                        }
                        xyseriesForRun[i][xyseriesForRun[i].length - offset] = new XYSeries("NVT Equilibration (run)");
                        for (int j = 0; j < arrayOfRunFileDatasetLengths[arrayOfRunFileDatasetLengths.length - offset]; ++j) {
                            try {
                                xyseriesForRun[i][xyseriesForRun[i].length - offset].add(
                                        Double.parseDouble(
                                        arrayOfRunFile[arrayOfRunFileStartingIndices[arrayOfRunFileStartingIndices.length - offset] + 1 + j][columnInRunFileForXAxis]),
                                        Double.parseDouble(
                                        arrayOfRunFile[arrayOfRunFileStartingIndices[arrayOfRunFileStartingIndices.length - offset] + 1 + j][columnsInRunFileForYAxis[i]]));
                            } catch (Exception e) {
                                System.out.println("Incorrect data point at: " + j + " for column " + i);
                                xyseriesForRun[i][xyseriesForRun[i].length - offset].add(
                                        Double.parseDouble(
                                        arrayOfRunFile[arrayOfRunFileStartingIndices[arrayOfRunFileStartingIndices.length - offset] + 1 + j][columnInRunFileForXAxis]),
                                        null);
                            }
                        }
                        xyseriescollection.addSeries(xyseriesForRun[i][xyseriesForRun[i].length - offset]);
                    }
                }
            }
            if (cbRAV.isEnabled() && cbRAV.isSelected()) {
                xyseriesForRav = new XYSeries[columnsInRavFileForYAxis.length][arrayOfRavFileDatasetLengths.length];
                errorseriesForRav = new XYSeries[columnsInRavFileForYAxis.length][2];
                int offset;
                for (int i = 0; i < columnsInRavFileForYAxis.length; ++i) {
                    if (cbPROD.isEnabled() && cbPROD.isSelected()) {
                        offset = 1;
                        xyseriesForRav[i][xyseriesForRav[i].length - offset] = new XYSeries("Production (rav)");
                        errorseriesForRav[i][0] = new XYSeries("(+) Error");
                        errorseriesForRav[i][1] = new XYSeries("(-) Error");
                        for (int j = 0; j < arrayOfRavFileDatasetLengths[arrayOfRavFileDatasetLengths.length - offset]; ++j) {
                            try {
                                xyseriesForRav[i][xyseriesForRav[i].length - offset].add(
                                        Double.parseDouble(
                                        arrayOfRavFile[arrayOfRavFileStartingIndices[arrayOfRavFileStartingIndices.length - offset] + 1 + j][columnInRavFileForXAxis]),
                                        Double.parseDouble(
                                        arrayOfRavFile[arrayOfRavFileStartingIndices[arrayOfRavFileStartingIndices.length - offset] + 1 + j][columnsInRavFileForYAxis[i]]));
                            } catch (Exception e) {
                                System.out.println("Incorrect data point at: " + j + " for column " + i);
                                xyseriesForRav[i][xyseriesForRav[i].length - offset].add(
                                        Double.parseDouble(
                                        arrayOfRavFile[arrayOfRavFileStartingIndices[arrayOfRavFileStartingIndices.length - offset] + 1 + j][columnInRavFileForXAxis]),
                                        null);
                            }
                        }
                        if (buttonRES.isEnabled() && !arrayOfErrorsPLUS[i].startsWith("+ ")) {
                            for (int j = 0; j < arrayOfRavFileDatasetLengths[arrayOfRavFileDatasetLengths.length - offset]; ++j) {
                                try {
                                    errorseriesForRav[i][0].add(
                                            Double.parseDouble(
                                            arrayOfRavFile[arrayOfRavFileStartingIndices[arrayOfRavFileStartingIndices.length - offset] + 1 + j][RangeIndex]),//changed this instead of 0
                                            Double.parseDouble(
                                            arrayOfRavFile[arrayOfRavFile.length - 1][columnsInRavFileForYAxis[i]]) + Double.parseDouble(arrayOfErrorsPLUS[i]));
                                    errorseriesForRav[i][1].add(
                                            Double.parseDouble(
                                            arrayOfRavFile[arrayOfRavFileStartingIndices[arrayOfRavFileStartingIndices.length - offset] + 1 + j][RangeIndex]),//changed this instead of 0
                                            Double.parseDouble(
                                            arrayOfRavFile[arrayOfRavFile.length - 1][columnsInRavFileForYAxis[i]]) + Double.parseDouble(arrayOfErrorsMINUS[i]));
                                } catch (Exception e) {
                                    System.out.println("Incorrect data point (for error) at: " + j + " for column " + i);
                                    errorseriesForRav[i][0].add(
                                            Double.parseDouble(
                                            arrayOfRavFile[arrayOfRavFileStartingIndices[arrayOfRavFileStartingIndices.length - offset] + 1 + j][columnInRavFileForXAxis]),
                                            null);
                                    errorseriesForRav[i][1].add(
                                            Double.parseDouble(
                                            arrayOfRavFile[arrayOfRavFileStartingIndices[arrayOfRavFileStartingIndices.length - offset] + 1 + j][columnInRavFileForXAxis]),
                                            null);
                                }
                            }
                        }
                        xyseriescollection.addSeries(xyseriesForRav[i][xyseriesForRav[i].length - offset]);
                        if (buttonRES.isEnabled() && !arrayOfErrorsPLUS[i].startsWith("+ ")) {
                            xyseriescollection.addSeries(errorseriesForRav[i][0]);
                            xyseriescollection.addSeries(errorseriesForRav[i][1]);
                        }
                    }
                    if (cbEQUI.isEnabled() && cbEQUI.isSelected()) {
                        offset = 2;
                        xyseriesForRav[i][xyseriesForRav[i].length - offset] = new XYSeries("Equilibration (rav)");
                        for (int j = 0; j < arrayOfRavFileDatasetLengths[arrayOfRavFileDatasetLengths.length - offset]; ++j) {
                            try {
                                xyseriesForRav[i][xyseriesForRav[i].length - offset].add(
                                        Double.parseDouble(
                                        arrayOfRavFile[arrayOfRavFileStartingIndices[arrayOfRavFileStartingIndices.length - offset] + 1 + j][columnInRavFileForXAxis]),
                                        Double.parseDouble(
                                        arrayOfRavFile[arrayOfRavFileStartingIndices[arrayOfRavFileStartingIndices.length - offset] + 1 + j][columnsInRavFileForYAxis[i]]));
                            } catch (Exception e) {
                                System.out.println("Incorrect data point at: " + j + " for column " + i);
                                xyseriesForRav[i][xyseriesForRav[i].length - offset].add(
                                        Double.parseDouble(
                                        arrayOfRavFile[arrayOfRavFileStartingIndices[arrayOfRavFileStartingIndices.length - offset] + 1 + j][columnInRavFileForXAxis]),
                                        null);
                            }
                        }
                        xyseriescollection.addSeries(xyseriesForRav[i][xyseriesForRav[i].length - offset]);
                    }
                    if (cbNVTE.isEnabled() && cbNVTE.isSelected()) {
                        if (cbEQUI.isEnabled()) {
                            offset = 3;
                        } else {
                            offset = 2;
                        }
                        xyseriesForRav[i][xyseriesForRav[i].length - offset] = new XYSeries("NVT Equilibration (rav)");
                        for (int j = 0; j < arrayOfRavFileDatasetLengths[arrayOfRavFileDatasetLengths.length - offset]; ++j) {
                            try {
                                xyseriesForRav[i][xyseriesForRav[i].length - offset].add(
                                        Double.parseDouble(
                                        arrayOfRavFile[arrayOfRavFileStartingIndices[arrayOfRavFileStartingIndices.length - offset] + 1 + j][columnInRavFileForXAxis]),
                                        Double.parseDouble(
                                        arrayOfRavFile[arrayOfRavFileStartingIndices[arrayOfRavFileStartingIndices.length - offset] + 1 + j][columnsInRavFileForYAxis[i]]));
                            } catch (Exception e) {
                                System.out.println("Incorrect data point at: " + j + " for column " + i);
                                xyseriesForRav[i][xyseriesForRav[i].length - offset].add(
                                        Double.parseDouble(
                                        arrayOfRavFile[arrayOfRavFileStartingIndices[arrayOfRavFileStartingIndices.length - offset] + 1 + j][columnInRavFileForXAxis]),
                                        null);
                            }
                        }
                        xyseriescollection.addSeries(xyseriesForRav[i][xyseriesForRav[i].length - offset]);
                    }
                }
            }
        }
        
        JFreeChart chart = ChartFactory.createXYLineChart(null,
                /*(String) comboXAxis.getSelectedItem()*/titleXAxis,
                titleForYAxis,
                xyseriescollection,
                org.jfree.chart.plot.PlotOrientation.VERTICAL,
                true,
                false,
                false);
        //JOptionPane.showMessageDialog(null, chart.getLegend());
               
        jPanelOfChart = new ChartPanel(chart, true);
        panelChartContainer.removeAll();
        panelChartContainer.add(jPanelOfChart);
        panelChartContainer.setVisible(false);
        panelChartContainer.setVisible(true);
        XYPlot xyplot = chart.getXYPlot();
        rangeAxis = (NumberAxis)xyplot.getRangeAxis();
        domainAxis = (NumberAxis)xyplot.getDomainAxis();
        rangeAxis.setAutoRangeIncludesZero(false);
        XYTextAnnotation textAnnotation;
        XYLineAndShapeRenderer rr = (XYLineAndShapeRenderer)xyplot.getRenderer();
        Color green = new Color(0, 255, 0);
        Color blue = new Color(0, 0, 255);
        Color red = new Color(255, 0, 0);
        Color black = new Color(0, 0, 0);
        BasicStroke stroke1 = new BasicStroke(1.0F, 1, 1, 1.0F, new float[] {2.0F, 6F}, 0.0F);
        BasicStroke stroke2 = new BasicStroke(2.0F);
        int k = 0;
        try {
            
            if (rbRTR.isSelected() || rbTHI.isSelected() || rbDCP.isSelected()) {///// i added rbTHI check with rtr
                if (acfIsNeeded) {
                    int j = arrayOfRtrFileDatasetLengths[arrayOfRtrFileDatasetLengths.length - 1];
                    textAnnotation = new XYTextAnnotation("1/e",
                            Double.parseDouble(
                            arrayOfRtrFile[arrayOfRtrFileStartingIndices[arrayOfRtrFileStartingIndices.length - 1] + j][0]),
                            0.367879441);
                    xyplot.addAnnotation(textAnnotation);
                    rr.setSeriesPaint(k, black, true);
                    ++k;
                    textAnnotation = new XYTextAnnotation("0",
                            Double.parseDouble(
                            arrayOfRtrFile[arrayOfRtrFileStartingIndices[arrayOfRtrFileStartingIndices.length - 1] + j][0]),
                            0);
                    xyplot.addAnnotation(textAnnotation);
                    rr.setSeriesPaint(k, black, true);
                    ++k;
                    textAnnotation = new XYTextAnnotation("ACF",
                            Double.parseDouble(acf),
                            (domainAxis.getLowerBound() + domainAxis.getUpperBound()) / 2);
                    xyplot.addAnnotation(textAnnotation);
                    rr.setSeriesPaint(k, black, true);
                    ++k;
                }
                for (int i = 0; i < columnsInRtrFileForYAxis.length; ++i) {
                    int j = arrayOfRtrFileDatasetLengths[arrayOfRtrFileDatasetLengths.length - 1];
                    textAnnotation = new XYTextAnnotation((String) selectedValues [i],
                            Double.parseDouble(
                            arrayOfRtrFile[arrayOfRtrFileStartingIndices[arrayOfRtrFileStartingIndices.length - 1] + j][columnInRtrFileForXAxis]),
                            Double.parseDouble(
                            arrayOfRtrFile[arrayOfRtrFileStartingIndices[arrayOfRtrFileStartingIndices.length - 1] + j][columnsInRtrFileForYAxis[i]]));
                    xyplot.addAnnotation(textAnnotation);
                    rr.setSeriesStroke(k, stroke2, true);
                    rr.setSeriesPaint(k, red, true);
                    ++k;
                    if (buttonRES.isEnabled() && !arrayOfErrorsPLUS[i].startsWith("+ ")) {
                    
                        textAnnotation = new XYTextAnnotation((String) selectedValues [i] + " (+ error)",
                                Double.parseDouble(
                                arrayOfRtrFile[arrayOfRtrFileStartingIndices[arrayOfRtrFileStartingIndices.length - 1] + j][0]),
                               Double.parseDouble(
                                arrayOfRtrFile[arrayOfRtrFile.length - 1][columnsInRtrFileForYAxis[i]]) + Double.parseDouble(arrayOfErrorsPLUS[i]));
                       xyplot.addAnnotation(textAnnotation);
                       rr.setSeriesPaint(k, black, true);
                       ++k;
                       textAnnotation = new XYTextAnnotation((String) selectedValues [i] + " (- error)",
                                Double.parseDouble(
                               arrayOfRtrFile[arrayOfRtrFileStartingIndices[arrayOfRtrFileStartingIndices.length - 1] + j][0]),
                                Double.parseDouble(
                               arrayOfRtrFile[arrayOfRtrFile.length - 1][columnsInRtrFileForYAxis[i]]) + Double.parseDouble(arrayOfErrorsMINUS[i]));
                        xyplot.addAnnotation(textAnnotation);
                        rr.setSeriesPaint(k, black, true);
                        ++k;
                    }
                }
            } else {
                if (cbRUN.isEnabled() && cbRUN.isSelected()) {
                    int offset;
                    for (int i = 0; i < columnsInRunFileForYAxis.length; ++i) {
                        if (cbPROD.isEnabled() && cbPROD.isSelected()) {
                            offset = 1;
                            int j = arrayOfRunFileDatasetLengths[arrayOfRunFileDatasetLengths.length - offset];
                            textAnnotation = new XYTextAnnotation((String) selectedValues [i],
                                    Double.parseDouble(
                                    arrayOfRunFile[arrayOfRunFileStartingIndices[arrayOfRunFileStartingIndices.length - offset] + j][columnInRunFileForXAxis]),
                                    Double.parseDouble(
                                    arrayOfRunFile[arrayOfRunFileStartingIndices[arrayOfRunFileStartingIndices.length - offset] + j][columnsInRunFileForYAxis[i]]));
                            xyplot.addAnnotation(textAnnotation);
                            rr.setSeriesStroke(k, stroke1, true);
                            rr.setSeriesPaint(k, red, true);
                            ++k;
                        }
                        if (cbEQUI.isEnabled() && cbEQUI.isSelected()) {
                            offset = 2;
                            int j = arrayOfRunFileDatasetLengths[arrayOfRunFileDatasetLengths.length - offset];
                            textAnnotation = new XYTextAnnotation((String) selectedValues [i],
                                    Double.parseDouble(
                                    arrayOfRunFile[arrayOfRunFileStartingIndices[arrayOfRunFileStartingIndices.length - offset] + j][columnInRunFileForXAxis]),
                                    Double.parseDouble(
                                    arrayOfRunFile[arrayOfRunFileStartingIndices[arrayOfRunFileStartingIndices.length - offset] + j][columnsInRunFileForYAxis[i]]));
                            xyplot.addAnnotation(textAnnotation);
                            rr.setSeriesStroke(k, stroke1, true);
                            rr.setSeriesPaint(k, blue, true);
                            ++k;
                        }
                        if (cbNVTE.isEnabled() && cbNVTE.isSelected()) {
                            if (cbEQUI.isEnabled()) {
                                offset = 3;
                            } else {
                                offset = 2;
                            }
                            int j = arrayOfRunFileDatasetLengths[arrayOfRunFileDatasetLengths.length - offset];
                            textAnnotation = new XYTextAnnotation((String) selectedValues [i],
                                    Double.parseDouble(
                                    arrayOfRunFile[arrayOfRunFileStartingIndices[arrayOfRunFileStartingIndices.length - offset] + j][columnInRunFileForXAxis]),
                                    Double.parseDouble(
                                    arrayOfRunFile[arrayOfRunFileStartingIndices[arrayOfRunFileStartingIndices.length - offset] + j][columnsInRunFileForYAxis[i]]));
                            xyplot.addAnnotation(textAnnotation);
                            rr.setSeriesStroke(k, stroke1, true);
                            rr.setSeriesPaint(k, green, true);
                            ++k;
                        }
                    }
                }
                if (cbRAV.isEnabled() && cbRAV.isSelected()) {
                    int offset;
                    for (int i = 0; i < columnsInRavFileForYAxis.length; ++i) {
                        if (cbPROD.isEnabled() && cbPROD.isSelected()) {
                            offset = 1;
                            int j = arrayOfRavFileDatasetLengths[arrayOfRavFileDatasetLengths.length - offset];
                            textAnnotation = new XYTextAnnotation((String) selectedValues [i],
                                    Double.parseDouble(
                                    arrayOfRavFile[arrayOfRavFileStartingIndices[arrayOfRavFileStartingIndices.length - offset] + j][columnInRavFileForXAxis]),
                                    Double.parseDouble(
                                    arrayOfRavFile[arrayOfRavFileStartingIndices[arrayOfRavFileStartingIndices.length - offset] + j][columnsInRavFileForYAxis[i]]));
                            xyplot.addAnnotation(textAnnotation);
                            rr.setSeriesStroke(k, stroke2, true);
                            rr.setSeriesPaint(k, red, true);
                            ++k;
                            if (buttonRES.isEnabled() && !arrayOfErrorsPLUS[i].startsWith("+ ")) {
                                textAnnotation = new XYTextAnnotation((String) selectedValues [i] + " (+ error)",
                                        Double.parseDouble(
                                        arrayOfRavFile[arrayOfRavFileStartingIndices[arrayOfRavFileStartingIndices.length - offset] + j][RangeIndex]),// changed from 0 to this
                                        Double.parseDouble(
                                        arrayOfRavFile[arrayOfRavFile.length - 1][columnsInRavFileForYAxis[i]]) + Double.parseDouble(arrayOfErrorsPLUS[i]));
                                xyplot.addAnnotation(textAnnotation);
                                rr.setSeriesPaint(k, black, true);
                                ++k;
                                textAnnotation = new XYTextAnnotation((String) selectedValues [i] + " (- error)",
                                        Double.parseDouble(
                                        arrayOfRavFile[arrayOfRavFileStartingIndices[arrayOfRavFileStartingIndices.length - offset] + j][RangeIndex]),// changed from 0 to this
                                        Double.parseDouble(
                                        arrayOfRavFile[arrayOfRavFile.length - 1][columnsInRavFileForYAxis[i]]) + Double.parseDouble(arrayOfErrorsMINUS[i]));
                                xyplot.addAnnotation(textAnnotation);
                                rr.setSeriesPaint(k, black, true);
                                ++k;
                            }
                        }
                        if (cbEQUI.isEnabled() && cbEQUI.isSelected()) {
                            offset = 2;
                            int j = arrayOfRavFileDatasetLengths[arrayOfRavFileDatasetLengths.length - offset];
                            textAnnotation = new XYTextAnnotation((String) selectedValues [i],
                                    Double.parseDouble(
                                    arrayOfRavFile[arrayOfRavFileStartingIndices[arrayOfRavFileStartingIndices.length - offset] + j][columnInRavFileForXAxis]),
                                    Double.parseDouble(
                                    arrayOfRavFile[arrayOfRavFileStartingIndices[arrayOfRavFileStartingIndices.length - offset] + j][columnsInRavFileForYAxis[i]]));
                            xyplot.addAnnotation(textAnnotation);
                            rr.setSeriesStroke(k, stroke2, true);
                            rr.setSeriesPaint(k, blue, true);
                            ++k;
                        }
                        if (cbNVTE.isEnabled() && cbNVTE.isSelected()) {
                            if (cbEQUI.isEnabled()) {
                                offset = 3;
                            } else {
                                offset = 2;
                            }
                            int j = arrayOfRavFileDatasetLengths[arrayOfRavFileDatasetLengths.length - offset];
                            textAnnotation = new XYTextAnnotation((String) selectedValues [i],
                                    Double.parseDouble(
                                    arrayOfRavFile[arrayOfRavFileStartingIndices[arrayOfRavFileStartingIndices.length - offset] + j][columnInRavFileForXAxis]),
                                    Double.parseDouble(
                                    arrayOfRavFile[arrayOfRavFileStartingIndices[arrayOfRavFileStartingIndices.length - offset] + j][columnsInRavFileForYAxis[i]]));
                            xyplot.addAnnotation(textAnnotation);
                            rr.setSeriesStroke(k, stroke2, true);
                            rr.setSeriesPaint(k, green, true);
                            ++k;
                        }
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        tfXMax.setText(Math.round(domainAxis.getRange().getUpperBound() * 100) / 100.0 + "");
        tfXMin.setText(Math.round(domainAxis.getRange().getLowerBound() * 100) / 100.0 + "");
        tfYMax.setText(Math.round(rangeAxis.getRange().getUpperBound() * 100) / 100.0 + "");
        tfYMin.setText(Math.round(rangeAxis.getRange().getLowerBound() * 100) / 100.0 + "");
        tfXMax.setEnabled(true);
        tfXMin.setEnabled(true);
        tfYMax.setEnabled(true);
        tfYMin.setEnabled(true);
        buttonViewData.setEnabled(true);
        buttonPNGsave.setEnabled(true);
        buttonXAxisChange.setEnabled(true);
        buttonYAxisChange.setEnabled(true);
        buttonClear.setEnabled(true);
    }
    
    /**
     * Reading the errors for specified title with the help of .res file.
     * 
     * @param titleSelected Selected title
     * 
     * @return Error
     */
    private String getError(String titleSelected) {
        String trailingTitle = null;
        String midline = null;
        System.out.println(titleSelected);
        //JOptionPane.showMessageDialog(null, titleSelected);
        if (titleSelected.startsWith("MUE") ||
                titleSelected.startsWith("FRACT") ||
                titleSelected.startsWith("VW") || titleSelected.startsWith("Int") ||
                /*titleSelected.startsWith("D_i")||*/ titleSelected.startsWith("HM")) {
            if(titleSelected.contains(" "))
            {
            StringTokenizer st = new StringTokenizer(titleSelected);
            
            titleSelected = st.nextToken();
            trailingTitle = st.nextToken();
            }
            if(titleSelected.equals("MUE_1") || titleSelected.equals("HM_1") || titleSelected.equals("VW_1"))
            {
                String[] st = titleSelected.split("_");
                titleSelected = st[0];
                trailingTitle = st[1];
            }
        }
        String reducedError = "";
        File inputFile;
        if (rbRTR.isSelected()) {
            inputFile = new File(selectedDirectory.getAbsolutePath() + java.io.File.separator + caseTitle + ".res");
        } else {
            inputFile = new File(selectedDirectory.getAbsolutePath() + java.io.File.separator + caseTitle + "_1.res");
        }
        try {
            BufferedReader inputStream = new BufferedReader(new FileReader(inputFile));
            
            inputStream.readLine();
            String fileInString;
            while(!(fileInString = inputStream.readLine()).contains("SIMULATION RESULT FILE"))
            {
                continue;
            }
           do {
                fileInString = inputStream.readLine();
            }while (!fileInString.endsWith("===="));
           if (rbRTR.isSelected()) {
                do {
                    fileInString = inputStream.readLine();
                }while (!fileInString.endsWith("===="));
            }
            String line;
            while (!(line = inputStream.readLine()).endsWith("====")) {
                fileInString += line + " ";
            }
            StringTokenizer st = new StringTokenizer(fileInString);
            String[] error2 = new String[0];
            String title = titleSelected;
            
            title = title.replaceAll("PRESS", "Pressure");
            title = title.replaceAll("DENSITY", "Density");
            title = title.replaceAll("TEMP", "Temperature");
            title = title.replaceAll("EPOT", "Potential");
            title = title.replaceAll("ENTLP", "Enthalpy");
            title = title.replaceAll("EP_Intra", "Intramolecular");
            title = title.replaceAll("EP_Bonds", "Bonds");
            title = title.replaceAll("EP_Angles", "Angles");
            title = title.replaceAll("EP_Dihed", "torsion");
            title = title.replaceAll("EP_14_15", "14/15");
            title = title.replaceAll("EP_Inter", "Intermolecular");
            title = title.replaceAll("MUE","Chemical");
            title = title.replaceAll("FRACT", "Molar");
            title = title.replaceAll("VW", "Partial");
            title = title.replaceAll("HM", "Partial");
            //title = title.replaceAll("int_Lij", "integrated");
            //title = title.replaceAll("CO", "conductivity"); requested for: Don't need error bar any more...
            if(titleSelected.equals("D_i"))
            {
            title = title.replaceAll("D_i", "Self-diff.");
            }
            if(titleSelected.equals("IntD_i"))
            {
            title = title.replaceAll("IntD_i", "Self-diff.");
            }
            if(titleSelected.equals("Int_Lij"))
            {
            title = title.replaceAll("Int_Lij", "integrated");
            }
            //title = title.replaceAll("VS", "Shear-Viscosity");
            if("VS".equals(trailingTitle))
            {
                title = title.replaceAll("Int", "Shear-Viscosity");
            }
            if("VB".equals(trailingTitle))
            {
                title = title.replaceAll("Int", "Bulk-Viscosity");
            }
            if("C".equals(trailingTitle))
            {
                title = title.replaceAll("Int", "conductivity");
            }
            //title = title.replaceAll("CO", "Thermal");
            if (title.equals(titleSelected)) {
                setLabelStatusMessage("The " + titleSelected + " column is selected in " + caseTitle + ".res");
                return " ";
            }
            //JOptionPane.showMessageDialog(null, title +" again | " + st.hasMoreTokens());
            int i =0;
            while (st.hasMoreTokens()) {
                line = st.nextToken();
                String delimiter;
                //JOptionPane.showMessageDialog(null, line + " : this is line, by the way");
                if (line.equals(title)) {
                    if (title.equals("Pressure") ||
                            title.equals("Density") ||
                            title.equals("Temperature") ||
                            title.equals("Potential") ||
                            title.equals("Enthalpy") ||
                            title.equals("Intramolecular") ||
                            title.equals("Bonds") ||
                            title.equals("Angles") ||
                            title.equals("torsion") ||
                            title.equals("14/15") ||
                            title.equals("Intermolecular") ||
                            title.equals("Shear-Viscosity") ||
                            title.equals("Bulk-Viscosity"))
                           // title.equals("integrated")) //title.equals("conductivity")) requested for: no need for error bars
                             {
                                 //JOptionPane.showMessageDialog(null, "equals");
                        do {
                            delimiter = st.nextToken();
                            //JOptionPane.showMessageDialog(null, delimiter +" 1");
                        }while (!delimiter.endsWith(":"));
                        if (rbRTR.isSelected()) {
                            do {
                                delimiter = st.nextToken();
                                //JOptionPane.showMessageDialog(null, delimiter+" 2");
                            }while (!delimiter.endsWith(":"));
                        }
                        st.nextToken();
                        reducedError = st.nextToken();
                        
                        break;
                    } else if (title.equals("Chemical") ||
                            title.equals("Molar") ||
                            title.equals("Partial") ||
                           // title.equals("integrated") ||
                            title.equals("Self-diff."))
                                    {
                        String[] tempArray = new String[error2.length + 1];
                        System.arraycopy(error2, 0, tempArray, 0, error2.length);
                        do {
                            delimiter = st.nextToken();
                            if(delimiter.equals("volume") || delimiter.equals("enthalpy"))
                            {
                                midline = delimiter;
                                //JOptionPane.showMessageDialog(null, "this line containes " + midline);
                            }
                            
                        }while (!(delimiter.endsWith(":")));
                        if (rbRTR.isSelected()) {
                            do {
                                delimiter = st.nextToken();
                               // JOptionPane.showMessageDialog(null, delimiter+" 3");
                            }while (!delimiter.endsWith(":"));
                        }
                        st.nextToken();
                        if(titleSelected.equals("VW") && midline.equals("volume"))
                        {
                       // tempArray[error2.length] = st.nextToken();
                           tempArray[i] = st.nextToken();
                           
                        //JOptionPane.showMessageDialog(null, "vm exec  "+ tempArray[i] +" : " + i);
                        i++;
                        }
                        else
                            if(titleSelected.equals("HM") && midline.equals("enthalpy"))
                        {
                            tempArray[i] = st.nextToken();
                            //JOptionPane.showMessageDialog(null, "hm exec  "+ tempArray[i]+" : " + i);
                            i++;
                        }
                        else
                                if(!titleSelected.equals("VW") && !titleSelected.equals("HM"))
                            {
                             tempArray[i] = st.nextToken();
                             //JOptionPane.showMessageDialog(null, "rest exec "+ tempArray[i]+" : " + i);
                             i++;
                            }
                        
                        error2 = tempArray;
                    }
                }
            }
            if (title.equals("Chemical") ||
                    title.equals("Molar") ||
                    title.equals("Partial") ||
                    //title.equals("integrated") ||
                    title.equals("Self-diff.")) {
                reducedError = error2[Integer.parseInt(trailingTitle) - 1];
                System.out.println(reducedError);
//                JOptionPane.showMessageDialog(null, Integer.parseInt(trailingTitle) - 1);
//                JOptionPane.showMessageDialog(null,"0:"+ error2[0].toString() + " | 1:" + error2[1].toString());
            }
        } catch (Exception e) {
            e.printStackTrace();
            return " ";
        }
        
        return reducedError;
    }
    
//    private String[] getErrorRTR(String titleSelected) {
//        String trailingTitle = null;
//        if (titleSelected.startsWith("D_i")) { //|| titleSelected.startsWith("Int D_i")) {
//            StringTokenizer st = new StringTokenizer(titleSelected);
//            titleSelected = st.nextToken();
//            trailingTitle = st.nextToken();
//        }
//        String[] reducedError = new String [2];
//        File inputFile = new File(selectedDirectory.getAbsolutePath() + java.io.File.separator + caseTitle + ".res");
//        try {
//            BufferedReader inputStream = new BufferedReader(new FileReader(inputFile));
//            inputStream.readLine();
//            String fileInString;
//            do
//            {
//                fileInString = inputStream.readLine();
//            }while (!fileInString.endsWith("===="));
//            do
//            {
//                fileInString = inputStream.readLine();
//            }while (!fileInString.endsWith("===="));
//            String line;
//            while (!(line = inputStream.readLine()).endsWith("====")) {
//                fileInString += line + " ";
//            }
//            StringTokenizer st = new StringTokenizer(fileInString);
//            String[] error2PLUS = new String[0];
//            String[] error2MINUS = new String[0];
//            String title = titleSelected;
//            title = title.replaceAll("D_i", "Self-diff.");
//            title = title.replaceAll("VS", "Shear-Viscosity");
//            title = title.replaceAll("VB", "Bulk-Viscosity");
//            title = title.replaceAll("CO", "Thermal");
//            if (title.equals(titleSelected) || title.startsWith("Int")) {
//                setLabelStatusMessage("No column for " + titleSelected + " in " + caseTitle + ".res");
//                reducedError [0] = " ";
//                reducedError [1] = " ";
//                return reducedError;
//            }
//            while (st.hasMoreTokens()) {
//                line = st.nextToken();
//                String delimiter;
//                if (line.equals(title)) {
//                    if (title.equals("Shear-Viscosity") ||
//                            title.equals("Bulk-Viscosity") ||
//                            title.equals("Thermal")) {
//                        do {
//                            delimiter = st.nextToken();
//                        }while (!delimiter.endsWith(":"));
//                        st.nextToken();
//                        reducedError [0] = st.nextToken();
//                        do {
//                            delimiter = st.nextToken();
//                        }while (!delimiter.endsWith(":"));
//                        st.nextToken();
//                        reducedError [1] = st.nextToken();
//                        break;
//                    } else if (title.equals("Self-diff.")) {
//                        String[] tempArrayPLUS = new String[error2PLUS.length + 1];
//                        String[] tempArrayMINUS = new String[error2MINUS.length + 1];
//                        System.arraycopy(error2PLUS, 0, tempArrayPLUS, 0, error2PLUS.length);
//                        System.arraycopy(error2MINUS, 0, tempArrayMINUS, 0, error2MINUS.length);
//                        do {
//                            delimiter = st.nextToken();
//                        }while (!(delimiter.endsWith(":")));
//                        st.nextToken();
//                        tempArrayPLUS[error2PLUS.length] = st.nextToken();
//                        do {
//                            delimiter = st.nextToken();
//                        }while (!(delimiter.endsWith(":")));
//                        st.nextToken();
//                        tempArrayMINUS[error2MINUS.length] = st.nextToken();
//                        error2PLUS = tempArrayPLUS;
//                        error2MINUS = tempArrayMINUS;
//                    }
//                }
//            }
//            if (title.equals("Self-diff.")) {
//                reducedError[0] = error2PLUS[Integer.parseInt(trailingTitle) - 1];
//                reducedError[1] = error2MINUS[Integer.parseInt(trailingTitle) - 1];
//            }
//        } catch (Exception e) {
//            e.printStackTrace();
//            reducedError [0] = " ";
//            reducedError [1] = " ";
//        }
//        return reducedError;
//    }
    
    /**
     * Read the value for the acf with the help of .res file
     * 
     * @return Value of acf
     */
    private String getAverageRTR() {
        String acf = " ";
        File inputFile = new File(selectedDirectory.getAbsolutePath() + java.io.File.separator + caseTitle + ".res");
        try {
            BufferedReader inputStream = new BufferedReader(new FileReader(inputFile));
            inputStream.readLine();
            String fileInString;
            while(!(fileInString = inputStream.readLine()).contains("SIMULATION RESULT FILE"))
            {
                continue;
            }
            
            do
            {
                fileInString = inputStream.readLine();
            }while (!fileInString.endsWith("===="));
            do
            {
                fileInString = inputStream.readLine();
            }while (!fileInString.endsWith("===="));
            String line;
            while (!(line = inputStream.readLine()).endsWith("====")) {
                fileInString += line + " ";
            }
            StringTokenizer st = new StringTokenizer(fileInString);
            while (st.hasMoreTokens()) {
                line = st.nextToken();
                if (line.equals("Time")) {
                    line = st.nextToken();
                    if (line.equals("span")) {
                        String delimiter;
                        do {
                            delimiter = st.nextToken();
                        }while (!delimiter.endsWith(":"));
                        do {
                            delimiter = st.nextToken();
                        }while (!delimiter.endsWith(":"));
                        acf = st.nextToken();
                        break;
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            acf = " ";
        }
        return acf;
    }
    
    private void setLabelStatusMessage(String statusMessage) {
        String newStatusMessage = statusMessage + " | " + labelStatus.getText();
        if (newStatusMessage.length() > 100) {
            newStatusMessage = newStatusMessage.substring(0, 100) + "...";
        }
        labelStatus.setText(newStatusMessage);
    }
    
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton Helpbtn;
    private javax.swing.JButton buttonBrowse;
    private javax.swing.JButton buttonClear;
    private javax.swing.JButton buttonDrawChart;
    private javax.swing.ButtonGroup buttonGroupFileType;
    private javax.swing.JButton buttonLOG;
    private javax.swing.JButton buttonPAR;
    private javax.swing.JButton buttonPNGsave;
    private javax.swing.JButton buttonQuit1;
    private javax.swing.JButton buttonRES;
    private javax.swing.JButton buttonRestart;
    private javax.swing.JButton buttonViewData;
    private javax.swing.JButton buttonXAxisChange;
    private javax.swing.JButton buttonYAxisChange;
    private javax.swing.JCheckBox cbAskDirectory;
    private javax.swing.JCheckBox cbAutoDraw;
    private javax.swing.JCheckBox cbEQUI;
    private javax.swing.JCheckBox cbNVTE;
    private javax.swing.JCheckBox cbPROD;
    private javax.swing.JCheckBox cbRAV;
    private javax.swing.JCheckBox cbRUN;
    private javax.swing.JComboBox comboCaseList;
    private javax.swing.JComboBox comboXAxis;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel10;
    private javax.swing.JLabel jLabel11;
    private javax.swing.JLabel jLabel12;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel6;
    private javax.swing.JLabel jLabel7;
    private javax.swing.JLabel jLabel8;
    private javax.swing.JLabel jLabel9;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JSeparator jSeparator1;
    private javax.swing.JSeparator jSeparator2;
    private javax.swing.JSeparator jSeparator3;
    private javax.swing.JLabel labelStatus;
    private javax.swing.JList listYAxises;
    private javax.swing.JPanel main_panel;
    private javax.swing.JScrollPane mainscroll;
    private javax.swing.JPanel panelChartContainer;
    private javax.swing.JPanel panelControl;
    private javax.swing.JPanel panelStatusBar;
    private javax.swing.JRadioButton rbDCP;
    private javax.swing.JRadioButton rbRTR;
    private javax.swing.JRadioButton rbRUN_RAV;
    private javax.swing.JRadioButton rbTHI;
    private javax.swing.JTextField tfError;
    private javax.swing.JTextField tfXMax;
    private javax.swing.JTextField tfXMin;
    private javax.swing.JTextField tfYMax;
    private javax.swing.JTextField tfYMin;
    private javax.swing.JButton vimtovmd_btn;
    // End of variables declaration//GEN-END:variables
    
    private boolean isErrorSeriesDrawable;
    
    boolean acfIsNeeded;
    
    private boolean isHavingCases;
    
    private ChartPanel jPanelOfChart;
    
    private File selectedDirectory;
    
    private static File startingDirectory;
    
    private int[] arrayOfRtrFileDatasetLengths;
    
    private int[] arrayOfRavFileDatasetLengths;
    
    private int[] arrayOfRunFileDatasetLengths;
    
    private int[] arrayOfRtrFileStartingIndices;
    
    private int[] arrayOfRavFileStartingIndices;
    
    private int[] arrayOfRunFileStartingIndices;
    
    private int columnInRtrFileForXAxis = -1;
    
    private int columnInRavFileForXAxis = -1;
    
    private int columnInRunFileForXAxis = -1;
    
    private int[] columnsInRtrFileForYAxis = new int[0];
    
    private int[] columnsInRavFileForYAxis = new int[0];
    
    private int[] columnsInRunFileForYAxis = new int[0];
    
    private NumberAxis domainAxis;
    
    private NumberAxis rangeAxis;
    
    private String[][] arrayOfRtrFile;
    
    private String[][] arrayOfRavFile;
    
    private String[][] arrayOfRunFile;
    
    private String caseTitle;
    
    String UserDirectory = null;
    
    private XYSeries[][] errorseriesForRav;
    
    private XYSeries[][] errorseriesForRtr;
    
    private XYSeries[][] xyseriesForRtr;
    
    private XYSeries[][] xyseriesForRav;
    
    private XYSeries[][] xyseriesForRun;
    
    private XYSeries averageSeries;
    private int RangeIndex = 0;
    private XYSeries ZeroLine;
    private XYSeries oneByESeries;

    /**
     *
     */
    public HashMap<String,String> helpMap = new HashMap<String,String>();
    
}
