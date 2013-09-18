package inputfilegen;

import java.io.*;
/*import java.lang.String;*/
import javax.swing.JFrame;
import javax.swing.JOptionPane;


public class PotentialModel {
    private File file;
    
    private String path;            
    
    private String filename;
    
    private String name;
    
    private String pmFileData = "";

    private double molarFrac;
    
    private String chemPotMethod = "none";
    
    private int nTest = 0;
    
    private String weightFactorsType = "Guess";
    
    private String weightFactors = "Weighting factors for gradual insertion";

    private double geMolarFrac = 0.0;

    private double geLiqChemPot = 0.0;

    private double geLiqPartMolVol = 0.0;

    private double geVarChemPot = 0.0;

    private double geVarPartMolVol = 0.0;
    


    public PotentialModel(File pmfile, Double molfrac) {
        
        this.file = pmfile;
        this.filename = this.file.getName();
        this.path = this.file.getPath();
        this.path = this.path.substring(0,(this.path.length()-this.filename.length()));
        this.name = this.filename.substring(0,this.filename.length()-3);
        try {
            readPmFileData(this.file);
        } catch (IOException ex) {
            JOptionPane.showMessageDialog(new JFrame(), "Could not read from file: "+this.filename);
        }
        this.molarFrac = molfrac;
    }

    public String getName() {
        return name;
    }

    public void setName(String val) {
        this.name = val;
    }
    
    public String getFileName() {
        return filename;
    }

    public void setFileName(String val) {
        this.filename = val;
    }
    
    public String getPmFileData() {
        return pmFileData;
    }
    
    public void setPmFileData(String val) {
        this.pmFileData = val;
    }
    
    public double getMolarFrac() {
        return molarFrac;
    }

    public void setMolarFrac(double val) {
        this.molarFrac = val;
    }
    
    public String getChemPotMethod() {
        return chemPotMethod;
    }
    
    public void setChemPotMethod(String val) {
        this.chemPotMethod = val;
    }
    
    public int getNTest() {
        return nTest;
    }
    
    public void setNTest(int val){
        this.nTest = val; 
    }
    
    public String getWeightFactorsType() {
        return weightFactorsType;
    }
    
    public void setWeightFactorsType(String val) {
        this.weightFactorsType = val;
    }
    
    public String getWeightFactors() {
        return weightFactors;
    }
    
    public void setWeightFactors(String val) {
        this.weightFactors = val;
    }

    public double getGeMolarFrac() {
        return geMolarFrac;
    }

    public void setGeMolarFrac(double val) {
        this.geMolarFrac = val;
    }

    public double getGeLiqChemPot() {
        return geLiqChemPot;
    }

    public void setGeLiqChemPot(double val) {
        this.geLiqChemPot = val;
    }

    public double getGeLiqPartMolVol() {
        return geLiqPartMolVol;
    }

    public void setGeLiqPartMolVol(double val) {
        this.geLiqPartMolVol = val;
    }

    public double getGeVarChemPot() {
        return geVarChemPot;
    }

    public void setGeVarChemPot(double val) {
        this.geVarChemPot = val;
    }

    public double getGeVarPartMolVol() {
        return geVarPartMolVol;
    }

    public void setGeVarPartMolVol(double val) {
        this.geVarPartMolVol = val;
    }
    
    public File getFile(){
        return file;
    }
    
    public void setFile( File val){
        this.file = val;
    }
    
    
    private void readPmFileData(File file) throws IOException {
        
        FileInputStream input = new FileInputStream(file);
        BufferedReader in = new BufferedReader(new InputStreamReader(input));
        while(in.ready()){
            this.pmFileData += in.readLine();
            this.pmFileData += "\n";
        }
    }
    
}
