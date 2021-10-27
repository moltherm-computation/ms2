package ms_newpars;

import java.io.*;
/*import java.lang.String;*/
import javax.swing.JFrame;
import javax.swing.JOptionPane;

/**
 * This class contains getters and setters for potential model
 * @author Syed Ahsan Ali
 */
public class PotentialModel {
    private File file;
    
    private String path;            
    
    private String filename;
    
    private String name;
    
    private String pmFileData = "";

    private double molarFrac;
    
    private String permeability = "No";
    
    private String chemPotMethod = "none";
    
    private int nTest = 0;
    
    private double lambdamin = 0.2;
    private double lambdamax = 1.0;
    private int nbins = 100;
    private double labdastepmax = 0.1;
    private int lambdaexp = 4;
    
    private String weightFactorsType = "Guess";
    
    private String weightFactors = "Weighting factors for gradual insertion";

    private double geMolarFrac = 0.0;

    private double geLiqChemPot = 0.0;

    private double geLiqPartMolVol = 0.0;

    private double geVarChemPot = 0.0;

    private double geVarPartMolVol = 0.0;
    
    private double partmolarenthalpy = -34.15698594;
    
    /**
     *
     * @param pmfile
     * @param molfrac
     * @param perm
     */
    public PotentialModel(File pmfile, Double molfrac, String perm) {
        
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
        this.permeability = perm;
    }

    /**
     *
     * @return
     */
    public String getName() {
        return name;
    }

    /**
     *
     * @param val
     */
    public void setName(String val) {
        this.name = val;
    }
    
    /**
     *
     * @return
     */
    public String getFileName() {
        return filename;
    }

    /**
     *
     * @param val
     */
    public void setFileName(String val) {
        this.filename = val;
    }
    
    /**
     *
     * @return
     */
    public String getPmFileData() {
        return pmFileData;
    }
    
    /**
     *
     * @param val
     */
    public void setPmFileData(String val) {
        this.pmFileData = val;
    }
    
    /**
     *
     * @return
     */
    public double getMolarFrac() {
        return molarFrac;
    }

    /**
     *
     * @param val
     */
    public void setMolarFrac(double val) {
        this.molarFrac = val;
    }

    /**
     *
     * @return
     */
    public String getpermeablity() {
        return this.permeability;
    }

    /**
     *
     * @param val
     */
    public void setpermeability(String val) {
        this.permeability = val;
    }
    
    /**
     *
     * @return
     */
    public double getpartmolarenthalpy() {
        return partmolarenthalpy;
    }
    
    /**
     *
     * @param val
     */
    public void setpartmolarenthalpy(double val) {
        this.partmolarenthalpy = val;
    }
    
    /**
     *
     * @return
     */
    public String getChemPotMethod() {
        return chemPotMethod;
    }
    
    /**
     *
     * @param val
     */
    public void setChemPotMethod(String val) {
        this.chemPotMethod = val;
    }
    
    /**
     *
     * @return
     */
    public int getNTest() {
        return nTest;
    }
    
    /**
     *
     * @param val
     */
    public void setNTest(int val){
        this.nTest = val; 
    }

    /**
     *
     * @param val
     */
    public void setlambdamin(double val){
        this.lambdamin = val;
    }

    /**
     *
     * @param val
     */
    public void setlambdamax(double val){
        this.lambdamax = val;
    }

    /**
     *
     * @param val
     */
    public void setnbins(int val){
        this.nbins = val;
    }

    /**
     *
     * @param val
     */
    public void setlabdastepmax(double val){
        this.labdastepmax = val;
    }

    /**
     *
     * @param val
     */
    public void setlambdaexp(int val){
        this.lambdaexp = val;
    }

    /**
     *
     * @return
     */
    public double getlambdamin()
    {
        return lambdamin;
    }

    /**
     *
     * @return
     */
    public double getlambdamax()
    {
        return lambdamax;
    }

    /**
     *
     * @return
     */
    public int getnbins()
    {
        return nbins;
    }

    /**
     *
     * @return
     */
    public double getlambdastepmax()
    {
        return labdastepmax;
    }

    /**
     *
     * @return
     */
    public int getlambdaexp()
    {
        return lambdaexp;
    }

    /**
     *
     * @return
     */
    public String getWeightFactorsType() {
        return weightFactorsType;
    }
    
    /**
     *
     * @param val
     */
    public void setWeightFactorsType(String val) {
        this.weightFactorsType = val;
    }
    
    /**
     *
     * @return
     */
    public String getWeightFactors() {
        return weightFactors;
    }
    
    /**
     *
     * @param val
     */
    public void setWeightFactors(String val) {
        this.weightFactors = val;
    }

    /**
     *
     * @return
     */
    public double getGeMolarFrac() {
        return geMolarFrac;
    }

    /**
     *
     * @param val
     */
    public void setGeMolarFrac(double val) {
        this.geMolarFrac = val;
    }

    /**
     *
     * @return
     */
    public double getGeLiqChemPot() {
        return geLiqChemPot;
    }

    /**
     *
     * @param val
     */
    public void setGeLiqChemPot(double val) {
        this.geLiqChemPot = val;
    }

    /**
     *
     * @return
     */
    public double getGeLiqPartMolVol() {
        return geLiqPartMolVol;
    }

    /**
     *
     * @param val
     */
    public void setGeLiqPartMolVol(double val) {
        this.geLiqPartMolVol = val;
    }

    /**
     *
     * @return
     */
    public double getGeVarChemPot() {
        return geVarChemPot;
    }

    /**
     *
     * @param val
     */
    public void setGeVarChemPot(double val) {
        this.geVarChemPot = val;
    }

    /**
     *
     * @return
     */
    public double getGeVarPartMolVol() {
        return geVarPartMolVol;
    }

    /**
     *
     * @param val
     */
    public void setGeVarPartMolVol(double val) {
        this.geVarPartMolVol = val;
    }
    
    /**
     *
     * @return
     */
    public File getFile(){
        return file;
    }
    
    /**
     *
     * @param val
     */
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
