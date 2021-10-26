package ms_newpars;

/*import java.lang.String;
import java.util.Collection;*/
import java.util.ArrayList;
import java.util.LinkedList;

/**
 *
 * @author Syed Ahsan Ali
 * 
 * This class has getters and setters for member variables to save all the
 * data privately for every Ensemble instance created seperately.
 * 
 * @see EnsemblePanel
 * 
 */
public class Ensemble {
    private double temperature = 273.15;
    
    private double hamiltonian = 9922.264335;
    
    private double enthalphy = -14606.018;

    private String temperatureUnit = "K";

    private double pressure = 0.0;

    private String pressureUnit = "MPa";

    private double density = 0.0;
    
    private double vapDensity;

    private double pistonMass = 0.0005;
    
    private String optPressure = "Yes";
    
    private String commonequi = "No" ;
    
    private String permeability = "No";

    private int nPart = 0;

    private int nTest = 0; //Used to read old par-Files
    
    /**
     *
     */
    public int numofcriteria = 0;
    
    /**
     *
     */
    public LinkedList<PotentialModel> potModelList;
    
    /**
     *
     */
    public ArrayList<ArrayList<Object>> hbondcriteria = null;

    private EtaXi etaxi = null; //Do not change default value!
    
    private EtaXi last = null; //Do not change default value!

    private double cutOff = -1.0;

    private double eps = 1.0E10;

    private double cutOffLJ = 0.0;

    private double cutOffDQ = 0.0;

    private double cutOffDD = 0.0;

    private double cutOffQQ = 0.0;

    private double gePressure0 = 0.0;

    private double geLiqDensity = 0.0;

    private double geLiqEnthalp = 0.0;

    private double geLiqBetaT = 0.0;

    private double geLiqdHdP = 0.0;

    private double geVarDensity = 0.0;

    private double geVarEnthalp = 0.0;

    private double geVarBetaT = 0.0;

    private double geVardHdP = 0.0;
    
    private int fluctFreq = 50;
    
    private int nFullFluct = 10;
    
    private int maxCounter = 0;

    private int corrLengthCF = 10000;
    
    private int resFreqCF = 5;

    private int spanFunCF = 300;
    
    private int stepFunCF = 5;

    private int viewFunCF = 100;
    
     private String transport = "Off";
     
    /**
     *
     */
    public String[] ComContainer;
   
    /**
     *
     */
    public ArrayList<String> Commentsdata;

    /**
     *
     */
    public ArrayList<String> tranportbox;

    /**
     *
     */
    public Ensemble() {
    }

    /**
     *
     * @return
     */
    public String getTransport() {
        return transport;
    }

    /**
     *
     * @param val
     */
    public void setTransport(String val) {
        this.transport = val;
    }

    /**
     *
     * @return
     */
    public double getTemperature() {
        return temperature;
    }

    /**
     *
     * @param val
     */
    public void setPermeability(String val) {
        this.permeability = val;
    }

    /**
     *
     * @return
     */
    public String getPermeability() {
        return permeability;
    }

    /**
     *
     * @param val
     */
    public void setTemperature(double val) {
        this.temperature = val;
    }
    
    /**
     *
     * @return
     */
    public double getHamiltonian() {
        return hamiltonian;
    }
    
    /**
     *
     * @param val
     */
    public void setHamiltonian(double val) {
        this.hamiltonian = val;
    }
    
    /**
     *
     * @return
     */
    public double getEnthalphy() {
        return enthalphy;
    }
    
    /**
     *
     * @param val
     */
    public void setEnthalphy(double val) {
        this.enthalphy = val;
    }

    /**
     *
     * @return
     */
    public String getTemperatureUnit() {
        return temperatureUnit;
    }

    /**
     *
     * @param val
     */
    public void setTemperatureUnit(String val) {
        this.temperatureUnit = val;
    }

    /**
     *
     * @param temperature
     * @return
     */
    public double CelsiusToKelvin(double temperature) {
        return temperature+273.15;
    }

    /**
     *
     * @param temperature
     * @return
     */
    public double KelvinToCelsius(double temperature) {
        return temperature-273.15;
    }

    /**
     *
     * @return
     */
    public double getPressure() {
        return pressure;
    }

    /**
     *
     * @param val
     */
    public void setPressure(double val) {
        this.pressure = val;
    }

    /**
     *
     * @return
     */
    public String getPressureUnit() {
        return pressureUnit;
    }

    /**
     *
     * @param val
     */
    public void setPressureUnit(String val) {
        this.pressureUnit = val;
    }

    /**
     *
     * @param pressure
     * @return
     */
    public double barToMPa(double pressure) {
        return pressure*0.1;
    }

    /**
     *
     * @param pressure
     * @return
     */
    public double MPaToBar(double pressure) {
        return pressure*10;
    }
    
    /**
     *
     * @return
     */
    public double getDensity() {
        return density;
    }

    /**
     *
     * @param val
     */
    public void setDensity(double val) {
        this.density = val;
    }
    
    /**
     *
     * @return
     */
    public double getVapDensity() {
        return vapDensity;
    }

    /**
     *
     * @param val
     */
    public void setVapDensity(double val) {
        this.vapDensity = val;
    }

    /**
     *
     * @return
     */
    public double getPistonMass() {
        return pistonMass;
    }

    /**
     *
     * @param val
     */
    public void setPistonMass(double val) {
        this.pistonMass = val;
    }

    /**
     *
     * @return
     */
    public int getNPart() {
        return nPart;
    }

    /**
     *
     * @param val
     */
    public void setNPart(int val) {
        this.nPart = val;
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
    public void setNTest(int val) {
        this.nTest = val;
    }
    
    /**
     *
     * @return
     */
    public LinkedList<PotentialModel> getPotModelList() {
        return potModelList;
    }

    /**
     *
     * @param val
     */
    public void setPotModelList(LinkedList<PotentialModel> val) {
        this.potModelList = val;
    }

    /**
     *
     * @return
     */
    public EtaXi getEtaxi() {
        return etaxi;
    }

    /**
     *
     * @param val
     */
    public void setEtaxi(EtaXi val) {
        this.etaxi = val;
    }
    
    /**
     *
     * @return
     */
    public EtaXi getLast() {
        return last;
    }

    /**
     *
     * @param val
     */
    public void setLast(EtaXi val) {
        this.last = val;
    }

    /**
     *
     * @return
     */
    public double getCutOff() {
        return cutOff;
    }

    /**
     *
     * @param val
     */
    public void setCutOff(double val) {
        this.cutOff = val;
    }

    /**
     *
     * @return
     */
    public double getEps() {
        return eps;
    }

    /**
     *
     * @param val
     */
    public void setEps(double val) {
        this.eps = val;
    }

    /**
     *
     * @return
     */
    public double getCutOffLJ() {
        return cutOffLJ;
    }

    /**
     *
     * @param val
     */
    public void setCutOffLJ(double val) {
        this.cutOffLJ = val;
    }

    /**
     *
     * @return
     */
    public double getCutOffDQ() {
        return cutOffDQ;
    }

    /**
     *
     * @param val
     */
    public void setCutOffDQ(double val) {
        this.cutOffDQ = val;
    }

    /**
     *
     * @return
     */
    public double getCutOffDD() {
        return cutOffDD;
    }

    /**
     *
     * @param val
     */
    public void setCutOffDD(double val) {
        this.cutOffDD = val;
    }

    /**
     *
     * @return
     */
    public double getCutOffQQ() {
        return cutOffQQ;
    }

    /**
     *
     * @param val
     */
    public void setCutOffQQ(double val) {
        this.cutOffQQ = val;
    }

    /**
     *
     * @return
     */
    public double getGePressure0() {
        return gePressure0;
    }

    /**
     *
     * @param val
     */
    public void setGePressure0(double val) {
        this.gePressure0 = val;
    }

    /**
     *
     * @return
     */
    public double getGeLiqDensity() {
        return geLiqDensity;
    }

    /**
     *
     * @param val
     */
    public void setGeLiqDensity(double val) {
        this.geLiqDensity = val;
    }

    /**
     *
     * @return
     */
    public double getGeLiqEnthalp() {
        return geLiqEnthalp;
    }

    /**
     *
     * @param val
     */
    public void setGeLiqEnthalp(double val) {
        this.geLiqEnthalp = val;
    }

    /**
     *
     * @return
     */
    public double getGeLiqBetaT() {
        return geLiqBetaT;
    }

    /**
     *
     * @param val
     */
    public void setGeLiqBetaT(double val) {
        this.geLiqBetaT = val;
    }

    /**
     *
     * @return
     */
    public double getGeLiqdHdP() {
        return geLiqdHdP;
    }

    /**
     *
     * @param val
     */
    public void setGeLiqdHdP(double val) {
        this.geLiqdHdP = val;
    }

    /**
     *
     * @return
     */
    public double getGeVarDensity() {
        return geVarDensity;
    }

    /**
     *
     * @param val
     */
    public void setGeVarDensity(double val) {
        this.geVarDensity = val;
    }

    /**
     *
     * @return
     */
    public double getGeVarEnthalp() {
        return geVarEnthalp;
    }

    /**
     *
     * @param val
     */
    public void setGeVarEnthalp(double val) {
        this.geVarEnthalp = val;
    }

    /**
     *
     * @return
     */
    public double getGeVarBetaT() {
        return geVarBetaT;
    }

    /**
     *
     * @param val
     */
    public void setGeVarBetaT(double val) {
        this.geVarBetaT = val;
    }

    /**
     *
     * @return
     */
    public double getGeVardHdP() {
        return geVardHdP;
    }

    /**
     *
     * @param val
     */
    public void setGeVardHdP(double val) {
        this.geVardHdP = val;
    }

    /**
     *
     * @return
     */
    public int getFluctFreq() {
        return fluctFreq;
    }

    /**
     *
     * @param val
     */
    public void setFluctFreq(int val) {
        this.fluctFreq = val;
    }
    
    /**
     *
     * @return
     */
    public int getNFullFluct() {
        return nFullFluct;
    }

    /**
     *
     * @param val
     */
    public void setNFullFluct(int val) {
        this.nFullFluct = val;
    }
    
    /**
     *
     * @return
     */
    public int getMaxCounter() {
        return maxCounter;
    }

    /**
     *
     * @param val
     */
    public void setMaxCounter(int val) {
        this.maxCounter = val;
    }

    /**
     *
     * @return
     */
    public int getresFreqCF() {
        return resFreqCF;
    }

    /**
     *
     * @param validateInputInt
     */
    public void setresFreqCF(int validateInputInt) {
	this.resFreqCF = validateInputInt;
    }

    /**
     *
     * @return
     */
    public int getcorrLengthCF() {
        return corrLengthCF;
    }

    /**
     *
     * @param validateInputInt
     */
    public void setcorrLengthCF(int validateInputInt) {
	this.corrLengthCF = validateInputInt;
    }

    /**
     *
     * @return
     */
    public int getspanFunCF() {
        return spanFunCF;
    }

    /**
     *
     * @param validateInputInt
     */
    public void setspanFunCF(int validateInputInt) {
	this.spanFunCF = validateInputInt;
    }
    
    /**
     *
     * @return
     */
    public int getstepFunCF() {
        return stepFunCF;
    }

    /**
     *
     * @param validateInputInt
     */
    public void setstepFunCF(int validateInputInt) {
	this.stepFunCF = validateInputInt;
    }

    /**
     *
     * @return
     */
    public int getviewFunCF() {
        return viewFunCF;
    }

    /**
     *
     * @param validateInputInt
     */
    public void setviewFunCF(int validateInputInt) {
	this.viewFunCF = validateInputInt;
    }
    
    /**
     *
     * @return
     */
    public String getOptPressure(){
        return optPressure;
    }
    
    /**
     *
     * @param val
     */
    public void setOptPressure (String val){
        this.optPressure = val;
    }

    /**
     *
     * @return
     */
    public String getcommonequi(){
        return commonequi;
    }
    
    /**
     *
     * @param val
     */
    public void setcommonequi (String val){
        this.commonequi = val;
    }
}