package inputfilegen;

/*import java.lang.String;
import java.util.Collection;*/
import java.util.LinkedList;

public class Ensemble {
    private double temperature = 273.15;

    private String temperatureUnit = "K";

    private double pressure = 0.0;

    private String pressureUnit = "MPa";

    private double density = 0.0;
    
    private double vapDensity;

    private double pistonMass = 0.0005;
    
    private String optPressure = "Yes";

    private int nPart = 0;

    private int nTest = 0; //Used to read old par-Files
    
    public LinkedList<PotentialModel> potModelList;

    private EtaXi etaxi = null; //Do not change default value!
    
    private EtaXi last = null; //Do not change default value!

    private double cutOff = 0.0;

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

    private int viewFunCF = 100;

    public Ensemble() {
    }

    public double getTemperature() {
        return temperature;
    }

    public void setTemperature(double val) {
        this.temperature = val;
    }

    public String getTemperatureUnit() {
        return temperatureUnit;
    }

    public void setTemperatureUnit(String val) {
        this.temperatureUnit = val;
    }

    public double CelsiusToKelvin(double temperature) {
        return temperature+273.15;
    }
        public double KelvinToCelsius(double temperature) {
        return temperature-273.15;
    }

    public double getPressure() {
        return pressure;
    }

    public void setPressure(double val) {
        this.pressure = val;
    }

    public String getPressureUnit() {
        return pressureUnit;
    }

    public void setPressureUnit(String val) {
        this.pressureUnit = val;
    }

    public double barToMPa(double pressure) {
        return pressure*0.1;
    }

    public double MPaToBar(double pressure) {
        return pressure*10;
    }
    
    public double getDensity() {
        return density;
    }

    public void setDensity(double val) {
        this.density = val;
    }
    
    public double getVapDensity() {
        return vapDensity;
    }

    public void setVapDensity(double val) {
        this.vapDensity = val;
    }

    public double getPistonMass() {
        return pistonMass;
    }

    public void setPistonMass(double val) {
        this.pistonMass = val;
    }

    public int getNPart() {
        return nPart;
    }

    public void setNPart(int val) {
        this.nPart = val;
    }

    public int getNTest() {
        return nTest;
    }

    public void setNTest(int val) {
        this.nTest = val;
    }
    
    public LinkedList<PotentialModel> getPotModelList() {
        return potModelList;
    }

    public void setPotModelList(LinkedList<PotentialModel> val) {
        this.potModelList = val;
    }

    public EtaXi getEtaxi() {
        return etaxi;
    }

    public void setEtaxi(EtaXi val) {
        this.etaxi = val;
    }
    
    public EtaXi getLast() {
        return last;
    }

    public void setLast(EtaXi val) {
        this.last = val;
    }
    public double getCutOff() {
        return cutOff;
    }

    public void setCutOff(double val) {
        this.cutOff = val;
    }

    public double getEps() {
        return eps;
    }

    public void setEps(double val) {
        this.eps = val;
    }

    public double getCutOffLJ() {
        return cutOffLJ;
    }

    public void setCutOffLJ(double val) {
        this.cutOffLJ = val;
    }

    public double getCutOffDQ() {
        return cutOffDQ;
    }

    public void setCutOffDQ(double val) {
        this.cutOffDQ = val;
    }

    public double getCutOffDD() {
        return cutOffDD;
    }

    public void setCutOffDD(double val) {
        this.cutOffDD = val;
    }

    public double getCutOffQQ() {
        return cutOffQQ;
    }

    public void setCutOffQQ(double val) {
        this.cutOffQQ = val;
    }

    public double getGePressure0() {
        return gePressure0;
    }

    public void setGePressure0(double val) {
        this.gePressure0 = val;
    }

    public double getGeLiqDensity() {
        return geLiqDensity;
    }

    public void setGeLiqDensity(double val) {
        this.geLiqDensity = val;
    }

    public double getGeLiqEnthalp() {
        return geLiqEnthalp;
    }

    public void setGeLiqEnthalp(double val) {
        this.geLiqEnthalp = val;
    }

    public double getGeLiqBetaT() {
        return geLiqBetaT;
    }

    public void setGeLiqBetaT(double val) {
        this.geLiqBetaT = val;
    }

    public double getGeLiqdHdP() {
        return geLiqdHdP;
    }

    public void setGeLiqdHdP(double val) {
        this.geLiqdHdP = val;
    }

    public double getGeVarDensity() {
        return geVarDensity;
    }

    public void setGeVarDensity(double val) {
        this.geVarDensity = val;
    }

    public double getGeVarEnthalp() {
        return geVarEnthalp;
    }

    public void setGeVarEnthalp(double val) {
        this.geVarEnthalp = val;
    }

    public double getGeVarBetaT() {
        return geVarBetaT;
    }

    public void setGeVarBetaT(double val) {
        this.geVarBetaT = val;
    }

    public double getGeVardHdP() {
        return geVardHdP;
    }

    public void setGeVardHdP(double val) {
        this.geVardHdP = val;
    }

    public int getFluctFreq() {
        return fluctFreq;
    }

    public void setFluctFreq(int val) {
        this.fluctFreq = val;
    }
    
    public int getNFullFluct() {
        return nFullFluct;
    }

    public void setNFullFluct(int val) {
        this.nFullFluct = val;
    }
    
    public int getMaxCounter() {
        return maxCounter;
    }

    public void setMaxCounter(int val) {
        this.maxCounter = val;
    }

    public int getresFreqCF() {
        return resFreqCF;
    }

    public void setresFreqCF(int validateInputInt) {
	this.resFreqCF = validateInputInt;
    }

    public int getcorrLengthCF() {
        return corrLengthCF;
    }

    public void setcorrLengthCF(int validateInputInt) {
	this.corrLengthCF = validateInputInt;
    }

    public int getspanFunCF() {
        return spanFunCF;
    }

    public void setspanFunCF(int validateInputInt) {
	this.spanFunCF = validateInputInt;
    }

    public int getviewFunCF() {
        return viewFunCF;
    }

    public void setviewFunCF(int validateInputInt) {
	this.viewFunCF = validateInputInt;
    }
    
    public String getOptPressure(){
        return optPressure;
    }
    
    public void setOptPressure (String val){
        this.optPressure = val;
    }

}