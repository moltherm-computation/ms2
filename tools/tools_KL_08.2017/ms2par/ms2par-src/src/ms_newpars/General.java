package ms_newpars;

/*import java.lang.String;*/

/**
 * This class contains the getters and setters of general data fields in the ms2par
 * @author Syed Ahsan Ali
 */



public class General {
    private int RSteps = 100;

    private double RMinRadius = 1.8;

    private double RMaxRadius = 6.0;

    private int NOrient = 1000;

    private String systemOfUnits = "SI";

    private double lengthUnit = 3.5;

    private double energyUnit = 100.0;

    private double massUnit = 40.0;

    private String typeOfEnsemble = "NVT";

    private String typeOfSimulation = "MD";

    private String integrator = "Gear";
    
    private double timeStep = 2.0;

    private String timeStepUnit = "femtosec";
    
    private int mcorSteps = 0;

    private double acceptance = 0.5;

    private int nvtSteps = 20000;

    private int nptSteps = 20000;
    
    private int muevtsteps = 20000;

    private int runSteps = 100000;

    private int resultFreq = 1000;

    private int errorFreq = 5000;

    private int visualFreq = 0;
    
    private int NumDinBins = 0;
    
    private int WallTime = 20160;
    
    private int rdffreq = 0;
    
    private int numshell = 200;
    
    private double kappa = 5.6;
    
    private int NVecMax = 1000;
    
    private int NsqMax = 27;
    
    private int nMax = 5;
    
    private String intDegFreedType = "Off";
    
    private String intraLJ_El = "Off";
    
    private String LJ_El_14 = "Off" ;
    
    private String longRangeCorrType = "ReactionField";
    
    private String cutOffType = "COM";
    
    private String qshaker = "Off";

    /**
     *
     */
    public General() {
    }

    /**
     *
     * @return
     */
    public String getSystemOfUnits() {
        return systemOfUnits;
    }

    /**
     *
     * @param val
     */
    public void setSystemOfUnits(String val) {
        this.systemOfUnits = val;
    }

    /**
     *
     * @return
     */
    public double getLengthUnit() {
        return lengthUnit;
    }

    /**
     *
     * @param val
     */
    public void setLengthUnit(double val) {
        this.lengthUnit = val;
    }

    /**
     *
     * @return
     */
    public double getEnergyUnit() {
        return energyUnit;
    }

    /**
     *
     * @param val
     */
    public void setEnergyUnit(double val) {
        this.energyUnit = val;
    }

    /**
     *
     * @return
     */
    public double getMassUnit() {
        return massUnit;
    }

    /**
     *
     * @param val
     */
    public void setMassUnit(double val) {
        this.massUnit = val;
    }

    /**
     *
     * @return
     */
    public String getTypeOfEnsemble() {
        return typeOfEnsemble;
    }

    /**
     *
     * @param val
     */
    public void setTypeOfEnsemble(String val) {
        this.typeOfEnsemble = val;
    }

    /**
     *
     * @return
     */
    public String getTypeOfSimulation() {
        return typeOfSimulation;
    }

    /**
     *
     * @param val
     */
    public void setTypeOfSimulation(String val) {
        this.typeOfSimulation = val;
    }

    /**
     *
     * @return
     */
    public String getIntegrator() {
        return integrator;
    }

    /**
     *
     * @param val
     */
    public void setIntegrator(String val) {
        this.integrator = val;
    }

    /**
     *
     * @return
     */
    public double getTimeStep() {
        return timeStep;
    }

    /**
     *
     * @param val
     */
    public void setTimeStep(double val) {
        this.timeStep = val;
    }
    
    /**
     *
     * @param val
     * @return
     */
    public double getReducedTimeStep(double val) {
        double reduced;
        reduced = val / 1000.0 * Math.sqrt(this.getEnergyUnit()/(this.getMassUnit()*this.getLengthUnit()*this.getLengthUnit()));
        return Math.round(reduced*10000.)/10000.;
    }

    /**
     *
     * @return
     */
    public String getTimeStepUnit() {
        return timeStepUnit;
    }

    /**
     *
     * @param val
     */
    public void setTimeStepUnit(String val) {
        this.timeStepUnit = val;
    }
    
    /**
     *
     * @return
     */
    public int getMcorSteps() {
        return mcorSteps;
    }

    /**
     *
     * @param val
     */
    public void setMcorSteps(int val) {
        this.mcorSteps = val;
    }

    /**
     *
     * @return
     */
    public double getAcceptance() {
        return acceptance;
    }

    /**
     *
     * @param val
     */
    public void setAcceptance(double val) {
        this.acceptance = val;
    }

    /**
     *
     * @return
     */
    public int getNvtSteps() {
        return nvtSteps;
    }

    /**
     *
     * @param val
     */
    public void setNvtSteps(int val) {
        this.nvtSteps = val;
    }

    /**
     *
     * @return
     */
    public int getNptSteps() {
        return nptSteps;
    }

    /**
     *
     * @param val
     */
    public void setNptSteps(int val) {
        this.nptSteps = val;
    }

    /**
     *
     * @return
     */
    public int getmuevttSteps() {
        return muevtsteps;
    }

    /**
     *
     * @param val
     */
    public void setmuevttSteps(int val) {
        this.muevtsteps = val;
    }

    /**
     *
     * @return
     */
    public int getRunSteps() {
        return runSteps;
    }

    /**
     *
     * @param val
     */
    public void setRunSteps(int val) {
        this.runSteps = val;
    }

    /**
     *
     * @return
     */
    public int getResultFreq() {
        return resultFreq;
    }

    /**
     *
     * @param val
     */
    public void setResultFreq(int val) {
        this.resultFreq = val;
    }

    /**
     *
     * @return
     */
    public int getNumDinBins() {
        return NumDinBins;
    }

    /**
     *
     * @param val
     */
    public void setNumDinBins(int val) {
        this.NumDinBins = val;
    }
    
    /**
     *
     * @return
     */
    public int getWallTime() {
        return WallTime;
    }

    /**
     *
     * @param val
     */
    public void setWallTime(int val) {
        this.WallTime = val;
    }

    /**
     *
     * @return
     */
    public int getErrorFreq() {
        return errorFreq;
    }

    /**
     *
     * @param val
     */
    public void setErrorFreq(int val) {
        this.errorFreq = val;
    }

    /**
     *
     * @return
     */
    public int getVisualFreq() {
        return visualFreq;
    }

    /**
     *
     * @param val
     */
    public void setVisualFreq(int val) {
        this.visualFreq = val;
    }

    /**
     *
     * @return
     */
    public int getrdffreq() {
        return rdffreq;
    }

    /**
     *
     * @param val
     */
    public void setrdffreq(int val) {
        this.rdffreq = val;
    }

    /**
     *
     * @return
     */
    public int getnumshell() {
        return numshell;
    }

    /**
     *
     * @param val
     */
    public void setnumshell(int val) {
        this.numshell = val;
    }

    /**
     *
     * @return
     */
    public String getCutOffType() {
        return cutOffType;
    }

    /**
     *
     * @param val
     */
    public void setCutOffType(String val) {
        this.cutOffType = val;
    }

    /**
     *
     * @return
     */
    public int getNOrient() {
	return NOrient;
    }

    /**
     *
     * @return
     */
    public double getRMaxRadius() {
	return RMaxRadius;
    }

    /**
     *
     * @return
     */
    public double getRMinRadius() {
	return RMinRadius;
    }

    /**
     *
     * @return
     */
    public int getRSteps() {
	return RSteps;
    }

    /**
     *
     * @param validateInputInt
     */
    public void setNOrient(int validateInputInt) {
	this.NOrient = validateInputInt;
    }

    /**
     *
     * @param validateInputDouble
     */
    public void setRMaxRadius(double validateInputDouble) {
	this.RMaxRadius = validateInputDouble;
    }

    /**
     *
     * @param validateInputDouble
     */
    public void setRMinRadius(double validateInputDouble) {
	this.RMinRadius = validateInputDouble;		
    }

    /**
     *
     * @param validateInputInt
     */
    public void setRSteps(int validateInputInt) {
	this.RSteps = validateInputInt;
    }

    /**
     *
     * @return
     */
    public String getLongRangeCorrType() {
        return longRangeCorrType;
    }

    /**
     *
     * @param val
     */
    public void setLongRangeCorrType(String val) {
        this.longRangeCorrType = val;
    }

    /**
     *
     * @param validateInputDouble
     */
    public void setKappa(double validateInputDouble) {
        this.kappa = validateInputDouble;
    }
    
    /**
     *
     * @return
     */
    public double getKappa() {
	return kappa;
    }

    /**
     *
     * @param validateInputInt
     */
    public void setNVecMax(int validateInputInt) {
        this.NVecMax = validateInputInt;
    }
    
    /**
     *
     * @return
     */
    public int getNVecMax() {
	return NVecMax;
    }

    /**
     *
     * @param validateInputInt
     */
    public void setNsqMax(int validateInputInt) {
        this.NsqMax = validateInputInt;
    }
    
    /**
     *
     * @return
     */
    public int getNsqMax() {
	return NsqMax;
    }

    /**
     *
     * @param validateInputInt
     */
    public void setNMax(int validateInputInt) {
        this.nMax = validateInputInt;
    }
    
    /**
     *
     * @return
     */
    public int getNMax() {
	return nMax;
    }
    
    /**
     *
     * @return
     */
    public String getIntDegFreedType() {
        return intDegFreedType;
    }

    /**
     *
     * @param val
     */
    public void setIntDegFreedType(String val) {
        this.intDegFreedType = val;
    }
    
    /**
     *
     * @param val
     */
    public void setqshaker(String val)
    {
        this.qshaker = val;
    }
    
    /**
     *
     * @return
     */
    public String getqshaker()
    {
        return qshaker;
    }

    /**
     *
     * @return
     */
    public String getIntraLJ_El() {
        return intraLJ_El;
    }

    /**
     *
     * @param val
     */
    public void setIntraLJ_Ej(String val){
        this.intraLJ_El = val;     
    }
    
    /**
     *
     * @return
     */
    public String getLJ_El_14() {
        return LJ_El_14;
    }

    /**
     *
     * @param val
     */
    public void setLJ_El_14(String val){
        this.LJ_El_14 = val;     
    }
    
}
