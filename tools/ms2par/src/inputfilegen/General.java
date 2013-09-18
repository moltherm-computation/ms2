package inputfilegen;

/*import java.lang.String;*/


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

    private int nvtSteps = 5000;

    private int nptSteps = 20000;

    private int runSteps = 100000;

    private int resultFreq = 100;

    private int errorFreq = 5000;

    private int visualFreq = 0;

    private String cutOffType = "COM";

    private String transport = "Off";

    public General() {
    }

    public String getSystemOfUnits() {
        return systemOfUnits;
    }

    public void setSystemOfUnits(String val) {
        this.systemOfUnits = val;
    }

    public double getLengthUnit() {
        return lengthUnit;
    }

    public void setLengthUnit(double val) {
        this.lengthUnit = val;
    }

    public double getEnergyUnit() {
        return energyUnit;
    }

    public void setEnergyUnit(double val) {
        this.energyUnit = val;
    }

    public double getMassUnit() {
        return massUnit;
    }

    public void setMassUnit(double val) {
        this.massUnit = val;
    }

    public String getTypeOfEnsemble() {
        return typeOfEnsemble;
    }

    public void setTypeOfEnsemble(String val) {
        this.typeOfEnsemble = val;
    }

    public String getTypeOfSimulation() {
        return typeOfSimulation;
    }

    public void setTypeOfSimulation(String val) {
        this.typeOfSimulation = val;
    }

    public String getIntegrator() {
        return integrator;
    }

    public void setIntegrator(String val) {
        this.integrator = val;
    }

    public double getTimeStep() {
        return timeStep;
    }

    public void setTimeStep(double val) {
        this.timeStep = val;
    }
    
    public double getReducedTimeStep(double val) {
        double reduced;
        reduced = val / 1000.0 * Math.sqrt(this.getEnergyUnit()/(this.getMassUnit()*this.getLengthUnit()*this.getLengthUnit()));
        return Math.round(reduced*10000.)/10000.;
    }

    public String getTimeStepUnit() {
        return timeStepUnit;
    }

    public void setTimeStepUnit(String val) {
        this.timeStepUnit = val;
    }
    
    public int getMcorSteps() {
        return mcorSteps;
    }

    public void setMcorSteps(int val) {
        this.mcorSteps = val;
    }
    public double getAcceptance() {
        return acceptance;
    }

    public void setAcceptance(double val) {
        this.acceptance = val;
    }

    public int getNvtSteps() {
        return nvtSteps;
    }

    public void setNvtSteps(int val) {
        this.nvtSteps = val;
    }

    public int getNptSteps() {
        return nptSteps;
    }

    public void setNptSteps(int val) {
        this.nptSteps = val;
    }

    public int getRunSteps() {
        return runSteps;
    }

    public void setRunSteps(int val) {
        this.runSteps = val;
    }

    public int getResultFreq() {
        return resultFreq;
    }

    public void setResultFreq(int val) {
        this.resultFreq = val;
    }

    public int getErrorFreq() {
        return errorFreq;
    }

    public void setErrorFreq(int val) {
        this.errorFreq = val;
    }

    public int getVisualFreq() {
        return visualFreq;
    }

    public void setVisualFreq(int val) {
        this.visualFreq = val;
    }

    public String getCutOffType() {
        return cutOffType;
    }

    public void setCutOffType(String val) {
        this.cutOffType = val;
    }

    public int getNOrient() {
	return NOrient;
    }

    public double getRMaxRadius() {
	return RMaxRadius;
    }

    public double getRMinRadius() {
	return RMinRadius;
    }

    public int getRSteps() {
	return RSteps;
    }

    public void setNOrient(int validateInputInt) {
	this.NOrient = validateInputInt;
    }

    public void setRMaxRadius(double validateInputDouble) {
	this.RMaxRadius = validateInputDouble;
    }

    public void setRMinRadius(double validateInputDouble) {
	this.RMinRadius = validateInputDouble;		
    }

    public void setRSteps(int validateInputInt) {
	this.RSteps = validateInputInt;
    }

    public String getTransport() {
        return transport;
    }

    public void setTransport(String val) {
        this.transport = val;
    }

}
