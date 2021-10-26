package ms_newpars;

/*import java.lang.String;*/

/**
 *
 * @author Syed Ahsan Ali
 * 
 * This class is used to store eta and xi value for each specified potential model
 * by the help of getters and setters.
 * 
 * @see EtaXiPanel
 */


public class EtaXi {
    private EtaXi next;
    
    private boolean rowEnd = true;

    private double eta = 1.0;

    private double xi = 1.0;

    private PotentialModel pm1;

    private PotentialModel pm2;

    /**
     *
     * @param rowPm
     * @param colPm
     */
    public EtaXi(PotentialModel rowPm, PotentialModel colPm) {
        this.pm1 = rowPm;
        this.pm2 = colPm;
    
    }

    /**
     *
     * @return
     */
    public EtaXi getNext() {
        return next;
    }

    /**
     *
     * @param val
     */
    public void setNext(EtaXi val) {
        this.next = val;
    }
    
    /**
     *
     * @return
     */
    public boolean getRowEnd() {
        return rowEnd;
    }
    
    /**
     *
     * @param val
     */
    public void setRowEnd(boolean val) {
        this.rowEnd = val;
    }

    /**
     *
     * @return
     */
    public double getEta() {
        return eta;
    }

    /**
     *
     * @param val
     */
    public void setEta(double val) {
        this.eta = val;
    }

    /**
     *
     * @return
     */
    public double getXi() {
        return xi;
    }

    /**
     *
     * @param val
     */
    public void setXi(double val) {
        this.xi = val;
    }

    /**
     *
     * @return
     */
    public PotentialModel getPm1() {
        return pm1;
    }

    /**
     *
     * @param val
     */
    public void setPm1(PotentialModel val) {
        this.pm1 = val;
    }

    /**
     *
     * @return
     */
    public PotentialModel getPm2() {
        return pm2;
    }

    /**
     *
     * @param val
     */
    public void setPm2(PotentialModel val) {
        this.pm2 = val;
    }
}
