package inputfilegen;

/*import java.lang.String;*/

public class EtaXi {
    private EtaXi next;
    
    private boolean rowEnd = true;

    private double eta = 1.0;

    private double xi = 1.0;

    private PotentialModel pm1;

    private PotentialModel pm2;

    public EtaXi(PotentialModel rowPm, PotentialModel colPm) {
        this.pm1 = rowPm;
        this.pm2 = colPm;
    
    }

    public EtaXi getNext() {
        return next;
    }

    public void setNext(EtaXi val) {
        this.next = val;
    }
    
    public boolean getRowEnd() {
        return rowEnd;
    }
    
    public void setRowEnd(boolean val) {
        this.rowEnd = val;
    }

    public double getEta() {
        return eta;
    }

    public void setEta(double val) {
        this.eta = val;
    }

    public double getXi() {
        return xi;
    }

    public void setXi(double val) {
        this.xi = val;
    }

    public PotentialModel getPm1() {
        return pm1;
    }

    public void setPm1(PotentialModel val) {
        this.pm1 = val;
    }

    public PotentialModel getPm2() {
        return pm2;
    }

    public void setPm2(PotentialModel val) {
        this.pm2 = val;
    }
}
