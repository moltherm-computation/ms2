/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ms_newpars;

import javax.swing.JFrame;
import javax.swing.JOptionPane;

/**
 *
 * @author Syed Ahsan Ali
 */
public class TherCondparam {
    
    private PotentialModel pm;
    
    /**
     *
     * @param curPm
     */
    public TherCondparam(PotentialModel curPm) {
        this.pm = curPm;
          try {
              Double input = validateInputDouble(JOptionPane.showInputDialog("part. molar enthalpy:",this.pm.getpartmolarenthalpy()));
              this.pm.setpartmolarenthalpy(input);
            } catch(Exception e) {
                    
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
    
}
