/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package molecular_modeling;

import java.util.ArrayList;

/**
 *
 * @author Syed Ahsan Ali
 * 
 * This class calculates the Bond, Angle and dihedral error
 * based on the information provided 
 * by the main class in the constructor.
 * 
 * bondlist: contains bonds e.g 1-2,2-3,3-4 etc
 * anglelist: contains angle.
 * dihedlist: contains dihedrals.
 * bondlength: bond length information from database
 * angle: information about angle from database
 * dihed, dihed const: Forceconstant, Gemma from database
 * coord: coordinates of all the connected atoms.
 */
public class BADError {
   
    private double bondError, angleError, dihedralError;
    private double[][] Coord;
    
    BADError(int[][] bondlist, int[][] anglelist, int[][] dihedlist, double[] bondlength,
            double[] angle, double[][] dihed,double[][] dihedconst, double[][] coord)
    {
        Coord = coord;
        
        CalBondError(bondlist, bondlength);
        CalAngleError(anglelist, angle);
        CalDihedralError(dihedlist, dihed, dihedconst);
    }
   
    /**
     * 
     * Calculate the error for Bond.
     * 
     * @param BondList  
     * @param BondLength
     * 
     */
   private void CalBondError(int[][] BondList, double[] BondLength)
   {
       int i,j;
        double rx,ry,rz,r2,r,r0,drE,BE = 0.0;
        
        for(int n=0;n<BondList.length;n++)
        {
            
            i = BondList[n][0] - 1;
            j = BondList[n][1] - 1;
            
            r0 = BondLength[n];
            
            rx = Coord[i][0] - Coord[j][0]; // X coordinate of i and j
            ry = Coord[i][1] - Coord[j][1]; // Y coordinate of i and j
            rz = Coord[i][2] - Coord[j][2]; // Z coordinate of i and j
            
            r2 = rx*rx + ry*ry + rz*rz;
            r  = Math.sqrt(r2);
            
            /*Extra steps for Error calculation*/
            drE = (r - r0)/r0;
            BE += Math.pow(drE, 2);
            /*End*/
            
        }
        
       bondError = BE / BondList.length;
   }
   
   /**
     * 
     * Calculate the error for Angle.
     * 
     * @param AngleList 
     * @param Angle
     * 
     */
   private void CalAngleError(int[][] AngleList, double[] Angle)
   {
       int i, j, k;
        double rijx, rijy, rijz, rkjx, rkjy, rkjz;
        double rij, rij2, rkj, rkj2;
        double angle, angle0, dangleE, AE = 0.0;
        double cosa;
        
        for(int n=0;n<AngleList.length;n++)
        {
            i = AngleList[n][0] - 1;
            j = AngleList[n][1] - 1;
            k = AngleList[n][2] - 1;
            
            angle0 = Angle[n];
            
            rijx = Coord[i][0] - Coord[j][0]; // X coordinate of i and j
            rijy = Coord[i][1] - Coord[j][1]; // Y coordinate of i and j
            rijz = Coord[i][2] - Coord[j][2]; // Z coordinate of i and j
            rkjx = Coord[k][0] - Coord[j][0]; // X coordinate of k and j
            rkjy = Coord[k][1] - Coord[j][1]; // Y coordinate of k and j
            rkjz = Coord[k][2] - Coord[j][2]; // Z coordinate of k and j
            
            rij2 = rijx*rijx + rijy*rijy + rijz*rijz;
            rkj2 = rkjx*rkjx + rkjy*rkjy + rkjz*rkjz;
            rij  = Math.sqrt(rij2);
            rkj  = Math.sqrt(rkj2);
            cosa = (rijx*rkjx+rijy*rkjy+rijz*rkjz)/(rij*rkj);
            
            if ( cosa > 1 ) cosa =  1;
            if ( cosa < -1) cosa = -1;
            angle = Math.acos(cosa);
            
            /*Extra steps for Error calculation*/
            dangleE = (angle - angle0)/angle0;
            AE += Math.pow(dangleE, 2);
            /*End*/
            
        }
        
       angleError = AE/AngleList.length;
   }
   
   /**
     * 
     * Calculate the dihedral error.
     * 
     * @param DihedList 
     * @param Dihed
     * @param DihedConst
     * 
     */
   private void CalDihedralError(int[][] DihedList, double[][] Dihed,double[][] DihedConst)
   {
       int i, j, k, l;
        double rijx, rijy, rijz, rjkx, rjky, rjkz, rklx, rkly, rklz;
        double rij, rij2, rkj, rjk2, rkl, rkl2;
        double rik2, rjl2, ril2;
        double axb, bxc, den, num, dden, dnum;
        double cosd, sind, signum;
        double f, fi, fj, fk, fl;
        double dihedral, c, DE = 0.0, dihedralE, minConst;
        
        for(int n=0;n<DihedList.length;n++)
        {
            i = DihedList[n][0] - 1;
            j = DihedList[n][1] - 1;
            k = DihedList[n][2] - 1;
            l = DihedList[n][3] - 1;
            molecularModeling.log.append(">>> Dihedral for error calculation is: "+i+"-"+j+"-"+k+"-"+l+"\n");
           double[] dihedral0 = new double[Dihed[n].length];
            double[] Forconst = new double[DihedConst[n].length];
        
            for(int d = 0; d < Dihed[n].length; d++)
            {
                dihedral0[d] =  Dihed[n][d];
            }
            for(int d = 0; d < DihedConst[n].length; d++)
            {
                Forconst[d] = DihedConst[n][d];
            }
            
            rijx = Coord[j][0] - Coord[i][0]; // X coordinate of i and j
            rijy = Coord[j][1] - Coord[i][1]; // Y coordinate of i and j
            rijz = Coord[j][2] - Coord[i][2]; // Z coordinate of i and j
            rjkx = Coord[k][0] - Coord[j][0]; // X coordinate of k and j
            rjky = Coord[k][1] - Coord[j][1]; // Y coordinate of k and j
            rjkz = Coord[k][2] - Coord[j][2]; // Z coordinate of k and j
            rklx = Coord[l][0] - Coord[k][0]; // X coordinate of k and l
            rkly = Coord[l][1] - Coord[k][1]; // Y coordinate of k and l
            rklz = Coord[l][2] - Coord[k][2]; // Z coordinate of k and l
            
            rij2 = rijx*rijx + rijy*rijy + rijz*rijz;
            rjk2 = rjkx*rjkx + rjky*rjky + rjkz*rjkz;
            rkl2 = rklx*rklx + rkly*rkly + rklz*rklz;
            rij  = Math.sqrt(rij2);
            rkj  = Math.sqrt(rjk2);
            rkl  = Math.sqrt(rkl2);
            
            rik2 = rijx*rjkx + rijy*rjky + rijz*rjkz;
            rjl2 = rjkx*rklx + rjky*rkly + rjkz*rklz;
            ril2 = rijx*rklx + rijy*rkly + rijz*rklz;

            axb = (rij2*rjk2) - (rik2*rik2);
            bxc = (rjk2*rkl2) - (rjl2*rjl2);
            den = axb*bxc;
            
            if ( den > 1E-20 )
      {
        den = Math.sqrt(den);

        num = (rik2*rjl2) - (ril2*rjk2);
        cosd = num/den;
        if ( cosd > 1 ) cosd =  1;
        if ( cosd < -1) cosd = -1;

        dihedral= Math.acos(cosd);
        
        signum = rijx*(rjky*rklz-rkly*rjkz)+rijy*(rjkz*rklx-rklz*rjkx)+rijz*(rjkx*rkly-rklx*rjky); 
        if (signum < 0) dihedral = -dihedral;
        
        /*Extra steps for Error calculation*/
            minConst = minDihedral(Forconst,dihedral0,dihedral);
            dihedralE = (dihedral - minConst)/minConst;
            DE += Math.pow(dihedralE, 2);
        /*End*/
      }
        }
        
       dihedralError = DE/DihedList.length;
   }

   /**
     * 
     * Calculate the minimum dihedral by which dihedral error is calculated.
     * 
     * @param forconst  
     * @param dihedral0
     * @param dihedral
     * 
     * @return double - Minimum Dihedral
     * 
     */
   private double minDihedral(double[] forconst, double[] dihedral0, double dihedral)
    {
        double Xmin = 0.0, E = 0.0, Emin = 0.0, ddx = 0.0;
        ArrayList<Double> XminE = new ArrayList<Double>();
        for(double x = -Math.PI; x <= Math.PI; x += Math.PI/1000)
        {
            E = 0.0;
            for(int m = 0; m < forconst.length; m++)
            {
               // E += forconst[m] * (1 + Math.cos((m * x) - dihedral0[m]));
                E += -m*forconst[m] * Math.sin((m * x) - dihedral0[m]);
            }
            
                        
            //if(x == -Math.PI || E < Emin)
            //if(x == -Math.PI || (E >= -Math.PI && E <= Math.PI) || E == 0)
            if(Math.round(E) == 0)
            {
                Emin = E;
               //XminE.clear();
                XminE.add(x);
            }
            /*if((E == Emin && Math.abs(x - XminE.get(XminE.size()-1)) > 0.05))
            {
                XminE.add(x);
            }*/
                       
        }
        
        if(!XminE.isEmpty())
        {
            double nearst = findNearest(XminE,dihedral);
            //molecularModeling.log.append("  Array size is: "+XminE.size()+" , array is :\n");
            
            for(int i=0;i<dihedral0.length;i++)
                molecularModeling.log.append("   ("+i+") "+dihedral0[i]+"\n");
            
            molecularModeling.log.append(">>> all the mins...\n");
            for(int i=0;i<XminE.size();i++)
                molecularModeling.log.append("   ("+i+") "+XminE.get(i)+"\n");
            
            Xmin = Math.round(nearst);
            molecularModeling.log.append("   Nearest value to "+dihedral+" is: "+nearst+"\n");
            
        }
         
            return Xmin;
    }

   /**
     * 
     * From many minimum dihedrals find the nearest one to dihedral
     * by which the dihedral error will be calculated.
     * 
     */
   private double findNearest(ArrayList<Double> Xmin,double dihedral)
   {
        double distance = Math.abs(Xmin.get(0) - dihedral);
        int idx = 0;
        for(int c = 1; c < Xmin.size(); c++){
            double cdistance = Math.abs(Xmin.get(c) - dihedral);
            if(cdistance < distance){
                idx = c;
                distance = cdistance;
                 }
        }
        
        return Xmin.get(idx);
   }
   
    /**
     * Calculates the following error:
     * Bonds Error, Angles Error, Dihedrals Error.
     * 
     * @return Array of Errors
     */
       public int[] getGeometryError()
   {
       int[] error =new int[3];
       error[0] = (int) (bondError * 100);
       error[1] = (int) (angleError * 100);
       error[2] = (int) (dihedralError * 100);
       
       return error;
   }
}
