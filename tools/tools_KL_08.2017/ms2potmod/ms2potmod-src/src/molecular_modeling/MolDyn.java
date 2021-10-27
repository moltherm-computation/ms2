/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package molecular_modeling;

import java.util.Arrays;
import java.util.Random;

/**
 *
 * @author Syed Ahsan Ali
 * This is very important class. which calculates the dynamics of
 * molecular structure e.g calculating the Geometry. on the basis
 * of:
 * BondList, AngleList, DihedList: list of all the existing bonds
 * Angles and Dihedrals in the structure.
 * BondLength: Array of bond_length attributes from the database file for list of matched bonds.
 * BondConst: Array of forceconst attributes from the database file for list of matched bonds.
 * Angle: Array of bending_angle attributes from the database file for list of matched angles.
 * AngleConst: Array of bending_Constant attributes from the database file for list of matched angles.
 * Dihed: Array of torsion attributes from the database file for list of matched dihedrals.
 * DihedConst: Array of Gemma attributes from the database file for list of matched dihedrals.
 */
public class MolDyn {
    
    double[][] NetForce, Dihed, DihedConst, Coord, Vel;
    int[][] BondList, AngleList, DihedList;
    double[] BondLength, BondConst, Angle, AngleConst, Mass;
    boolean isBondSet = false, isAngleSet = false, isDihedSet = false;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
                
    }
    
    /**
     * Initialize the Bond Potential Arrays.
     * 
     * @see molecularModeling#CalculateGeometryActionPerformed(java.awt.event.ActionEvent)
     * 
     * @param bondlist Array of the bonds
     * @param bondlength Array of the bond length
     * @param bondconst Array of the bond constant
     * 
     */
    public void SetBondPot(int[][] bondlist, double[] bondlength, double[] bondconst)
    {
        BondList = bondlist;
        BondLength = bondlength;
        BondConst = bondconst;
        
        isBondSet = true;
    }

    /**
     * 
     * Initialize the Angle potential Arrays.
     * 
     * @see molecularModeling#CalculateGeometryActionPerformed(java.awt.event.ActionEvent)
     * 
     * @param anglelist Array of three connected atoms
     * @param angle Array of angle values from database
     * @param angleconst Array of angle constant
     * 
     */
    public void SetAnglePot(int[][] anglelist, double[] angle, double[] angleconst)
    {
        AngleList = anglelist;
        Angle = angle;
        AngleConst = angleconst;
        
        isAngleSet = true;
    }

    /**
     * Initialize the Dihedrals potential Arrays.
     * 
     * @see molecularModeling#CalculateGeometryActionPerformed(java.awt.event.ActionEvent)
     * 
     * @param dihedlist Array of four connected atoms
     * @param dihed Array of dihedrals value from database file
     * @param dihedconst Array of dihedral constant
     * 
     */
    public void SetDihedPot(int[][] dihedlist, double[][] dihed, double[][] dihedconst)
    {
        DihedList = dihedlist;
        Dihed = dihed;
        DihedConst = dihedconst;
        
        isDihedSet = true;
    }

    /**
     *
     * From here down is the algorithm for calculating geometry of the 
     * provided molecular structure.
     * 
     * @see molecularModeling#CalculateGeometryActionPerformed(java.awt.event.ActionEvent)
     * 
     * @param cord Array of Coordinates from the table
     * @param mass Array of mass values from the table
     * 
     * @return Array of new coordinates
     */
    
    public double[][] CalculateGeometry(double[][] cord, double[] mass)
    {
        if(isBondSet && isAngleSet && isDihedSet)
        {
            isBondSet = isAngleSet = isDihedSet = false;
            
            int NSteps=1000000;
           
            Coord = cord;
            Mass = mass;
            Vel = new double[cord.length][3];
            NetForce = new double[cord.length][3];

            for(int n=0;n<Angle.length;n++)
            {
            Angle[n] = Angle[n] * (Math.PI/180);
            }
            
            for(int n=0;n<Dihed.length;n++)
            {   
               for(int d = 0; d < Dihed[n].length; d++)
               {
               Dihed[n][d] =  Dihed[n][d] * (Math.PI/180);
               }
            }

            for(int i=0; i<Vel.length; i++)
            {
            Vel[i][0] = 0.5* Math.pow(-1,i)* randDouble(-1,1);
            Vel[i][1] = 0.5* Math.pow(-1,i)* randDouble(-1,1);
            Vel[i][2] = 0.5* Math.pow(-1,i)* randDouble(-1,1);
            }

            for (int i=0; i<BondConst.length;i++)
            {
               if (BondConst[i] == 0) BondConst[i] = 10000000;
            }

            for (int i=1; i<=NSteps; i++)
            {
               
                Interaction();
                Integrator();
            }
            
            return Coord;
        }
        else
            return null;
    }
    private void Interaction()
    {
       for (int i=0; i<NetForce.length; i++)
        {
            NetForce[i][0] = 0.0;
            NetForce[i][1] = 0.0;
            NetForce[i][2] = 0.0;
        }
        BondPot();
        AnglePot();
        DihedralPot();
    }

    /**
     * Generate random value b/w two numbers.
     * 
     * @param min 
     * @param max
     * 
     * @return Random value
     */
    public double randDouble(double min, double max) {

    Random random = new Random();

    double range = max - min;
    double scaled = random.nextDouble() * range;
    double shifted = scaled + min;
    
    return shifted;
}
    private void Integrator()
    {
        
        double TimeStep = 0.0001;
        for (int i=0; i<NetForce.length; i++)
        {
            if(Double.isNaN(NetForce[i][0]) || Double.isNaN(NetForce[i][1]) || Double.isNaN(NetForce[i][2]))
            {
                System.out.println("");
            }
           Vel[i][0]=0.9*Vel[i][0]+(NetForce[i][0]/Mass[i])*TimeStep;
           Vel[i][1]=0.9*Vel[i][1]+(NetForce[i][1]/Mass[i])*TimeStep;
           Vel[i][2]=0.9*Vel[i][2]+(NetForce[i][2]/Mass[i])*TimeStep;
           Coord[i][0]=Coord[i][0]+Vel[i][0]*TimeStep;
           Coord[i][1]=Coord[i][1]+Vel[i][1]*TimeStep;
           Coord[i][2]=Coord[i][2]+Vel[i][2]*TimeStep;
           if(Double.isNaN(Coord[i][0]) || Double.isNaN(Coord[i][1]) || Double.isNaN(Coord[i][2]))
            {
                System.out.println("");
            }
        }
    }
    private void BondPot()
    {
        int i,j;
        double rx,ry,rz,r2,r,dr,r0,Forconst,ff,BE = 0.0;
        
        for(int n=0;n<BondList.length;n++)
        {
            
            i = BondList[n][0] - 1;
            j = BondList[n][1] - 1;
            
            r0 = BondLength[n];
            Forconst = BondConst[n];
            
            rx = Coord[i][0] - Coord[j][0]; // X coordinate of i and j
            ry = Coord[i][1] - Coord[j][1]; // Y coordinate of i and j
            rz = Coord[i][2] - Coord[j][2]; // Z coordinate of i and j
            
            r2 = rx*rx + ry*ry + rz*rz;
            r  = Math.sqrt(r2);
            dr = r - r0;
            
            ff = -2.0*Forconst*dr/r;

            NetForce[i][0] += ff*rx; // X of i
            NetForce[i][1] += ff*ry; // Y of i
            NetForce[i][2] += ff*rz; // z of i
            NetForce[j][0] -= ff*rx; // X of j
            NetForce[j][1] -= ff*ry; // Y of j
            NetForce[j][2] -= ff*rz; // Z of j
            
        }
        
    }
    
    private void AnglePot()
    {
        int i, j, k;
        double rijx, rijy, rijz, rkjx, rkjy, rkjz;
        double rij, rij2, rkj, rkj2;
        double angle, angle0, dangle, Forconst, AE = 0.0;
        double fi, fj, fk;
        double cosa, sina, sab, tab;
        
        for(int n=0;n<AngleList.length;n++)
        {
            i = AngleList[n][0] - 1;
            j = AngleList[n][1] - 1;
            k = AngleList[n][2] - 1;
            
            angle0 = Angle[n];
            Forconst = AngleConst[n];
            
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
            dangle = angle - angle0;
            
            //TotalPotEn += Forconst*(dangle*dangle);

            sina = Math.sqrt(1-cosa*cosa);
            if ( sina < 1E-20 ) sina = 1E-20;
            sab = -2.0*Forconst*dangle/sina;
            tab = sab*cosa;

            fi = tab/rij2;
            fk = sab/(rij*rkj);
            fj = tab/rkj2;

            NetForce[i][0] += fi*rijx - fk*rkjx;
            NetForce[i][1] += fi*rijy - fk*rkjy;
            NetForce[i][2] += fi*rijz - fk*rkjz;
            NetForce[j][0] += (fk-fi)*rijx + (fk-fj)*rkjx;
            NetForce[j][1] += (fk-fi)*rijy + (fk-fj)*rkjy;
            NetForce[j][2] += (fk-fi)*rijz + (fk-fj)*rkjz;
            NetForce[k][0] += fj*rkjx - fk*rijx;
            NetForce[k][1] += fj*rkjy - fk*rijy;
            NetForce[k][2] += fj*rkjz - fk*rijz;
        }
      
    }
    private void DihedralPot()
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

        sind = Math.sin(dihedral);
        if ( sind > -1E-20 && sind < 1E-20)
        {
          if (sind < 0)
          { sind = -1E-20;}
          else
          { sind =  1E-20;}
        }

        f = 0.0;
        for (int m=0; m<Forconst.length; m++)
        {
           //c = Math.cos(m*dihedral-dihedral0[m]);
           //TotalPotEn += Forconst[m]* (1+c);
            c = Math.sin(m*dihedral-dihedral0[m]);
           f -= m*Forconst[m]*c;
        }

        axb = axb/den*cosd;
        bxc = bxc/den*cosd;
        f = f/den/sind;

        /* Forces to X-direction */
        dnum =  rklx*rjk2 - rjkx*rjl2;
        dden = ( rik2*rjkx - rijx*rjk2 ) * bxc;
        fi = (dnum - dden) * f;
        dnum = ((rjkx-rijx)*rjl2 - rik2*rklx ) + (2.0*ril2*rjkx - rklx*rjk2);
        dden = axb*(rjl2*rklx-rjkx*rkl2) + (rijx*rjk2-rij2*rjkx-rik2*(rjkx-rijx))*bxc;
        fj = (dnum - dden) * f;
        dnum = rik2*rjkx - rijx*rjk2;
        dden = axb*( rjk2*rklx - rjl2*rjkx );
        fl = (dnum - dden) * f;
        fk = -(fi+fj+fl);

        NetForce[i][0] += fi;
        NetForce[j][0] += fj;
        NetForce[k][0] += fk;
        NetForce[l][0] += fl;

       /* Forces to Y-direction */
        dnum =  rkly*rjk2 - rjky*rjl2;
        dden = ( rik2*rjky - rijy*rjk2 ) * bxc;
        fi = (dnum - dden) * f;
        dnum = ((rjky-rijy)*rjl2 - rik2*rkly ) + (2.0*ril2*rjky - rkly*rjk2);
        dden = axb*(rjl2*rkly-rjky*rkl2) + (rijy*rjk2-rij2*rjky-rik2*(rjky-rijy))*bxc;
        fj = (dnum - dden) * f;
        dnum = rik2*rjky - rijy*rjk2;
        dden = axb*( rjk2*rkly - rjl2*rjky );
        fl = (dnum - dden) * f;
        fk = -(fi+fj+fl);

        NetForce[i][1] += fi;
        NetForce[j][1] += fj;
        NetForce[k][1] += fk;
        NetForce[l][1] += fl;

        /* Forces to Z-direction */
        dnum =  rklz*rjk2 - rjkz*rjl2;
        dden = ( rik2*rjkz - rijz*rjk2 ) * bxc;
        fi = (dnum - dden) * f;
        dnum = ((rjkz-rijz)*rjl2 - rik2*rklz ) + (2.0*ril2*rjkz - rklz*rjk2);
        dden = axb*(rjl2*rklz-rjkz*rkl2) + (rijz*rjk2-rij2*rjkz-rik2*(rjkz-rijz))*bxc;
        fj = (dnum - dden) * f;
        dnum = rik2*rjkz - rijz*rjk2;
        dden = axb*( rjk2*rklz - rjl2*rjkz );
        fl = (dnum - dden) * f;
        fk = -(fi+fj+fl);

        NetForce[i][2] += fi;
        NetForce[j][2] += fj;
        NetForce[k][2] += fk;
        NetForce[l][2] += fl;
      }
            
        }
        
    }
    
    
}
