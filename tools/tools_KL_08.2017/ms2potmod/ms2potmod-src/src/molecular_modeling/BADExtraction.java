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
 * This is responsible for calculating the bonds, angles and dihedrals
 * from the IDs of the atom and links from the table...
 * After calculating the list of bonds it then calculates the angle
 * and dihedral from this list.
 * 
 */
public class BADExtraction {

private int[][] g_bonds;

    /**
     * 
     * @see molecularModeling#CalculateGeometryActionPerformed(java.awt.event.ActionEvent)
     * @see molecularModeling#WriteDBPropertiesInPmFile(java.io.BufferedWriter)
     * 
     * @param iDlinks Array of ids(index) and there links from the table
     * @return int[][] Array of calculated Bonds
     * 
     */
    public int[][] getBondsList(int[][] iDlinks)
{
   
    int length = iDlinks.length;
    
    int[][] bonds = new int[length - 1][2];
    int inex = 0;
    for(int i = 0; i < length; i++)
    {
        int id = iDlinks[i][0];
//        bonds.add(i,new ArrayList<Integer>());
        for(int j = 1; j <iDlinks[i].length; j++)
        {
            
            if(id < iDlinks[i][j])
            {
                bonds[inex][0] = id;
                bonds[inex][1] = iDlinks[i][j];
                inex++;
//                bonds.get(i).add(iDlinks[i][j]);
            }
        }
    }
    
    g_bonds = bonds;
    System.out.println("Bonds:\n");
    molecularModeling.log.append(">>>Bonds: "+g_bonds.length+"\n");
    for(int i = 0; i < g_bonds.length; i++)
    {
        System.out.println((g_bonds[i][0]-1)+"-"+(g_bonds[i][1]-1)+"  ");
        molecularModeling.log.append("\t"+(g_bonds[i][0]-1)+"-"+(g_bonds[i][1]-1)+"\n");
    }
    
    return bonds;
}
    /**
     * 
     * @see molecularModeling#CalculateGeometryActionPerformed(java.awt.event.ActionEvent)
     * @see molecularModeling#WriteDBPropertiesInPmFile(java.io.BufferedWriter)
     * 
     * @param iDlinks Array of ids(index) and there links from the table
     * @return int[][] Array of calculated Angles
     * 
     */
    public int[][] getAnglesList(int[][] iDlinks)
    {
        ArrayList<ArrayList<Integer>> Angles = new ArrayList<ArrayList<Integer>>();
        
        int length = g_bonds.length;
        int index = 0;
        
        for(int i = 0; i < length; i++)
        {
            
            int e1 = g_bonds[i][0];
            int e2 = g_bonds[i][1];
            
            int j = i + 1;
            while(j< length)
            {
                if(e2 == g_bonds[j][0] )//|| e1 == g_bonds[j][0]) // Either e1 equals or e2 for angle: newly added
                {
                    Angles.add(index,new ArrayList<Integer>());
                    Angles.get(index).add(e1);
                    Angles.get(index).add(e2);
                    Angles.get(index).add(g_bonds[j][1]);
                    index++;
                }
                else
                  if(e1 == g_bonds[j][0])
                  {
                    Angles.add(index,new ArrayList<Integer>());
                    Angles.get(index).add(e2);
                    Angles.get(index).add(e1);
                    Angles.get(index).add(g_bonds[j][1]);
                    index++;
                  }
                
                j++;
            }
        }
        
        System.out.println("\nAngles:\n");
        molecularModeling.log.append(">>> Angles: "+Angles.size()+"\n");
        for(int i = 0; i < Angles.size(); i++)
        {
            System.out.println((Angles.get(i).get(0)-1)+"-"+(Angles.get(i).get(1)-1)+"-"+(Angles.get(i).get(2)-1)+" ");
            molecularModeling.log.append("\t"+(Angles.get(i).get(0)-1)+"-"+(Angles.get(i).get(1)-1)+"-"+(Angles.get(i).get(2)-1)+"\n");
        }
        
        int[][] angle = new int[Angles.size()][3];
        for (int i = 0; i < Angles.size(); i++) {
            
                angle[i][0] = Angles.get(i).get(0);
                angle[i][1] = Angles.get(i).get(1);
                angle[i][2] = Angles.get(i).get(2);
            }
        
        return angle;
    }

    /**
     * 
     * @see molecularModeling#CalculateGeometryActionPerformed(java.awt.event.ActionEvent)
     * @see molecularModeling#WriteDBPropertiesInPmFile(java.io.BufferedWriter)
     * 
     * @param iDlinks Array of ids(index) and there links from the table
     * @return int[][] Array of calculated Dihedrals.
     * 
     */
    public int[][] getDihedralList(int[][] iDlinks)
    {
        ArrayList<ArrayList<Integer>> Diheds = new ArrayList<ArrayList<Integer>>();
        
        int length = g_bonds.length;
        int index = 0;
        for(int i = 0; i < length; i++)
        {
            int e1 = g_bonds[i][0];
            int e2 = g_bonds[i][1];
            
            int j = i + 1;
            //for(int j = i + 1; j < length; j++)
            while(j< length)
            {
                if(e2 == g_bonds[j][0] )
                {
                    for(int k= j + 1; k < length; k++)
                    {
                        if(g_bonds[j][1] == g_bonds[k][0])
                        {
                            Diheds.add(index,new ArrayList<Integer>());
                            Diheds.get(index).add(e1);
                            Diheds.get(index).add(e2); 
                            Diheds.get(index).add(g_bonds[k][0]);
                            Diheds.get(index).add(g_bonds[k][1]);

                            break;
                        }
                    }
                    
                }
                else
                if(e1 == g_bonds[j][0]){
                    
                    for(int k= j + 1; k < length; k++)
                    {
                        if(g_bonds[j][1] == g_bonds[k][0])
                        {
                            Diheds.add(index,new ArrayList<Integer>());
                            Diheds.get(index).add(e2);
                            Diheds.get(index).add(e1); 
                            Diheds.get(index).add(g_bonds[k][0]);
                            Diheds.get(index).add(g_bonds[k][1]);

                            break;
                        }
                    }
                    
                    for(int k= j + 1; k < length; k++)
                    {
                        if(e2 == g_bonds[k][0])
                        {
                            Diheds.add(index,new ArrayList<Integer>());
                            Diheds.get(index).add(g_bonds[j][1]);
                            Diheds.get(index).add(e1); 
                            Diheds.get(index).add(e2);
                            Diheds.get(index).add(g_bonds[k][1]);

                            break;
                        }
                    }
                    
                }
                
                j++;
            }
        }
        //System.out.println("\n Size of diheds Before:\t"+Diheds.size());
        System.out.println("\nDihedrals:\n");
        molecularModeling.log.append(">>>Dihedrals: "+Diheds.size()+"\n");
        for(int i = 0; i < Diheds.size(); i++)
        {
            for(int j = 0; j < Diheds.size(); j++)
            {
                if((Diheds.get(i).get(1) == Diheds.get(j).get(1)) && (Diheds.get(i).get(2) == Diheds.get(j).get(2))
                        && (Diheds.get(i).get(3) == Diheds.get(j).get(3)) && (i != j))
                {
                    if(Diheds.get(i).get(0) > Diheds.get(j).get(0))
                        Diheds.remove(i);
                    else
                        Diheds.remove(j);
                }
            } 
        }
       // System.out.println("\n Size of diheds After:\t"+Diheds.size());
        for(int i = 0; i < Diheds.size(); i++)
        {
            System.out.println((Diheds.get(i).get(0)-1)+"-"+(Diheds.get(i).get(1)-1)+"-"+(Diheds.get(i).get(2)-1)+"-"+(Diheds.get(i).get(3)-1)+" ");
            molecularModeling.log.append("\t"+(Diheds.get(i).get(0)-1)+"-"+(Diheds.get(i).get(1)-1)+"-"+(Diheds.get(i).get(2)-1)+"-"+(Diheds.get(i).get(3)-1)+"\n");
        }
        
        int[][] di = new int[Diheds.size()][4];
        for (int i = 0; i < Diheds.size(); i++) {
            
                di[i][0] = Diheds.get(i).get(0);
                di[i][1] = Diheds.get(i).get(1);
                di[i][2] = Diheds.get(i).get(2);
                di[i][3] = Diheds.get(i).get(3);
            }
        
        return di;
        
    }
}
