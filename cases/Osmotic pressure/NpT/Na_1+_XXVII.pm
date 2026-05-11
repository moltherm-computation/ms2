
#************************************************************************
#*              Publishing with the MolMod Database                     *
#* Every user agrees to cite the MolMod Database upon usage as follows  *
#* -------------------------------------------------------------------- *
#* Simon Stephan, Martin T. Horsch, Jadran Vrabec & Hans Hasse:         *
#* MolMod – an Open Access Database of Force Fields for Molecular       *
#* Simulations of Fluids, Molecular Simulation, 45, 10, 806-814 (2019)  *
#************************************************************************
#
#
#                             Model data                               
# -------------------------------------------------------------------- 
# Substance     Na_1+_XXVII
# CAS-No.       17341-25-2
# Reference     Y. Qiu, Y. Jiang, Y. Zhang and H. Zhang: Rational Design 
#               of Nonbonded Point Charge Models for Monovalent Ions with
#                Lennard-Jones 12-6 Potential, The Journal of Physical Ch
#               emistry B 125, 49, 13502-13518 (2021) 10.1021/acs.jpcb.1c
#               09103 
#               

NSiteTypes  =  2

SiteType   =  LJ126
NSites   =  1


#Na_1+
x   =  0.0
y   =  0.0
z   =  0.0
sigma   =  1.39506
epsilon   =  7.76
mass   =  22.9892

SiteType   =  Charge
NSites   =  1


#e(1)
x   =  0.0
y   =  0.0
z   =  0.0
charge   =  +1.0
mass   =  0.0
shielding   =  1.0

NRotAxes   =   auto
