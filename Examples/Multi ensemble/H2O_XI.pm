
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
# Substance     H2O_XI
# CAS-No.       7732-18-5
# Reference     J. L. Abascal and C. Vega: A General Purpose Model for th
#               e Condensed Phases of Water: TIP4P/2005, The Journal of P
#               hysical Chemistry 123, 23, 234505 (2005) 10.1063/1.212168
#               7 
#               

NSiteTypes  =  2

SiteType   =  LJ126
NSites   =  1


#O
x   =  0.0
y   =  0.0
z   =  0.0
sigma   =  3.1589
epsilon   =  93.2
mass   =  16.00

SiteType   =  Charge
NSites   =  3


#e(1)
x   =  0.1546
y   =  0.0
z   =  0.0
charge   =  -1.1128
mass   =  0.0
shielding   =  0.0

#H[e(2)]
x   =  0.5858822766
y   =  0.7569503273
z   =  0.0
charge   =  +0.5564
mass   =  1.008
shielding   =  0.0

#H[e(3)]
x   =  0.5858822766
y   =  -0.7569503273
z   =  0.0
charge   =  +0.5564
mass   =  1.008
shielding   =  0.0

NRotAxes   =   auto
