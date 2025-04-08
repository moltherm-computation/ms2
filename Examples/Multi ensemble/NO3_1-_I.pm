
#************************************************************************
#*              Publishing with the MolMod Database                     *
#* Every user agrees to cite the MolMod Database upon usage as follows  *
#* -------------------------------------------------------------------- *
#* Simon Stephan, Martin T. Horsch, Jadran Vrabec & Hans Hasse:         *
#* MolMod â€“ an Open Access Database of Force Fields for Molecular       *
#* Simulations of Fluids, Molecular Simulation, 45, 10, 806-814 (2019)  *
#************************************************************************
#
#
#                             Model data                               
# -------------------------------------------------------------------- 
# Substance     NO3_1-_I
# CAS-No.       14797-55-8
# Reference     H. Krienke and G. Schmeer: Hydration of Molecular Anions 
#               with Oxygen Sites â€“ a Monte Carlo Study, Zeitschrift fÃ
#               ¼r Physikalische Chemie 218, 6, (2004) 10.1524/zpch.218.6
#               .749.33456 
#               

NSiteTypes  =  2

SiteType   =  LJ126
NSites   =  4


#N
x   =  0.0
y   =  0.0
z   =  0.0
sigma   =  3.9
epsilon   =  100.644
mass   =  14.007

#O(1)
x   =  0.0
y   =  1.27
z   =  0.0
sigma   =  3.154
epsilon   =  77.99662
mass   =  15.999

#O(2)
x   =  1.099852263
y   =  -0.635
z   =  0.0
sigma   =  3.154
epsilon   =  77.99662
mass   =  15.999

#O(3)
x   =  -1.099852263
y   =  -0.635
z   =  0.0
sigma   =  3.154
epsilon   =  77.99662
mass   =  15.999

SiteType   =  Charge
NSites   =  4


#N
x   =  0.0
y   =  0.0
z   =  0.0
charge   =  +0.8604
mass   =  0.0
shielding   =  0.71092

#O(1)
x   =  0
y   =  1.27
z   =  0
charge   =  -0.6201
mass   =  0.0
shielding   =  0.71092

#O(2)
x   =  -1.099852263
y   =  -0.635
z   =  0
charge   =  -0.6201
mass   =  0.0
shielding   =  0.71092

#O(3)
x   =  1.099852263
y   =  -0.635
z   =  0
charge   =  -0.6201
mass   =  0.0
shielding   =  0.71092

NRotAxes   =   auto
