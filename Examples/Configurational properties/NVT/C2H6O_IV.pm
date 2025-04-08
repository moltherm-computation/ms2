
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
# Substance     C2H6O_IV
# CAS-No.       64-17-5
# Reference     T. Schnabel, J. Vrabec and H. Hasse: Henry’s Law Consta
#               nts of Methane, Nitrogen, Oxygen and Carbon Dioxide in Et
#               hanol from 273 to 498 K: Prediction from Molecular Simula
#               tion, Fluid Phase Equilibria 233, 2, 134-143 (2005) 10.10
#               16/j.fluid.2005.04.016 
#               

NSiteTypes  =  2

SiteType   =  LJ126
NSites   =  3


#CH3
x   =  -1.470767
y   =  -0.338351
z   =  0.0
sigma   =  3.6072
epsilon   =  120.15
mass   =  15.035

#CH2
x   =  +0.092772
y   =  +0.883285
z   =  0.0
sigma   =  3.4612
epsilon   =  86.291
mass   =  14.027

#O
x   =  +1.171546
y   =  -0.450976
z   =  0.0
sigma   =  3.149559
epsilon   =  85.053449
mass   =  15.999

SiteType   =  Charge
NSites   =  3


#e (1)
x   =  +0.092772
y   =  +0.883285
z   =  0.0
charge   =  +0.2556
mass   =  0.0
shielding   =  1.38448

#e (2)
x   =  +1.171546
y   =  -0.450976
z   =  0.0
charge   =  -0.697107
mass   =  0.0
shielding   =  1.25982

#H [e(3)]
x   =  +2.049156
y   =  -0.085872
z   =  0.0
charge   =  +0.441507
mass   =  1.008
shielding   =  0.0

NRotAxes   =   auto
