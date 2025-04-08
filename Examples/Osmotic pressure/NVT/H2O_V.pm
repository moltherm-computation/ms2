
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
# Substance     H2O_V
# CAS-No.       7732-18-5
# Reference     H. J. Berendsen, J. R. Grigera and T. P. Straatsma: The M
#               issing Term in Effective Pair Potentials, The Journal of 
#               Physical Chemistry 91, 24, 6269-6271 (1987) 10.1021/j1003
#               08a038 
#               

NSiteTypes  =  2

SiteType   =  LJ126
NSites   =  1


#O
x   =  0.0
y   =  0.0
z   =  0.0
sigma   =  3.165557883763
epsilon   =  78.197418511978
mass   =  16.00

SiteType   =  Charge
NSites   =  3


#e(1)
x   =  0.0
y   =  0.0
z   =  0.0
charge   =  -0.8476
mass   =  0.0
shielding   =  0.0

#H[e(2)]
x   =  0.816490431
y   =  0.577358967
z   =  0.0
charge   =  +0.4238
mass   =  1.008
shielding   =  0.0

#H[e(3)]
x   =  -0.816490431
y   =  0.577358967
z   =  0.0
charge   =  +0.4238
mass   =  1.008
shielding   =  0.0

NRotAxes   =   auto
