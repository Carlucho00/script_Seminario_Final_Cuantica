# script_Seminario_Final_Cuantica
# HF.py
# Simple Hartree-Fock self-consistent-field routine based on
# Sec. II of this paper.
import numpy as np
# Specification of the atomic number:
Z ¼ 2
# Number of SCF iterations:
num_iter ¼ 20
# Atomic parameters from Table I. First the two orbital energies:
eps_1 ¼ Z**2/2
eps_2 ¼ Z**2/8
# and then the six Coulomb integrals:
I_1111 ¼ 5/8 * Z
I_1112 ¼ 2**12*np.sqrt(2)/27/7**4 * Z
I_1122 ¼ 16/9**3 * Z
I_1212 ¼ 17/3**4 * Z
I_1222 ¼ 2**9*np.sqrt(2)/27/5**5 * Z
I_2222 ¼ 77/2**9 * Z
# Open three files, one for the theta values calculated in each
# iteration, one for the two-electron (“total”) energy, and
# one for the associated orbital energy:
f_theta ¼ open(’theta.dat’, ’w’)
f_total_energy ¼ open(’total_energy.dat’, ’w’)
f_orbital_energy ¼ open(’orbital_energy.dat’, ’w’)
# Initialization of the HF SCF procedure:
theta ¼ 0
f_theta.write(“{} {}\n”.format(0, theta))
# Loop for carrying out the iterative SCF procedure:
for i in range(num_iter):
c_1 ¼ np.cos(theta) # Eq. (16)
c_2 ¼ np.sin(theta) # Eq. (17)
# Calculation of the two-electron energy [Eq. (13)]:
total_energy ¼ 2*(eps_1*c_1**2 þ eps_2*c_2**2)\
þ I_1111*c_1**4 þ 4*I_1112*c_1**3*c_2\
þ 2*(2*I_1122þI_1212)*c_1**2*c_2**2\
þ 4*I_1222*c_1*c_2**3 þ I_2222*c_2**4
f_total_energy.write(“{} {}\n”.format(i, total_energy))
434
Am. J. Phys., Vol. 89, No. 4, April 2021
R. Santra and M. Obermeyer
434# Calculation of the elements of the Fock matrix:
F_11 ¼ eps_1 þ I_1111*c_1**2\
þ 2*I_1112*c_1*c_2 þ I_1212*c_2**2 # Eq. (23)
F_12 ¼ I_1112*c_1**2 þ 2*I_1122*c_1*c_2\
þ I_1222*c_2**2 # Eq. (24)
F_21 ¼ F_12 # Eq. (25)
F_22 ¼ eps_2 þ I_1212*c_1**2\
þ 2*I_1222*c_1*c_2 þ I_2222*c_2**2 # Eq. (26)
# Calculation of the lower of the two roots of the characteristic
# polynomial of the Fock matrix, using the quadratic formula:
orbital_energy ¼ 0.5*(F_11 þ F_22)\
 np.sqrt(0.25*(F_11  F_22)**2 þ F_12*F_21)
f_orbital_energy.write(“{} {}\n”.format(iþ1, orbital_energy))
# For a given Fock matrix and orbital energy, Eq. (22) corresponds
# to a simple linear system for c_1 and c_2. By analytically solving
# this system, an explicit expression for the polar angle theta is
# obtained. This expression is used in the following for updating
# theta:
theta ¼ np.arctan ((orbital_energy  F_11)/F_12)
f_theta.write(“{} {}\n”.format(iþ1, theta))
f_theta.close()
f_total_energy.close()
f_orbital_energy.close()
