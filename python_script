# HF.py
# Simple Hartree-Fock self-consistent-field routine based on
# Sec. II of this paper.
import numpy as np
# Specification of the atomic number:
Z = 2
# Number of SCF iterations:
num_iter = 20
# Atomic parameters from Table I. First the two orbital energies:
eps_1 = -Z**2/2
eps_2 = -Z**2/8
# and then the six Coulomb integrals:
I_1111 = 5/8 * Z
I_1112 = 2**12*np.sqrt(2)/27/7**4 * Z
I_1122 = 16/9**3 * Z
I_1212 = 17/3**4 * Z
I_1222 = 2**9*np.sqrt(2)/27/5**5 * Z
I_2222 = 77/2**9 * Z
# Open three files, one for the theta values calculated in each
# iteration, one for the two-electron (“total”) energy, and
# one for the associated orbital energy:
f_theta = open("theta.dat", "w")
f_total_energy = open("total_energy.dat", "w")
f_orbital_energy = open("orbital_energy.dat", "w")
# Initialization of the HF SCF procedure:
theta = 0
f_theta.write("{} {}\n".format(0, theta))
# Loop for carrying out the iterative SCF procedure:
for i in range(num_iter):
	c_1 = np.cos(theta) # Eq. (16)
	c_2 = np.sin(theta) # Eq. (17)
# Calculation of the two-electron energy [Eq. (13)]:
	total_energy = 2*(eps_1*c_1**2 + eps_2*c_2**2)\
		+ I_1111*c_1**4 + 4*I_1112*c_1**3*c_2\
		+ 2*(2*I_1122+I_1212)*c_1**2*c_2**2\
		+ 4*I_1222*c_1*c_2**3 + I_2222*c_2**4

	f_total_energy.write("{} {}\n".format(i, total_energy))

	F_11 = eps_1 + I_1111*c_1**2\
			+ 2*I_1112*c_1*c_2 + I_1212*c_2**2 # Eq. (23)
	F_12 = I_1112*c_1**2 + 2*I_1122*c_1*c_2\
			+ I_1222*c_2**2 # Eq. (24)
	F_21 = F_12 # Eq. (25)
	F_22 = eps_2 + I_1212*c_1**2\
			+ 2*I_1222*c_1*c_2 + I_2222*c_2**2 # Eq. (26)
	# Calculation of the lower of the two roots of the characteristic
	# polynomial of the Fock matrix, using the quadratic formula:
	orbital_energy = 0.5*(F_11 + F_22)\
						- np.sqrt(0.25*(F_11 - F_22)**2 + F_12*F_21)
	f_orbital_energy.write("{} {}\n".format(i+1, orbital_energy))
	# For a given Fock matrix and orbital energy, Eq. (22) corresponds
	# to a simple linear system for c_1 and c_2. By analytically solving
	# this system, an explicit expression for the polar angle theta is
	# obtained. This expression is used in the following for updating
	# theta:
	theta = np.arctan ((orbital_energy - F_11)/F_12)
	f_theta.write("{} {}\n".format(i+1, theta))
f_theta.close()
f_total_energy.close()
f_orbital_energy.close()
