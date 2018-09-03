#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# Paramètres de CHARMM choisis depuis par_all36_prot.prm 
# @ http://mackerell.umaryland.edu/charmm_ff.shtml#charmm
#
# ATOMS
# MASS  -1  C         12.01100 ! carbonyl C, peptide backbone
#
# BONDS
# !
# !V(bond) = Kb(b - L0)**2
# !
# !Kb: kcal/mole/A**2
# !L0: A
# !#!atom type Kb          b0
# !
# C    C     600.000     1.3350 ! ALLOW ARO HEM
# ! Heme vinyl substituent (KK, from propene 


import math as m

# Constantes pour le calcul de l'Énergie potentielle et des Forces:
M = 12.011e-3 # en kg/mol (1 amu = 1 g/mol)
L0 = 1.3350 # en Å
K = 600.000 # en kcal/mol/Å^2

# Constantes pour la MD:
KB = 1.9872041e-3 # Constante de Boltzmann en kcal/mol/K
NDDL = 3 # 3 degrés de liberté (xB, xC et yC)
NSTEPS = 1000 # Nombre de pas
DELTA_t = 0.0005 # en ps
DELTA = 0.00001

def calc_Ep(xB, xC, yC):
	# Champs de forces empiriques contenant 3 oscillateurs harmoniques (Énergie potentielle en kcal/mol)
    Ep = 0.5 * K * ( ( xB - L0 )**2 + ( m.sqrt( xC**2 + yC**2 ) - L0 )**2 + ( m.sqrt( (xC - xB)**2 + yC**2 ) - L0 )**2 )
    return(Ep)

def init_v():
    # Vitesse initiale en Å/ps pour T = 300 K
    return(( ( (300*KB) / M )*0.418 )**.5)

def do_MD(xB, xC, yC):

	file = open("MD_cyclopropane.dat","w")
	header = "%4s %6s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n" % \
	("step", "t", "T", "xB" , "xC", "yC", "xG", "yG", "F_xB", "F_xC", "F_yC", "E_pot", "E_cin", "E_tot")
	file.write(header)

	for step in range(NSTEPS):
		t = step * DELTA_t
		if step == 0:
			# Détermination des coordonnées précédentes en Å (avant t = 0)
            # (Expansion de Taylor tronquée après le premier terme):
			xB_prev = xB - init_v() * DELTA_t
			xC_prev = xC - init_v() * DELTA_t
			yC_prev = yC - init_v() * DELTA_t

		# Forces en kcal/mol/Å (dérivées partielles analytiques de l'Énergie potentielle) :
		F_xB = - ( (calc_Ep(xB + DELTA, xC, yC) - calc_Ep(xB - DELTA, xC, yC) ) / (2 * DELTA)) 
		F_xC = - ( (calc_Ep(xB, xC + DELTA, yC) - calc_Ep(xB, xC - DELTA, yC) ) / (2 * DELTA))
		F_yC = - ( (calc_Ep(xB, xC, yC + DELTA) - calc_Ep(xB, xC, yC - DELTA) ) / (2 * DELTA))

		# Accélérations en Å/ps^-2 (2ème loi de Newton):
		acc_xB = (F_xB / M) * 0.418
		acc_xC = (F_xC / M) * 0.418
		acc_yC = (F_yC / M) * 0.418

		# Nouvelles coordonnées en Å (Algorithme de Verlet):
		xB_new = (2 * xB) - xB_prev + (DELTA_t**2 * acc_xB)
		xC_new = (2 * xC) - xC_prev + (DELTA_t**2 * acc_xC)
		yC_new = (2 * yC) - yC_prev + (DELTA_t**2 * acc_yC)

		# Vitesses en Å/ps (Algorithme de Verlet):
		v_xB = ( xB_new - xB_prev ) / (2 * DELTA_t)
		v_xC = ( xC_new - xC_prev ) / (2 * DELTA_t)
		v_yC = ( yC_new - yC_prev ) / (2 * DELTA_t)

		# Énergie cinétique en kcal/mol: 
		E_cin = .5 * M * ( (v_xB*100.0)**2 + (v_xC*100.0)**2 + (v_yC*100.0)**2  )
		E_cin /= (1000.0 * 4.18)

		# Énergie potentielle en kcal/mol:
		E_pot = calc_Ep(xB, xC, yC)

		# Énergie totale en kcal/mol:
		E_tot = E_pot + E_cin

		# Température en K:
		T = (2 * E_cin) / (NDDL * KB)

		# Coordonnées du barycentre en Å:
		xG = (xB + xC)/3
		yG = yC / 3

		out = "%4i %6.4f %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f\n" % \
		(step, t, T, xB , xC, yC, xG, yG, F_xB, F_xC, F_yC, E_pot, E_cin, E_tot)
		file.write(out)

		xB_prev, xC_prev, yC_prev = xB, xC, yC
		xB, xC, yC = xB_new, xC_new, yC_new

	file.close()

if __name__ == '__main__':
	# Initialisation des coordonnées (en Å) selon un triangle équilatéral:
	xB = L0
	xC = L0 / 2
	yC = (m.sqrt(3) / 2) * L0
	do_MD(xB, xC, yC)