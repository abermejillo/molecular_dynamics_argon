import numpy as np
import matplotlib.pyplot as plt
import time

"""
Initial attempt: 3 Ar particles in 2D
=====================================
STEPS:
1. Initialize uniformly random velocity and position 
(it should be fcc for positions and Maxwell distribution for velocities)
2. Force calculation using the Lennard-Jones potential
3. Euler method for time evolution
4. Periodic boundary condition
5. Check momenta and energy conservation
(another check can be time reversibility)

UNITS: in International System
"""

# Input parameters
N = 5 # number of particles
L = 1E-6 # box length
T = 300 # temperature
At = 1E-11 # time-step
run_time = 10E-9 # run time of the simulation

# Global constants
KB = 1.3806E-23 # Boltzmann constant 
#SIGMA = 3.405E-10 # parameter of LJ potential
SIGMA = 0.5E-7
EPSILON = 119.8*KB # parameter of LJ potential
MASS = 6.6335209E-26 # Argon particle mass

# Functions to compute the energy of the system
def potential_energy(distance): # Potential energy function
    return 4*EPSILON*(SIGMA**12/distance**12-SIGMA**6/distance**6)

def total_energy(rel_dist, vel): # 
    kin_energy = np.sum(0.5*MASS*np.linalg.norm(vel)**2)
    mask = np.ones(rel_dist.shape,dtype=bool)
    mask[np.diag_indices(N)] = False 
    pot_energy = np.sum(potential_energy(rel_dist),where=mask)/2 # The diagonal terms are skipped (no self-interaction)
    return kin_energy + pot_energy

# Velocity and Position initialization (uniformly random)
pos = L*np.random.rand(N, 2)
vel = np.sqrt(KB*T/MASS) - 2*np.sqrt(KB*T/MASS)*np.random.rand(N, 2)

# Total energy storage vector
tot_energy=[]

# Euler method
for k, t in enumerate(np.arange(0, run_time + At, At)):

    # Compute relative positions
    pos1 = np.repeat(pos[:, np.newaxis, :], N, axis=1)
    pos2 = np.repeat(pos[np.newaxis, :, :], N, axis=0)
    rel_pos = pos1 - pos2

    rel_pos_ = rel_pos.copy() # for plotting

    # minimum image convention (for the periodic boundary conditions)
    wrong_pairs = np.where(np.abs(rel_pos) > L/2)
    rel_pos[wrong_pairs] = rel_pos[wrong_pairs] - np.sign(rel_pos[wrong_pairs])*L

    # Compute relative distance
    rel_dist = np.linalg.norm(rel_pos, axis=2) # axis 2 contains the cartesian coordinates
    rel_dist = rel_dist[:,:,np.newaxis] # add axis for LJ force calculation
    rel_dist[np.diag_indices(N)] = 1 # avoiding division by zero in the diagonal when calculating LJ force
    # Force calculation using the Lennard-Jones potential  
    force = 24*EPSILON*rel_pos*(2*SIGMA**12/rel_dist**14 - SIGMA**6/rel_dist**8)

    # Total force exerted on each particle
    total_force = np.sum(force,axis=0)

    # Update velocities and positions
    pos = pos + vel*At
    vel = vel + total_force*At/MASS

    # Check periodic boundary conditions
    pos = pos - np.floor(pos/L)*L

    # Save the total energy
    tot_energy += [total_energy(rel_dist, vel)]

    # plotting for checking the interactions
    if True:
        for i in range(N):
            # plot normal box and its eight neighbours
            plt.plot(pos[i,0], pos[i,1], "r.")
            plt.plot(pos[i,0]+L, pos[i,1], "r.")
            plt.plot(pos[i,0], pos[i,1]+L, "r.")
            plt.plot(pos[i,0]+L, pos[i,1]+L, "r.")
            plt.plot(pos[i,0]-L, pos[i,1], "r.")
            plt.plot(pos[i,0], pos[i,1]-L, "r.")
            plt.plot(pos[i,0]-L, pos[i,1]-L, "r.")
            plt.plot(pos[i,0]-L, pos[i,1]+L, "r.")
            plt.plot(pos[i,0]+L, pos[i,1]-L, "r.")

            for j in range(i):
                plt.plot([pos[j,0], pos[j,0] + rel_pos_[i,j,0]], [pos[j,1], pos[j,1] + rel_pos_[i,j,1]], "y-") # no BC correction
                plt.plot([pos[j,0], pos[j,0] + rel_pos[i,j,0]], [pos[j,1], pos[j,1] + rel_pos[i,j,1]], "b--") # correction for BC

        plt.plot([0,0,L,L,0],[0,L,L,0,0], "g-") # plot square for normal box
        plt.xlim(-L/2, 3*L/2)
        plt.ylim(-L/2, 3*L/2)
        plt.savefig("{}.png".format(k))
        plt.cla()
plt.close()

# Plotting the energy

plt.figure(2)
plt.plot(np.linspace(0, run_time, num = int(run_time/At)+2),tot_energy)
plt.show()
