import numpy as np

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
N = 3 # number of particles
L = 100*1E-10 # box length
T = 300 # temperature

# Global constants
KB = 1.3806E-23 # Boltzmann constant 
SIGMA = 3.405E-10 # parameter of LJ potential
EPSILON = 119.8*KB # parameter of LJ potential
MASS = 6.6335209E-26 # Argon particle mass

# Velocity and Position initialization (uniformly random)
pos = L*np.random.rand(N, 2)
vel = np.sqrt(KB*T/MASS)*np.random.rand(N, 2)

# Compute relative positions and distances
rel_pos = np.zeros((N,N,2))
rel_dist = np.zeros((N,N))

for i in range(N):
    for j in range(i):
        rel_pos[i,j] = pos[i] - pos[j]
        rel_dist[i,j] = np.linalg.norm(rel_pos[i,j])

# Force calculation using the Lennard-Jones potential       
force = np.zeros((N,N,2)) # Forces between pairs

for i in range(N):
    for j in range(i):
        force[i,j] = 24*EPSILON*rel_pos[i,j] * (2*SIGMA**12/rel_dist[i,j]**14 - SIGMA**6/rel_dist[i,j]**8)

force = force + force.transpose((1,0,2)) # Fill upper triangle of the matrix 

# Total force exerted on each particle
total_force = force.sum(1) 