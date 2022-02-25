import numpy as np

import simulate as sim
import plotting as plot

# Input parameters
N = 4*(5**3) # number of particles
d = 3 # dimensionality of the box
lattice_const = 1.54478 # lattice constant in units of sigma
T = 2 # temperature in units of kB/epsilon
num_tsteps = 700 # number of steps of the simulation
run_time = 0.7 # run time of the simulation in units of sqrt(mass*sigma^2/epsilon)
algorithm_method = "verlet" # method for numerical time evolution (options: "verlet" or "euler")

# Velocity and Position initialization (uniformly random)
init_pos, L = sim.fcc_lattice(N, lattice_const)
init_vel = sim.init_velocity(N, T)

# Run simulation
print("Nparticles={}, lattice_const={}, box_dimension={}, temperature={}, run_time={}, num_tsteps={}, algorithm_method={}".format(N, lattice_const, L, T, run_time, num_tsteps, algorithm_method))
sim.simulate(init_pos, init_vel, num_tsteps, run_time/num_tsteps, L, T, "output.csv", method=algorithm_method)

# Check energy conservation
plot.E_vs_t("output.csv", L, kinetic=True, potential=False, total=False)