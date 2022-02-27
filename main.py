import numpy as np
import sys

import simulate as sim
import plotting as plot

# Input parameters
N = 4*(5**3) # number of particles
d = 3 # dimensionality of the box
lattice_const = 1.54478 # lattice constant in units of sigma
T = 2 # temperature in units of kB/epsilon
num_tsteps = 900 # number of steps of the simulation
run_time = 0.9 # run time of the simulation in units of sqrt(mass*sigma^2/epsilon)
algorithm_method = "verlet" # method for numerical time evolution (options: "verlet" or "euler")
rescale_time = 0.1 # interval between rescalings
T_error = 0.05 # error in the temperature wen rescaling

# Velocity and Position initialization (uniformly random)
init_pos, L = sim.fcc_lattice(N, lattice_const)
init_vel = sim.init_velocity(N, T)

# Run simulation
print("Nparticles={}, lattice_const={}, box_dimension={}, temperature={}, run_time={}, num_tsteps={}, algorithm_method={}".format(N, lattice_const, L, T, run_time, num_tsteps, algorithm_method))

# Get system to equilibrium
print("BRINGING SYSTEM TO EQUILIBRIUM...")
eq_reached = sim.get_equilibrium(init_pos, init_vel, num_tsteps, run_time/num_tsteps, L, T, "output_eq.csv", method=algorithm_method, resc_thr=[T_error, rescale_time])
if not eq_reached:
	print("ERROR: equilibrium not reached")
	sys.exit(0)
print("DONE")

# Run actual simulation under equilibrium
print("RUNNING SIMULATION UNDER EQUILIBRIUM...")
eq_pos, eq_vel = sim.load_final_data("output_eq.csv")
print("actual temperature = {:0.5f}".format(sim.temperature(eq_vel)))
sim.simulate(eq_pos, eq_vel, num_tsteps, run_time/num_tsteps, L, T, "output.csv", method=algorithm_method)
print("DONE")

# Check energy conservation
plot.E_vs_t("output.csv", L, kinetic=True, potential=False, total=False)