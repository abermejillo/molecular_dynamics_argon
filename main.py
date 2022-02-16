import numpy as np

import simulate as sim
import plotting as plot

# Input parameters
N = 10 # number of particles
d = 3 # dimensionality of the box
L = 3 * N**(1/d) # box length in units of sigma
T = 300 # temperature in SI
num_tsteps = 1000 # number of steps of the simulation
run_time = 1 # run time of the simulation in units of sqrt(mass*sigma^2/epsilon)
algorithm_method = "verlet" # method for numerical time evolution (options: "verlet" or "euler")

# Velocity and Position initialization (uniformly random)
init_pos = L*np.random.rand(N, d)
init_vel = np.sqrt(sim.KB*T/sim.EPSILON) - 2*np.sqrt(sim.KB*T/sim.EPSILON)*np.random.rand(N, d)

# Run simulation
sim.simulate(init_pos, init_vel, num_tsteps, run_time/num_tsteps, L, "output.csv", method=algorithm_method)

# GIF for movement of particles
plot.GIF_2D("movie.gif", "output.csv", 300, L) # all timesteps plotted with num_frames = num_tsteps+1

# Check energy conservation
plot.E_vs_t("output.csv", L)
plot.E_conservation("output.csv", L)