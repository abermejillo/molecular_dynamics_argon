import numpy as np

import simulate as sim
import plotting as plot

# Input parameters
N = 10 # number of particles
L = 1E-6 # box length in SI
T = 300 # temperature in SI
num_tsteps = 100 # number of steps of the simulation
run_time = 1E-9 # run time of the simulation in SI

# Velocity and Position initialization (uniformly random)
init_pos = L*np.random.rand(N, 2)
init_vel = np.sqrt(sim.KB*T/sim.MASS) - 2*np.sqrt(sim.KB*T/sim.MASS)*np.random.rand(N, 2)

# Run simulation
sim.simulate(init_pos, init_vel, num_tsteps, run_time/num_tsteps, L, "output.csv")

# GIF for movement of particles
plot.GIF_2D("movie.gif", "output.csv", num_tsteps+1, L) # num_frames = num_tsteps+1

# Check energy conservation
plot.E_vs_t("output.csv", L)
plot.E_conservation("output.csv", L)