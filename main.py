import numpy as np

import simulate as sim
import plotting as plot

# Input parameters
N = 10 # number of particles
L = 3.5*np.sqrt(N) # box length in units of sigma
T = 300 # temperature in SI
num_tsteps = 500 # number of steps of the simulation
run_time = 1 # run time of the simulation in units of sqrt(m sigma^2/epsilon)

# Velocity and Position initialization (uniformly random)
init_pos = L*np.random.rand(N, 2)
init_vel = np.sqrt(sim.KB*T/sim.EPSILON) - 2*np.sqrt(sim.KB*T/sim.EPSILON)*np.random.rand(N, 2)
print(init_vel)
# Run simulation
sim.simulate(init_pos, init_vel, num_tsteps, run_time/num_tsteps, L, "output.csv")

# GIF for movement of particles
plot.GIF_2D("movie.gif", "output.csv", num_tsteps+1, L) # num_frames = num_tsteps+1

# Check energy conservation
plot.E_vs_t("output.csv", L)
plot.E_conservation("output.csv", L)