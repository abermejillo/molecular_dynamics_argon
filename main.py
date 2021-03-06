import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

import simulate as sim
import plotting as plot
import observables as obs

##########################################################

particle_num = 4*(5**3) 
dim = 3
lattice_const = 2*1.5471 
temperature = 0.1 
temperature_error = 0.01 
rescale_time = 0.1 

run_time = 2
num_tsteps = 2000 
algorithm_method = "verlet"

simulation = ["equilibrium", "simulation"] 
observables = ["pair_correlation", "specific_heat", "pressure", "diffusion"] 
plotting = ["gif", "Evst"] 

##########################################################

box_length = (int(((particle_num - 1)/4)**(1/3)) + 1)*lattice_const
print("Nparticles={}, lattice_const={}, box_dimension={}, temperature={}, run_time={}, num_tsteps={}, algorithm_method={}".format(particle_num, lattice_const, box_length, temperature, run_time, num_tsteps, algorithm_method))

if "equilibrium" in simulation:
	print("BRINGING SYSTEM TO EQUILIBRIUM...")
	init_pos, box_length = sim.fcc_lattice(particle_num, lattice_const)
	init_vel = sim.init_velocity(particle_num, temperature)
	eq_reached = sim.get_equilibrium(init_pos, init_vel, num_tsteps, run_time/num_tsteps, box_length, temperature, "output_eq.csv", method=algorithm_method, resc_thr=[temperature_error, rescale_time])
	if not eq_reached:
		print("ERROR: equilibrium not reached")
		sys.exit(0)
	print("DONE")

if "simulation" in simulation:
	print("RUNNING SIMULATION UNDER EQUILIBRIUM...")
	eq_pos, eq_vel = sim.load_final_data("output_eq.csv")
	temperature_eq = sim.temperature(eq_vel)
	print("Temperature = {:0.5f}".format(temperature_eq))
	sim.simulate(eq_pos, eq_vel, num_tsteps, run_time/num_tsteps, box_length, "output.csv", method=algorithm_method)
	print("DONE")

if "gif" in plotting:
	plot.GIF_3D("movie_FCClattice.gif", "output.csv", 100, box_length)

if "Evst" in plotting:
	plot.E_vs_t("output_eq.csv", box_length, kinetic=True, potential=False, total=False, T=temperature, T_error=temperature_error)

if "pair_correlation" in observables:
	print("CALCULATING PAIR CORRELATION FUNCTION...")
	r, g = obs.pair_correlation_function("output.csv", 0.01, box_length, r_max=3)
	plot.plot_pair_correlation_function(r, g)
	print("The relative maxima are located in the following positions (in units of sigma)")
	print(r[argrelextrema(g, np.greater)])
	print("DONE")

if "specific_heat" in observables:
	print("CALCULATING SPECIFIC HEAT PER ATOM...")
	c, Ac_autocorr, Ac_datablock = obs.specific_heat_error("output.csv")
	eq_pos, eq_vel = sim.load_final_data("output_eq.csv")
	temperature_eq = sim.temperature(eq_vel)
	print("specific heat per atom (T = {:0.5f}) = {:0.5f} +-autocorrelation {:0.5f} or +- databloking {:0.5f}".format(temperature_eq, c, Ac_autocorr, Ac_datablock))
	print("DONE")

if "displacement" in observables:
	print("CALCULATING MEAN-SQUARED DISPLACEMENT...")
	time_steps, Ax2 = obs.mean_squared_displacement("output.csv")
	plot.plot_Ax2(time_steps, Ax2)
	print("DONE")

if "diffusion" in observables:
	print("CALCULATING DIFFUSION COEFFICIENT...")
	D, AD = obs.diffusion_error("output.csv")
	print("Diffusion coefficient = {:0.5f} +- {:0.5f}".format(D, AD))
	print("DONE")

if "pressure" in observables:
	print("CALCULATING PRESSURE")
	eq_pos, eq_vel = sim.load_final_data("output_eq.csv")
	temperature_eq = sim.temperature(eq_vel)
	P, AP_autocorr, AP_datablock = obs.pressure_error("output.csv", temperature_eq, box_length)
	print("pressure (T = {:0.5f}) = {:0.5f} +-autocorrelation {:0.5f} or +- databloking {:0.5f}".format(temperature_eq, P, AP_autocorr, AP_datablock))
	print("DONE")