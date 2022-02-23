# Project 1: Molecular dynamics simulation of Argon atoms

Classical dynamics of $`N`$ Argon atoms inside a $`d`$-dimensional box with periodic boundary conditions. 

The atom-pair interaction is modelled as a Lennard-Jones potential and the numerical evolution can be done using the Euler method or the velocity-Verlet algorithm. 


## Setup

Clone this repo and run `pip install -r requirements.txt` to install its dependencies.


## Usage

Open `main.py` and specify the input parameters:
- `N` : number of particles, always in the shape 4*(k**3)
- `d` : dimensionality of the box
- `T` : temperature in units of $`k_{B}/\epsilon`$
- `num_tsteps` : number of time steps of the simulation
- `run_time` : run time of the simulation in units of $`(m \sigma^2 / \epsilon )^{1/2}`$
- `algorithm_method` : algorithm to calculate the temporal evolution (`verlet` or `euler`)



Run `main.py` to:
1. Execute the simulation and store the results in `output.csv`:

    `sim.simulate(init_pos, init_vel, num_tsteps, run_time/num_tsteps, L, "output.csv", method=algorithm_method)`

2. Create a 3D gif of the evolution:  

    `plot.GIF_3D("movie.gif", "output.csv", 300, L)`

3. Plot the total energy as a function of time (as well as kinetic and potential energies separately):

    `plot.E_vs_t("output.csv", L, kinetic_potential=True)`

4. Plot the relative energy deviation as a function of time:

    `plot.E_conservation("output.csv", L)`

5. Plot the relative distance between two particles as a function of time:

    `plot.reldist_vs_t("output.csv", 0, 1, L)`

6. Create a GIF that shows how the energy is translated from kinetic to potential and viceversa (only significant with two particles). 

    `plot.GIF_potential_energy("movie2.gif", "output.csv", 300, 0 , 1, L)`

7. Plot the probability density function for the initial velocities (Maxwell distribution). 

    `plot.plot_maxwell_distribution(init_vel,T)`



## Authors 
- Álvaro Bermejillo Seco
- Dani Bedialauneta Rodríguez
- Marc Serra Peralta
