# Project 1: Molecular dynamics simulation of Argon atoms

Classical dynamics of $`N`$ Argon atoms inside a $`d`$-dimensional box with periodic boundary conditions. 

The atom-pair interaction is modelled as a Lennard-Jones potential and the numerical evolution can be done using the Euler method or the velocity-Verlet algorithm. The observables than can be calculated are (1) pair correlation function, (2) specific heat, (3) mean-squared displacement, (4) diffusion constant and (5) pressure. 


## Setup

Clone this repo and run `pip install -r requirements.txt` to install its dependencies.


## Usage

Open `main.py` and specify the input parameters:
- Box and particles
    - `particle_num` : number of particles, always in the shape $`4k^3`$ with $`k \in \mathbb{N}`$ to completely fill the box
    - `dim` : dimensionality of the box (FCC lattice only implemented for `dim=3`)
    - `lattice_const` : lattice constant of the FCC in units of $`\sigma`$
- Temperature
    - `temperature` : temperature in units of $`\epsilon / k_{B}`$
    - `temperature_error` : maximum error in the temperature when rescaling in units of $`\epsilon / k_{B}`$
    - `rescale_time` : time interval between rescalings in units of $`(m \sigma^2 / \epsilon )^{1/2}`$
- Molecular simulation
    - `run_time` : run time of the simulation in units of $`(m \sigma^2 / \epsilon )^{1/2}`$
    - `num_tsteps` : number of time steps of the simulation
    - `algorithm_method` : algorithm used to calculate the temporal evolution (`verlet` or `euler`)
    - `simulation` : list of steps to do `["equilibrium", "simulation"]`
- Post-processing of the simulation
    - `observables` : list of observables to calculate `["pair_correlation", "specific_heat", "displacement", "diffusion", "pressure"]`
    - `plotting` : list of plots to perform `["gif", "Evst"]`

For a more detailed information of the available plots (including e.g. histogram of velocities and plots of observables), see `plotting.py`. 


## Authors 
- Álvaro Bermejillo Seco
- Dani Bedialauneta Rodríguez
- Marc Serra Peralta
