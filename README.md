# Project 1: Molecular dynamics simulation of Argon atoms

Classical dynamics of $`N`$ Argon atoms inside a $`d`$-dimensional box with periodic boundary conditions. 

The atom-pair interaction is modelled as a Lennard-Jones potential and the numerical evolution can be done using the Euler method or the velocity-Verlet algorithm. 


## Setup

Clone this repo and run `pip install -r requirements.txt` to install its dependencies.


## Usage

Open `main.py` and specify the input parameters:
- `N` : number of particles
- `L` : box length in units of $`\sigma`$
- `d` : dimensionality of the box
- `T` : temperature in SI \[K\]
- `num_tsteps` : number of time steps of the simulation
- `run_time` : run time of the simulation in units of $`(m \sigma^2 / \epsilon )^{1/2}`$
- `algorithm_method` : algorithm to calculate the temporal evolution (`verlet` or `euler`)

Run `main.py` to:
1. execture the simulation and store the results in `output.csv`
2. create a gif of the evolution
3. plot the total energy as a function of time (as well as kinetic and potential energies separately)
4. relative energy deviation as a function of time
5. relative distance between two particles as a funciton of time 


## Authors 
- Álvaro Bermejillo
- Dani Bedialauneta
- Marc Serra Peralta
