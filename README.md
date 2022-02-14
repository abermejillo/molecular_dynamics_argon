# Project 1: Molecular dynamics simulation of Argon atoms

Classical dynamics of $`N`$ Argon atoms inside a square box in 2D with periodic boundary conditions. 

The atom-pair interaction is modelled as a Lennard-Jones potential and the numerical evolution is done using the Euler method. 


## Setup

Clone this repo and run `pip install -r requirements.txt` to install its dependencies.


## Usage

Open `main.py` and specify your input parameters:
- `N` : number of particles
- `L` : box length in SI \[m\]
- `T` : temperature in SI \[K\]
- `num_tsteps` : number of time steps of the simulation
- `run_time` : run time of the simulation in SI \[s\]

Run `python main.py` to (1) execture the simulation and store the results in `output.csv`, 
(2) create a gif of the evolution, and (3) plot the total energy as a function of time. 


## Authors 
- Álvaro Bermejillo
- Dani Bedialauneta
- Marc Serra Peralta
