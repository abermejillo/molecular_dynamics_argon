"""
CALCULATION OF THE TIMING OF EACH FUNCTION USED IN THE SIMULATION (VERLET ALGORITHM)
====================================================================================
1. Launch a Python shell
2. Launch IPython shell (for the %timeit)
>>> from IPython import embed
>>> embed()
3. Copy the code below to the IPython shell and press Enter
4. The result looks like this:
```
N=200
atomic_distances
4.83 ms ± 39.2 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)
lj_force
3.82 ms ± 55 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)
update position
4.43 µs ± 42.1 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)
atomic_distances
3.93 ms ± 30.2 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)
lj_force
3.77 ms ± 15 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)
update velocity
4.15 µs ± 12 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)
update periodic BC
5.4 µs ± 6 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)

TOTAL SIMULATION OF ONE STEP
20.3 ms ± 4.61 ms per loop (mean ± std. dev. of 7 runs, 10 loops each)
```
"""

import simulate as sim
import numpy as np

# Input parameters
N = 200 # number of particles
d = 3 # dimensionality of the box
L = 3 * N**(1/d) # box length in units of sigma
T = 300 # temperature in SI
timestep = 1E-4

pos = L*np.random.rand(N, d)
vel = np.sqrt(sim.KB*T/sim.EPSILON) - 2*np.sqrt(sim.KB*T/sim.EPSILON)*np.random.rand(N, d)

print("N={}".format(N))

rel_pos, rel_dist = sim.atomic_distances(pos, L)
print("atomic_distances")
%timeit sim.atomic_distances(pos, L)

total_force = sim.lj_force(rel_pos, rel_dist)
print("lj_force")
%timeit sim.lj_force(rel_pos, rel_dist)

pos = pos + vel*timestep + 0.5*timestep**2*total_force
print("update position")
%timeit pos + vel*timestep + 0.5*timestep**2*total_force

rel_pos, rel_dist = sim.atomic_distances(pos, L)
print("atomic_distances")
%timeit sim.atomic_distances(pos, L)

total_force_next = sim.lj_force(rel_pos, rel_dist)
print("lj_force")
%timeit sim.lj_force(rel_pos, rel_dist)

vel = vel + 0.5*(total_force + total_force_next)*timestep
print("update velocity")
%timeit vel + 0.5*(total_force + total_force_next)*timestep

pos = pos - np.floor(pos/L)*L
print("update periodic BC")
%timeit pos - np.floor(pos/L)*L

print("\nTOTAL SIMULATION OF ONE STEP")
%timeit sim.simulate_step_verlet(pos, vel, timestep, L)