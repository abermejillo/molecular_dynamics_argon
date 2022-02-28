import numpy as np

import simulate as sim


def pair_correlation_function(file_name, dr, L, r_max=None):
    """
    Returns the pair correlation function averaged over time. 

    Parameters
    ----------
    file_name : str
        Name of the CSV file in which the data is stored
    dr : float
		distance between bins in the histogram
	L : float
		Box size
    r_max : float
        Maximum value for r (distance between particles)

    Returns
    -------
    r : np.ndarray(int(r_max/dr)+1)
        Distance between pairs of particles
    g : np.ndarray(int(r_max/dr)+1)
        Pair correlation function as a function of r
    """

    if r_max is None: r_max = L
    r = np.arange(dr, r_max+dr, dr) # start != 0, because g = 1/0 is not defined
    n = np.zeros(len(r))
    time, pos, _ = sim.load_data(file_name)
    N = pos.shape[0]

    for k, t in enumerate(time):
        print("\r{}/{}".format(k+1, len(time)), end="")
        rel_pos, rel_dist = sim.atomic_distances(pos[k], L)
        for i, r_ in enumerate(r):
            n[i] += len(np.where((rel_dist >= r_) & (rel_dist < r_ + dr))[0])
    
    g = 2*L**3 / (N*(N-1)) * (n/len(time)) / (4*np.pi*r_**2 * dr)

    print("")

    return r, g

def specific_heat(data_file,starting_time_step):
    """
    Computes the specific heat per atom of a system.

    Parameters
    ----------
    data_file : str
        Name of the CSV file in which the data is stored

    Returns
    -------
    c : float
        Specific heat per atom of the system.
    """

    time, pos, vel = sim.load_data(data_file)
    time = time[starting_time_step:-1]
    pos = pos[starting_time_step:-1,:,:]
    vel = vel[starting_time_step:-1,:,:]
    M = len(time) # Number of time steps
    N = np.shape(vel)[1] # Number of particles
    square_of_mean_kin = (sim.kinetic_energy(vel)/M)**2
    mean_of_square_kin = (((0.5*(vel**2).sum(2)).sum(1))**2).sum()/M
    
    r = (mean_of_square_kin - square_of_mean_kin)/square_of_mean_kin # relative fluctuations in kinetic energy

    c = 1.5/(1-3*N*r/2)

    return c