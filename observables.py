import numpy as np

import simulate as sim


def pair_correlation_function(file_name, dr, box_length, r_max=None):
    """
    Returns the pair correlation function averaged over time. 

    Parameters
    ----------
    file_name : str
        Name of the CSV file in which the data is stored
    dr : float
		distance between bins in the histogram
	box_length : float
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

    if r_max is None: r_max = box_length
    r = np.arange(dr, r_max+dr, dr) # start != 0, because g = 1/0 is not defined
    n = np.zeros(len(r))
    time, pos, _ = sim.load_data(file_name)
    particle_num = pos.shape[0]

    for k, t in enumerate(time):
        print("\r{}/{}".format(k+1, len(time)), end="")
        rel_pos, rel_dist = sim.atomic_distances(pos[k], box_length)
        for i, r_i in enumerate(r):
            n[i] += len(np.where((rel_dist >= r_i) & (rel_dist < r_i + dr))[0])
    
    g = 2*box_length**3 / (particle_num*(particle_num-1)) * (n/len(time)) / (4*np.pi*r**2 * dr)

    print("")

    return r, g


def specific_heat(file_name, starting_time_step=0):
    """
    Computes the specific heat per atom of a system.

    Parameters
    ----------
    file_name : str
        Name of the CSV file in which the data is stored
    starting_time_step : int
        Number of time steps to skip from the beginning of the file_name

    Returns
    -------
    c : float
        Specific heat per atom of the system
    """

    time, pos, vel = sim.load_data(file_name)
    time = time[starting_time_step:]
    pos = pos[starting_time_step:,:,:]
    vel = vel[starting_time_step:,:,:]
    num_tsteps = len(time) 
    particle_num = np.shape(vel)[1] 
    square_of_mean_kin = (sim.kinetic_energy(vel)/num_tsteps)**2
    mean_of_square_kin = (((0.5*(vel**2).sum(2)).sum(1))**2).sum()/num_tsteps
    
    r = (mean_of_square_kin - square_of_mean_kin)/square_of_mean_kin # relative fluctuations in kinetic energy

    c = 1.5/(1-3*particle_num*r/2)

    return c