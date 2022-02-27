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