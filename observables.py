import numpy as np
from scipy import optimize

import simulate as sim

#---------------------------------------------------------
# FIRST FUNCTION BLOCK: OBSERVABLES WITHOUT ERROR
#---------------------------------------------------------

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


def specific_heat(file_name):
    """
    Computes the specific heat per atom of a system.

    Parameters
    ----------
    file_name : str
        Name of the CSV file in which the data is stored

    Returns
    -------
    c : float
        Specific heat per atom of the system
    """

    time, pos, vel = sim.load_data(file_name)
    num_tsteps = len(time) 
    particle_num = np.shape(vel)[1] 

    total_kin = (0.5*(vel**2).sum(2)).sum(1)
    mean_kin = total_kin.sum()/num_tsteps
    square_of_mean_kin = mean_kin**2
    mean_of_square_kin = (total_kin**2).sum()/num_tsteps
    mean_of_fourth_kin = (total_kin**4).sum()/num_tsteps
    r = mean_of_square_kin/square_of_mean_kin - 1 # relative fluctuations in kinetic energy

    c = 1.5/(1-3*particle_num*r/2)

    return c


def diffusion(file_name):
    """
    Computes the diffussion constant D. 

    Parameters
    ----------
    file_name : str
        Name of the CSV file in which the data is stored

    Returns
    -------
    D : float
        Diffussion constant
    """

    time, pos, _ = sim.load_data(file_name)
    particle_num = pos.shape[0]

    dist = (pos[-1] - pos[0])
    dist = (dist*dist).sum(axis=1)**0.5
    D = (dist.sum()/particle_num)/(6*time[-1])

    return D


def pressure(file_name, T, box_length):
    """
    Computes the pressure of the system.

    Parameters
    ----------
    file_name : str
        Name of the CSV file in which the data is stored
    T : float
        Temperature
    box_length : float
		Box size

    Returns
    -------
    P : float
        Pressure
    """
    time, pos, vel = sim.load_data(file_name)

    N = np.shape(pos)[1]
    M = len(time)

    second_term_instantenous = np.zeros(M)

    for k, t in enumerate(time):
        print("\r{}/{}".format(k+1, len(time)), end="")

        rel_pos, rel_dist = sim.atomic_distances(pos[k], box_length)

        rel_dist = rel_dist[:,:,np.newaxis] # add axis for LJ force calculation (so that it agrees with rel_pos dimensions)
        rel_dist[np.diag_indices(np.shape(rel_dist)[0])] = 1 # avoiding division by zero in the diagonal when calculating LJ force

        matrix = (1/(6*N*T))*24*(2/rel_dist**12-1/rel_dist**7)
        matrix[np.diag_indices(np.shape(matrix)[0])] = 0 # diagonal terms should be zero by definition

        second_term_instantenous[k] = matrix.sum()

    print("")

    BP_rho =  1+np.average(second_term_instantenous)

    P = BP_rho*T*N/box_length**3

    return P


def mean_squared_displacement(file_name, time_steps=None):
    """
    Returns the mean-squared displacement as a function of time. 

    Parameters
    ----------
    file_name : str
        Name of the CSV file in which the data is stored
    time_steps : np.ndarray
        Array of times in which to calculate the mean-squared displacement
        If None, it uses time from sim.load_data as time_steps

    Returns
    -------
    time_steps : np.ndarray
        See Parameters
    Ax2 : np.ndarray(len(time))
        Mean-squared displacement as a function of time
    """

    time, pos, _ = sim.load_data(file_name)
    particle_num = pos.shape[0]
    if time_steps is None: time_steps = time
    Ax2 = np.zeros(len(time_steps))

    for k, t in enumerate(time_steps):
        dist = (pos[k] - pos[0])
        dist = (dist*dist).sum(axis=1)**0.5
        Ax2[k] = dist.sum()/particle_num

    return time_steps, Ax2


#---------------------------------------------------------
# SECOND FUNCTION BLOCK: AUXILIARY FUNCTIONS FOR OBSERVABLES AND ERROR COMPUTATION
#---------------------------------------------------------

def autocorrelation_function(data):
    """
    Returns the autocorrelation function for a given variable as a function of time. 

    Parameters
    ----------
    data : np.ndarray(len(time))
        Variable as a function of time

    Returns
    -------
    Xa : np.ndarray(len(time)-1)
        Autocorrelation function for the given variable as a function of time
    """

    N = len(data)
    Xa = np.zeros(N-1)

    for t in range(N-1):
        n = N - t
        Xa[t] = ((N-t)*(data[:n]*data[t:n+t]).sum() - (data[:n]).sum()*(data[t:n+t]).sum()) /  \
                (np.sqrt((N-t)*(data[:n]**2).sum() - (data[:n].sum())**2) * np.sqrt((N-t)*(data[t:n+t]**2).sum() - (data[t:n+t].sum())**2))

    return Xa


def correlation_time(data, runtime, num_tsteps):
    """
    Returns the autocorrelation time for a given data set evenly distributed along the simulation runtime. 

    Parameters
    ----------
    data : np.ndarray(len(time))
        Variable as a function of time
    runtime : float
        Runtime of the simulation
    num_tsteps : float
        Number of timesteps in which the simulation is divided.

    Returns
    -------
    tau : float
        Autocorrelation time for the given variable.
    """

    Xa = autocorrelation_function(data)
    time = np.arange(num_tsteps)
    N_max = np.where(np.diff(Xa) > 0)[0][0] # find first point that Xa increases
    function = lambda x,tau: np.exp(-x/tau) # function for fitting the error
    popt, _ = optimize.curve_fit(function,  time[:N_max],  Xa[:N_max])
    tau = popt[0]

    return tau


def data_blocking(data, b_range):
    """
    Returns error from data blocking as a function of block length size (b_range) and also the total one. 

    Parameters
    ----------
    data : np.ndarray(len(time))
        Variable as a function of time
    b_range : np.ndarray
        Block lenght sizes

    Returns
    -------
    error : float
        Error of the given variable
    sigma : np.ndarray(len(b_range))
        Error of the given variable as a function of block length size
    """

    N = len(data)
    sigma = np.zeros(len(b_range))

    for k, b in enumerate(b_range):
        Nb = int(N/b) 
        data_ = data[:b*Nb]
        blocks = data_.reshape(Nb, b)
        average_blocks = np.average(blocks, axis=1)
        sigma[k] = np.sqrt(((average_blocks**2).sum()/Nb - (average_blocks.sum()/Nb)**2) / (Nb-1)) 

    try: 
        function = lambda x,a,b,c: b-c*np.exp(-a*x) # function for fitting the error
        popt, _ = optimize.curve_fit(function,  b_range,  sigma)
        error = popt[1]
    except:
        error = np.max(sigma)

    return error, sigma


#---------------------------------------------------------
# THIRD FUNCTION BLOCK: OBSERVABLES WITH ERROR COMPUTATION
#---------------------------------------------------------

def specific_heat_error(file_name):
    """
    Computes the specific heat per atom of a system. Gives the error with autocorrelation function method.

    Parameters
    ----------
    file_name : str
        Name of the CSV file in which the data is stored

    Returns
    -------
    c : float
        Specific heat per atom of the system
    Ac_autocorr: float
        Error of the specific heat per atom with the autocorrelation method
    Ac_datablock: float
        Error of the specific heat per atom with the datablocking method
    """

    time, pos, vel = sim.load_data(file_name)
    num_tsteps = len(time) 
    particle_num = np.shape(vel)[1] 

    total_kin = (0.5*(vel**2).sum(2)).sum(1)
    mean_kin = total_kin.sum()/num_tsteps
    square_of_mean_kin = mean_kin**2
    mean_of_square_kin = (total_kin**2).sum()/num_tsteps
    mean_of_fourth_kin = (total_kin**4).sum()/num_tsteps
    r = mean_of_square_kin/square_of_mean_kin - 1 # relative fluctuations in kinetic energy

    c = 1.5/(1-3*particle_num*r/2)

    # Computation of the error with the autocorrelation function method
    tau_mean_of_square_kin = correlation_time(total_kin**2, time[-1], num_tsteps)
    tau_mean_kin = correlation_time(total_kin, time[-1], num_tsteps)

    err_mean_of_square_kin_autocorr = np.sqrt(2*tau_mean_of_square_kin/num_tsteps*(mean_of_fourth_kin - mean_of_square_kin**2))
    err_mean_of_kin_autocorr = np.sqrt(2*tau_mean_kin/num_tsteps*(mean_of_square_kin - square_of_mean_kin))

    Ac_autocorr = np.sqrt( (particle_num*(1/square_of_mean_kin)/(2/3*particle_num + 1 - mean_of_square_kin/square_of_mean_kin)**2)**2*err_mean_of_square_kin_autocorr**2 + \
         particle_num*mean_of_square_kin/square_of_mean_kin**(3/2)/(2/3*particle_num + 1 - mean_of_square_kin/square_of_mean_kin)**2*err_mean_of_kin_autocorr )

    # Computation of the error with the data-blocking method
    b_range = np.arange(1, min([500, int(len(time)/2)]), dtype=int)

    err_mean_of_square_kin_datablock = data_blocking(total_kin**2, b_range)[0]
    err_mean_of_kin_datablock = data_blocking(total_kin, b_range)[0]

    Ac_datablock = np.sqrt( (particle_num*(1/square_of_mean_kin)/(2/3*particle_num + 1 - mean_of_square_kin/square_of_mean_kin)**2)**2*err_mean_of_square_kin_datablock**2 + \
         particle_num*mean_of_square_kin/square_of_mean_kin**(3/2)/(2/3*particle_num + 1 - mean_of_square_kin/square_of_mean_kin)**2*err_mean_of_kin_datablock )
    
    return c, Ac_autocorr, Ac_datablock 


def pair_correlation_function_error(file_name, dr, box_length, r_max=None):
    """
    Returns the pair correlation function averaged over time and its respective error. 

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
    Ag_datablock : np.ndarray(int(r_max/dr)+1)
        Error of the pair correlation function as a function of r using data blocking method
    Ag_autocorr : np.ndarray(int(r_max/dr)+1)
        Error of the pair correlation function as a function of r using autocorrelation method
    """

    if r_max is None: r_max = box_length
    time, pos, _ = sim.load_data(file_name)
    particle_num = pos.shape[1]
    num_tsteps = len(time)
    r = np.arange(dr, r_max+dr, dr) # start != 0, because g = 1/0 is not defined
    n_time = np.zeros((num_tsteps, len(r)))

    for k, t in enumerate(time):
        print("\rGet data : {}/{}".format(k+1, len(time)), end="")
        rel_pos, rel_dist = sim.atomic_distances(pos[k], box_length)
        for i, r_i in enumerate(r):
            n_time[k,i] = len(np.where((rel_dist >= r_i) & (rel_dist < r_i + dr))[0])

    n = n_time.sum(0)
    g = 2*box_length**3 / (particle_num*(particle_num-1)) * (n/num_tsteps) / (4*np.pi*r**2 * dr)

    # Computation of the error with the data-blocking method
    n_error = np.zeros(len(r))
    b_range = np.arange(1, min([500, int(len(time)/2)]), dtype=int)
    for i, r_i in enumerate(r):
        print("\rCompute error (data blocking): {}/{}".format(i+1, len(r)), end="")
        n_error[i] = data_blocking(n_time[:,i], b_range)[0]

    Ag_datablock = 2*box_length**3 / (particle_num*(particle_num-1)) * (n_error) / (4*np.pi*r**2 * dr)

    # Computation of the error with the autocorrelation function method
    n_error = np.zeros(len(r))
    for i, r_i in enumerate(r):
        print("\rCompute error (autocorrelation): {}/{}".format(i+1, len(r)), end="")
        if n_time[:,i].sum() != 0: # if all of them are zeros, it is not possible to calculate autocorrelation function
            tau = correlation_time(n_time[:,i], time[-1], num_tsteps)
            mean = np.average(n_time[:,i])
            mean_squared = np.average(n_time[:,i]**2)
            n_error[i] = np.sqrt(2*tau/num_tsteps*(mean_squared - mean**2))

    Ag_autocorr = 2*box_length**3 / (particle_num*(particle_num-1)) * (n_error) / (4*np.pi*r**2 * dr)

    print("")

    return r, g, Ag_datablock, Ag_autocorr


def pressure_error(file_name, T, box_length):
    """
    Computes the pressure of the system. Gives the errors computed with the
    autocorrelation function and the data-blocking method

    Parameters
    ----------
    file_name : str
        Name of the CSV file in which the data is stored
    T : float
        Temperature
    box_length : float
		Box size

    Returns
    -------
    P : float
        Pressure
    AP_autocorr : float
        Error computed with the autocorrelation function
    AP_data_block : float
        Error computed with the data-blocking method
    """

    # Calculation of pressure
    time, pos, vel = sim.load_data(file_name)

    N = np.shape(pos)[1]
    M = len(time)

    second_term_instantenous = np.zeros(M)

    for k, t in enumerate(time):
        print("\r{}/{}".format(k+1, len(time)), end="")

        rel_pos, rel_dist = sim.atomic_distances(pos[k], box_length)

        rel_dist = rel_dist[:,:,np.newaxis] # add axis for LJ force calculation (so that it agrees with rel_pos dimensions)
        rel_dist[np.diag_indices(np.shape(rel_dist)[0])] = 1 # avoiding division by zero in the diagonal when calculating LJ force

        matrix = (1/(6*N*T))*24*(2/rel_dist**12-1/rel_dist**7)
        matrix[np.diag_indices(np.shape(matrix)[0])] = 0 # diagonal terms should be zero by definition

        second_term_instantenous[k] = matrix.sum()

    print("")

    second_term_average = np.average(second_term_instantenous)
    second_term_squared_average = np.average(second_term_instantenous**2)
    BP_rho =  1+second_term_average

    P = BP_rho*T*N/box_length**3

    # Computation of the error with the autocorrelation function method
    tau_second_term = correlation_time(second_term_instantenous, time[-1], M)
    sigma_second_term  = np.sqrt(2*tau_second_term*(second_term_squared_average-second_term_average**2)/M)
    AP_autocorr = N*T*sigma_second_term/box_length**3 # Through propagation of errors

    # Computation of the error with the data-blocking method
    b_range = np.arange(1, min([500, int(len(time)/2)]), dtype=int)
    error_second_term = data_blocking(second_term_instantenous, b_range)[0]
    AP_data_block = N*T*error_second_term/box_length**3 # Through propagation of errors

    return P, AP_autocorr, AP_data_block