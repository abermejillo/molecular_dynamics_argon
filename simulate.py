import numpy as np

KB = 1.3806E-23 # Boltzmann constant in SI
SIGMA = 3.405E-10 # $\sigma$ of L-J potential for Argon in SI
EPSILON = 119.8*KB # $\epsilon$ of L-J potential for Argon in SI
MASS = 6.6335209E-26 # Argon particle mass in SI

#--------------------------------------------------
# 1. Functions for the time evolution of the system
#--------------------------------------------------

def simulate_step_euler(pos, vel, timestep, box_dim):
    """
    Computes positions and velocities of particles at next time step
    using Euler method.
    N = number of particles
    d = dimensionality of the box

    Parameters
    ----------
    pos : np.ndarray(N,d)
        Positions of the particles in Cartesian space
    vel : np.ndarray(N,d)
        Velocities of the particles in Cartesian space
    timestep : float
        Duration of a single simulation step
    box_dim : float
        Dimensions of the simulation box

    Returns
    -------
    pos : np.ndarra(N,d)
        Positions of the particles in Cartesian space after one time step
    vel : np.ndarra(N,d)
        Velocities of the particles in Cartesian space after one time step
    """
    
    rel_pos, rel_dist = atomic_distances(pos, box_dim)
    total_force = lj_force(rel_pos, rel_dist)

    pos = pos + vel*timestep
    vel = vel + total_force*timestep

    # Update positions due to periodic BC
    pos = pos - np.floor(pos/box_dim)*box_dim

    return pos, vel


def simulate_step_verlet(pos, vel, timestep, box_dim):
    """
    Computes positions and velocities of particles at next time step
    using velocity-Verlet method.
    N = number of particles
    d = dimensionality of the box

    Parameters
    ----------
    pos : np.ndarray(N,d)
        Positions of the particles in Cartesian space
    vel : np.ndarray(N,d)
        Velocities of the particles in Cartesian space
    timestep : float
        Duration of a single simulation step
    box_dim : float
        Dimensions of the simulation box

    Returns
    -------
    pos : np.ndarra(N,d)
        Positions of the particles in Cartesian space after one time step
    vel : np.ndarra(N,d)
        Velocities of the particles in Cartesian space after one time step
    """
    
    rel_pos, rel_dist = atomic_distances(pos, box_dim)
    total_force = lj_force(rel_pos, rel_dist)

    pos = pos + vel*timestep + 0.5*timestep**2*total_force
    rel_pos, rel_dist = atomic_distances(pos, box_dim)
    total_force_next = lj_force(rel_pos, rel_dist)
    vel = vel + 0.5*(total_force + total_force_next)*timestep

    # Update positions due to periodic BC
    pos = pos - np.floor(pos/box_dim)*box_dim

    return pos, vel


def simulate(init_pos, init_vel, num_tsteps, timestep, box_dim, file_name, method="verlet"):
    """
    Molecular dynamics simulation using the Euler algorithm
    to integrate the equations of motion. 
    Saves data for each iteration using 'save_data' function. 
    N = number of particles
    d = dimensionality of the box

    Parameters
    ----------
    init_pos : np.ndarray(N,d)
        Initial positions of the particles in Cartesian space
    init_vel : np.ndarray(N,d)
        Initial velocities of the particles in Cartesian space
    num_tsteps : int
        Total number of simulation steps
    timestep : float
        Duration of a single simulation step
    box_dim : float
        Dimensions of the simulation box
    file_name : str
        Name of the CSV file to store the simulation data
    method : "verlet" or "euler"
        Selects method for the time evolution algorithm

    Returns
    -------
    None
    """

    pos, vel = init_pos, init_vel
    f = open(file_name, "w")
    f.write("N={} d={}\n".format(pos.shape[0], pos.shape[1])) # header
    save_data(f, 0, pos, vel) # saves initial position
    
    for k in np.arange(1, num_tsteps+1):
        print("\rtime step : {:7d}/{}".format(k, num_tsteps), end="")

        if method == "verlet":
            pos, vel = simulate_step_verlet(pos, vel, timestep, box_dim)
        if method == "euler":
            pos, vel = simulate_step_euler(pos, vel, timestep, box_dim)

        save_data(f, k*timestep, pos, vel)

    print("")
    f.close()

    return 


#--------------------------------
# 2. Thermalization of the system
#--------------------------------

def rescale_velocity(vel, T, T_current=None):
    """
    Rescales velocities to match the desired temperature
    N = number of particles
    d = dimensionality of the box

    Parameters
    ----------
    vel : np.ndarray(N,d)
        Velocities of particles
    T : float
        Desired temperature of the system
    T_current : float or None type
        Temperature of the system before the reescaling

    Returns
    -------
    resc_vel : np.ndarray(N,d)
        Rescaled velocities of particles
    T_current : float
        Temperature of the system before the reescaling
    """

    if T_current is None:
        T_current = temperature(vel)
    resc_factor = np.sqrt(T/T_current)
    resc_vel = resc_factor*vel

    return resc_vel, T_current


def get_equilibrium(init_pos, init_vel, max_num_tsteps, timestep, box_dim, T, file_name, method="verlet", resc_thr=[1E-2, 0.1]):
    """
    Given a temperature, initial positions and velocities, sets the system under equilibrium
    (up to some error given by the parameters in 'resc_thr'). 
    N = number of particles
    d = dimensionality of the box

    Parameters
    ----------
    init_pos : np.ndarray(N,d)
        Initial positions of the particles in Cartesian space
    init_vel : np.ndarray(N,d)
        Initial velocities of the particles in Cartesian space
    max_num_tsteps : int
        Maximum number of simulation steps before reaching equilibrium
    timestep : float
        Duration of a single simulation step
    box_dim : float
        Dimensions of the simulation box
    T : float
        Temperature of the system
    file_name : str
        Name of the CSV file to store the simulation data
    method : "verlet" or "euler"
        Selects method for the time evoluiton algorithm
    resc_thr : [float, float]
        First argument is the error for the temperature of the system (to do rescalings)
        Second argument is the minimum time between reescalings (to let the system equilibrate)

    Returns
    -------
    eq_reached : bool
        If True, equilibrium has been reached
    """

    pos, vel = init_pos, init_vel
    f = open(file_name, "w")
    f.write("N={} d={}\n".format(pos.shape[0], pos.shape[1])) # header
    save_data(f, 0, pos, vel) # saves initial position

    rescaling = [[1,0]] # [temperature of the system, time of the reescaling]
    eq_reached = False # if False, continues to check rescaling
    
    for k in np.arange(1, max_num_tsteps+1):
        print("\rtime step : {:7d}/{} | temperature={:0.5f}/{:0.5f} ({:0.5f})".format(k, max_num_tsteps, rescaling[-1][0], T, np.abs(1-np.sqrt(T/rescaling[-1][0]))), end="")

        if method == "verlet":
            pos, vel = simulate_step_verlet(pos, vel, timestep, box_dim)
        if method == "euler":
            pos, vel = simulate_step_euler(pos, vel, timestep, box_dim)

        # rescaling of vel for temperature equilibrium
        if (not eq_reached) and ((k*timestep - rescaling[-1][1]) > resc_thr[1]):
            T_current = temperature(vel) 
            if np.abs(T_current - T) > resc_thr[0]:
                vel, T_prev = rescale_velocity(vel, T, T_current=T_current)
                rescaling += [[T_prev, k*timestep]]
            else:
                eq_reached = True

        save_data(f, k*timestep, pos, vel)

        if eq_reached:
            break

    print("\rtime step : {:7d}/{} | temperature={:0.5f}/{:0.5f}".format(k, max_num_tsteps, T_current, T))

    f.close()

    return eq_reached


#------------------------------------
# 3. Functions for the time evolution
# ----------------------------------- 

def atomic_distances(pos, box_dim):
    """
    Calculates relative positions and distances between particles.
    N = number of particles
    d = dimensionality of the box

    parameters
    ----------
    pos : np.ndarray(N,d)
        Positions of the particles in cartesian space
    box_dim : float
        Dimension of the simulation box

    returns
    -------
    rel_pos : np.ndarray(N,N,d)
        Relative positions of particles without taking into account periodic boundary conditions
        Let pos be the following matrix [x1, x2, x3, ... ] where xi are d-dimensional vectors
        Then, rel_pos is
            0       x1-x2       x1-x3       ...        
            x2-x1   0           x2-x3       ...      
            x3-x1   x3-x2       0           ...      
            .
            .
            .
        where the separation between x and y relative positions is across the third axis of the matrix.

    rel_dist : np.ndarray(N,N)
        The distance between particles
        Each element is sqrt(sum_alpha (xi_alpha-xj_alpha)**2), 
        where i denotes row and j denotes column and alpha the cartesian coordinate
    """

    particle_num = np.shape(pos)[0] 

    # Compute relative positions
    pos1 = np.repeat(pos[:, np.newaxis, :], particle_num, axis=1)
    pos2 = np.repeat(pos[np.newaxis, :, :], particle_num, axis=0)
    rel_pos = pos1 - pos2
    # impose minimum distance between pairs due to periodic BC
    wrong_pairs = np.where(np.abs(rel_pos) > box_dim/2)
    rel_pos[wrong_pairs] = rel_pos[wrong_pairs] - np.sign(rel_pos[wrong_pairs])*box_dim

    # Compute relative distance
    rel_dist = np.sqrt(np.einsum('ijk, ijk->ij', rel_pos, rel_pos, optimize="optimal"))

    return  rel_pos, rel_dist


def lj_force(rel_pos, rel_dist):
    """
    Calculates the net forces on each particle interacting with a Lennard-Jones potential.
    N = number of particles
    d = dimensionality of the box

    Parameters
    ----------
    rel_pos : np.ndarray(N,N,d)
        Relative particle positions as obtained from atomic_distances
    rel_dist : np.ndarray(N,N)
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    total_force : np.ndarray(N,d)
        Net force acting on each particle due to all other particles
    """

    rel_dist = rel_dist[:,:,np.newaxis] # add axis for LJ force calculation so that it agrees with rel_pos dimensions
    rel_dist[np.diag_indices(np.shape(rel_dist)[0])] = 1 # avoids division by zero in the diagonal when calculating LJ force

    # Force calculation 
    force = 24*rel_pos*(2/rel_dist**14 - 1/rel_dist**8)
    total_force = force.sum(1)

    return total_force


#----------------------------------------
# 4. Position and velocity initialization
#----------------------------------------

def fcc_lattice(particle_num, lattice_const):
    """
    Initializes a system of particles on an FCC lattice.

    Parameters
    ----------
    particle_num : int
        The number of particles in the system
    lattice_const : float
        The lattice constant for an FCC lattice

    Returns
    -------
    pos : np.ndarray(particle_num,3)
        Coordinates of particles
    L : float
        Length of the box that encloses the lattice
    """
    b = int(((particle_num-1)/4)**(1/3)) + 1 # Minimum number of nodes that the lattice will have in each direction | (N-1) to solve when particle_num = 4*b^3

    # create position vectors for all cells
    xij = np.arange(0, b)
    combinations_vectors = np.transpose(np.meshgrid(xij, xij, xij), (2,1,3,0)).reshape(-1, 3)

    # fill cells with particles in FCC positions
    director_vectors = np.array([[0,0,0], [1,1,0], [1,0,1], [0,1,1]])*1/2
    pos = (director_vectors[0] + combinations_vectors)*lattice_const
    pos = np.append(pos, (director_vectors[1] + combinations_vectors)*lattice_const, axis=0)
    pos = np.append(pos, (director_vectors[2] + combinations_vectors)*lattice_const, axis=0)
    pos = np.append(pos, (director_vectors[3] + combinations_vectors)*lattice_const, axis=0)

    # shifts FCC away from the boundary
    pos = pos + lattice_const/4 

    # delete extra particles
    pos = pos[:particle_num]

    # box length
    L = b*lattice_const

    return pos, L


def random_gaussian_vector(num, sigma):
    """
    Returns a vector of length num with a Gaussian probability of zero mean and standard deviation sigma.

    Parameters
    ----------
    num : int
        Number of items in the vector
    sigma : float
        Standard deviation of the gaussian distribution

    Returns
    -------
    gaussian_dist : np.array(num,1)
        Array of gaussian distributed numbers with standard deviation sigma
    """

    u1 = np.random.rand(int(np.ceil(num/2)))
    u2 = np.pi*2*np.random.rand(int(np.ceil(num/2)))
    R = np.sqrt(-2*np.log(u1))
    x = R*np.cos(u2)
    y = R*np.sin(u2)
    gaussian_dist = np.concatenate((x,y))[0:num]
    gaussian_dist = sigma*(gaussian_dist - np.mean(gaussian_dist))/np.std(gaussian_dist)

    return gaussian_dist


def init_velocity(particle_num, T):
    """
    Initializes the system with Gaussian distributed velocities.
    N = number of particles
    d = dimensionality of the box

    Parameters
    ----------
    particle_num : int
        Number of particles in the system
    T : float
        Temperature of the system

    Returns
    -------
    vel : np.ndarray(N,d)
        Array of particle velocities
    """

    vel = np.array([random_gaussian_vector(particle_num, np.sqrt(T)),\
                    random_gaussian_vector(particle_num, np.sqrt(T)),\
                    random_gaussian_vector(particle_num, np.sqrt(T))])
    vel = vel.T
    
    return vel


#------------------------------------------------------------------
# 5. System's characteristics computation: Energies and temperature 
#------------------------------------------------------------------

def kinetic_energy(vel):
    """
    Computes the kinetic energy of the system.
    N = number of particles
    d = dimensionality of the box

    Parameters
    ----------
    vel: np.ndarray(N,d)
        Velocities of particles

    Returns
    -------
    float
        Total kinetic energy of the system
    """

    return 0.5*(vel**2).sum()


def potential_energy(pos, box_dim):
    """
    Computes the potential energy of the system.
    N = number of particles
    d = dimensionality of the box

    Parameters
    ----------
    pos : np.ndarray(N,d)
        Positions of the particles in cartesian space
    box_dim : float
        Dimension of the simulation box

    Returns
    -------
    pot_energy : float
        Total potential energy of the system
    """
    _, rel_dist = atomic_distances(pos, box_dim)
    rel_dist[np.diag_indices(np.shape(rel_dist)[0])] = 1 # avoids division by zero in the diagonal when calculating potential energy

    mask = np.array(1 - np.identity(rel_dist.shape[0]), dtype=bool) # mask for skipping diagonal terms (no self_interaction)
    pot_energy = 4*(1/rel_dist**12-1/rel_dist**6)
    pot_energy = np.sum(pot_energy, where=mask)/2

    return pot_energy


def total_energy(pos, vel, box_dim): 
    """
    Computes the total energy of the system.
    N = number of particles
    d = dimensionality of the box

    Parameters
    ----------
    pos : np.ndarray(N,d)
        Positions of the particles in cartesian space
    vel: np.ndarray(N,d)
        Velocities of particles
    box_dim : float
        Dimension of the simulation box

    Returns
    -------
    float
        The total energy of the system
    """

    kin_energy = kinetic_energy(vel)
    pot_energy = potential_energy(pos, box_dim)

    return kin_energy + pot_energy


def temperature(vel):
    """
    Returns the temperature of the system given its velocities. 
    N = number of particles
    d = dimensionality of the box

    Parameters
    ----------
    vel : np.ndarray(N,d)
        Velocities of particles

    Returns
    -------
    T : float
        (Temperature of the system
    """

    particle_num = vel.shape[0]
    T = (vel**2).sum()/(3*(particle_num-1))

    return T


#--------------------------------------------------
# 6. Saving and loading data to and from text files
#--------------------------------------------------

def save_data(file_class, time, pos, vel):
    """
    Writes to a CSV file the following information:
    time | pos_x1_alpha1, pos_x1_alpha2..., pos_x2_alpha1, pos_x2_alpha2, ... | vel_x1_alpha1, vel_x1_alpha2..., vel_x2_alpha1, vel_x2_alpha2, ... 
    (number of columns is 1 + d*N + d*N, where N = number of particles and d = dimensionality of the space)

    Parameters
    ----------
    file_class : TextIOWrapper 
        File in which to write the data, e.g. file_class = open("output.csv", "w")
    time : float
        Time of the current positions and velocities
    pos : np.ndarray(N,d)
        Positions of the particles
    vel: np.ndarray(N,d)
        Velocities of particles

    Returns
    -------
    file_class : TextIOWrapper 
        File in which to the data has been written
    """

    particle_num = pos.shape[0] 
    box_dim = pos.shape[1] 
    pos, vel = pos.reshape(-1), vel.reshape(-1) # reshape as: pos_x1, pos_y1, ... pos_x2, pos_y2, ... | vel_x1, vel_y1, ... vel_x2, vel_y2, ...

    data = "{:0.25e}".format(time)
    data += (",{:0.25e}"*box_dim*particle_num).format(*pos)
    data += (",{:0.25e}"*box_dim*particle_num).format(*vel)
    data += "\n"

    file_class.write(data)

    return file_class


def load_data(file_name):
    """
    Loads simulation data from CSV file that has the same structure
    as specified in 'save_data' function. 
    M = number of time steps
    N = number of particles
    d = dimensionality of the box

    Parameters
    ----------
    file_name : str
        Name of the CSV file in which the data is stored

    Returns
    -------
    time : np.ndarray(M)
        Time steps of the simulation
    pos : np.ndarray(M,N,d)
        Positions of the particles for all the time steps of the simulation
    vel: np.ndarray(M,N,d)
        Velocities of the particles for all the time steps of the simulation
    """

    # load header information
    f = open(file_name, "r")
    header = f.readline()
    f.close()
    N, d = header.split(" ")[:2]
    N, d = int(N.replace("N=", "")), int(d.replace("d=", ""))

    # load data
    data = np.loadtxt(file_name, delimiter=",", skiprows=1)
    M = data.shape[0]

    time = data[:,0]
    pos = data[:,1:d*N+1]
    pos = pos.reshape(M,N,d)
    vel = data[:,d*N+1:2*d*N+1]
    vel = vel.reshape(M,N,d)

    return time, pos, vel


def load_initial_data(file_name):
    """
    Loads initial data from CSV file that has the same structure
    as specified in 'save_data' function. 
    N = number of particles
    d = dimensionality of the box

    Parameters
    ----------
    file_name : str
        Name of the CSV file in which the data is stored

    Returns
    -------
    init_pos : np.ndarray(N,d)
        Initial positions of the particles
    init_vel: np.ndarray(N,d)
        Initial velocities of the particles
    """

    # load header information
    f = open(file_name, "r")
    header = f.readline()
    initial_data = f.readline()
    f.close()
    N, d = header.split(" ")[:2]
    N, d = int(N.replace("N=", "")), int(d.replace("d=", ""))

    # load data
    data = np.array([float(i) for i in initial_data.split(",")])
    init_pos = data[1:d*N+1]
    init_pos = init_pos.reshape(N,d)
    init_vel = data[d*N+1:2*d*N+1]
    init_vel = init_vel.reshape(N,d)

    return init_pos, init_vel


def load_final_data(file_name):
    """
    Loads final data from CSV file that has the same structure
    as specified in 'save_data' function. 
    N = number of particles
    d = dimensionality of the box

    Parameters
    ----------
    file_name : str
        Name of the CSV file in which the data is stored

    Returns
    -------
    final_pos : np.ndarray(N,d)
        Final positions of the particles
    final_vel: np.ndarray(N,d)
        Final velocities of the particles
    """

    # load header information
    f = open(file_name, "r")
    header = f.readline()
    data = f.read()
    final_data = data.split("\n")[-2] # last element is empty line
    f.close()
    N, d = header.split(" ")[:2]
    N, d = int(N.replace("N=", "")), int(d.replace("d=", ""))

    # load data
    data = np.array([float(i) for i in final_data.split(",")])
    final_pos = data[1:d*N+1]
    final_pos = final_pos.reshape(N,d)
    final_vel = data[d*N+1:2*d*N+1]
    final_vel = final_vel.reshape(N,d)

    return final_pos, final_vel