import numpy as np

KB = 1.3806E-23 # Boltzmann constant in SI
SIGMA = 3.405E-10 # parameter of LJ potential for Argon atoms in SI
EPSILON = 119.8*KB # parameter of LJ potential for Argon atoms in SI
MASS = 6.6335209E-26 # Argon particle mass in SI


def simulate(init_pos, init_vel, num_tsteps, timestep, box_dim):
    """
    Molecular dynamics simulation using the Euler algorithm
    to integrate the equations of motion. 

    Parameters
    ----------
    init_pos : np.ndarra(N,2)
        Initial positions of the atoms in Cartesian space
    init_vel : np.ndarray(N,2)
        Initial velocities of the atoms in Cartesian space
    num_tsteps : int
        Total number of simulation steps
    timestep : float
        Duration of a single simulation step
    box_dim : float
        Dimensions of the simulation box

    Returns
    -------
    pos : np.ndarra(N,2)
        Positions of the atoms in Cartesian space after the simulation
    vel : np.ndarra(N,2)
        Velocities of the atoms in Cartesian space after the simulation
    """

    pos = init_pos
    vel = init_vel
    
    for k in range(num_tsteps):
        pos, vel = simulate_step(pos, vel, timestep, box_dim)

    return pos, vel


def simulate_step(pos, vel, timestep, box_dim):
    """
    Computes positions and velocities of particles at next time step
    using Euler method.

    Parameters
    ----------
    pos : np.ndarray(N,2)
        Positions of the atoms in Cartesian space
    vel : np.ndarray(N,2)
        Velocities of the atoms in Cartesian space
    timestep : float
        Duration of a single simulation step
    box_dim : float
        Dimensions of the simulation box

    Returns
    -------
    pos : np.ndarra(N,2)
        Positions of the atoms in Cartesian space after one time step
    vel : np.ndarra(N,2)
        Velocities of the atoms in Cartesian space after one time step
    """
    
    rel_pos, rel_dist = atomic_distances(pos, box_dim)
    total_force = lj_force(rel_pos, rel_dist)

    # Update velocities and positions with Euler method
    pos = pos + vel*timestep
    vel = vel + total_force*timestep/MASS

    # Update positions due to periodic BC
    pos = pos - np.floor(pos/box_dim)*box_dim

    return pos, vel
    

def atomic_distances(pos, box_dim):
    """
    Calculates relative positions and distances between particles.

    parameters
    ----------
    pos : np.ndarray(N,2)
        Positions of the particles in cartesian space
    box_dim : float
        Dimension of the simulation box

    returns
    -------
    rel_pos : np.ndarray(N,N,2)
        Relative positions of particles without taking into account periodic boundary conditions
        Let pos be the following matrix
            x1 y1
            x2 y2
            x3 y3
            .  .
            .  .
            .  .
        Then, rel_pos is
            0       x1-x2       x1-x3       ...         |       0       y1-y2       y1-y3       ...
            x2-x1   0           x2-x3       ...         |       y2-y1   0           y2-y3       ...
            x3-x1   x3-x2       0           ...         |       y3-y1   y3-y2       0           ...
            .
            .
            .
        where the separation between x and y relative positions is across the third axis of the matrix.

    rel_dist : np.ndarray(N,N)
        The distance between particles
        Each element is sqrt((xi-xj)**2+(yi-yj)**2), where i denotes row and j denotes column
    """

    N = np.shape(pos)[0] # Number of particles

    # Compute relative positions
    pos1 = np.repeat(pos[:, np.newaxis, :], N, axis=1)
    pos2 = np.repeat(pos[np.newaxis, :, :], N, axis=0)
    rel_pos = pos1 - pos2
    # check if using the minimum distance between pairs due to BC
    wrong_pairs = np.where(np.abs(rel_pos) > box_dim/2)
    rel_pos[wrong_pairs] = rel_pos[wrong_pairs] - np.sign(rel_pos[wrong_pairs])*box_dim

    # Compute relative distance
    rel_dist = np.linalg.norm(rel_pos, axis=2) # axis 2 contains the cartesian coordinates

    return  rel_pos, rel_dist


def lj_force(rel_pos, rel_dist):
    """
    Calculates the net forces on each atom.

    Parameters
    ----------
    rel_pos : np.ndarray(N,N,2)
        Relative particle positions as obtained from atomic_distances
    rel_dist : np.ndarray(N,N)
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    np.ndarray(N,2)
        Net force acting on each particle due to all other particles
    """

    rel_dist = rel_dist[:,:,np.newaxis] # add axis for LJ force calculation (so that it agrees with rel_pos dimensions)
    rel_dist[np.diag_indices(np.shape(rel_dist)[0])] = 1 # avoiding division by zero in the diagonal when calculating LJ force

    # Force calculation using the Lennard-Jones potential
    force = 24*EPSILON*rel_pos*(2*SIGMA**12/rel_dist**14 - SIGMA**6/rel_dist**8)
    total_force = force.sum(1)

    return total_force


def fcc_lattice(num_atoms, lat_const):
    """
    Initializes a system of atoms on an fcc lattice.

    Parameters
    ----------
    num_atoms : int
        The number of particles in the system
    lattice_const : float
        The lattice constant for an fcc lattice

    Returns
    -------
    pos_vec : np.ndarray
        Array of particle coordinates
    """

    return


def kinetic_energy(vel):
    """
    Computes the kinetic energy of an atomic system.

    Parameters
    ----------
    vel: np.ndarray(N,2)
        Velocities of particles

    Returns
    -------
    float
        Total kinetic energy of the system
    """

    return 0.5*MASS*(vel**2).sum()


def potential_energy(rel_dist):
    """
    Computes the potential energy of an atomic system.

    Parameters
    ----------
    rel_dist : np.ndarray(N,N)
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    float
        Total potential energy of the system
    """

    mask = np.array(1 - np.identity(rel_dist.shape[0]), dtype=bool)[:,:,np.newaxis] # mask for skipping diagonal terms
    pot_energy = 4*EPSILON*(SIGMA**12/rel_dist**12-SIGMA**6/rel_dist**6)
    pot_energy = np.sum(pot_energy, where=mask)/2 # The diagonal terms are skipped (no self-interaction)

    return pot_energy


def total_energy(rel_dist, vel): 
    """
    Computes the total energy of an atomic system

    Parameters
    ----------
    rel_dist : np.ndarray(N,N)
        Relative particle distances as obtained from atomic_distances
    vel: np.ndarray(N,2)
        Velocities of particles

    Returns
    -------
    float
        The total energy of the system
    """

    kin_energy = kinetic_energy(vel)
    pot_energy = potential_energy(rel_dist)

    return kin_energy + pot_energy


def init_velocity(num_atoms, temp):
    """
    Initializes the system with Gaussian distributed velocities.

    Parameters
    ----------
    num_atoms : int
        The number of particles in the system.
    temp : float
        The (unitless) temperature of the system.

    Returns
    -------
    vel_vec : np.ndarray
        Array of particle velocities
    """

    return
