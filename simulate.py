"""
This is a suggestion for structuring your simulation code properly.
However, it is not set in stone. You may modify it if you feel like
you have a good reason to do so.
"""

import numpy as np


def simulate(init_pos, init_vel, num_tsteps, timestep, box_dim):
    """
    Molecular dynamics simulation using the Euler or Verlet's algorithms
    to integrate the equations of motion. Calculates energies and other
    observables at each timestep.

    Parameters
    ----------
    init_pos : np.ndarray
        The initial positions of the atoms in Cartesian space
    init_vel : np.ndarray
        The initial velocities of the atoms in Cartesian space
    num_tsteps : int
        The total number of simulation steps
    timestep : float
        Duration of a single simulation step
    box_dim : float
        Dimensions of the simulation box

    Returns
    -------
    Any quantities or observables that you wish to study.
    """

    N = len(init_pos(:,0)) # Number of particles
    pos = init_pos[:] # Create a new object with the positions, that will be updated in the iterations
    vel = init_vel[:] # Create a new object with the velocities, that will be updated in the iterations
    run_time = num_tsteps * timestep
    L = box_dim #Shorten name of variable
    
    for k, t in enumerate(np.arange(0, run_time + timestep, timestep)):

        #Compute relative positions and distances
        rel_pos_, rel_pos, rel_dist = atomic_distances(pos, L)

        # Total force exerted on each particle
        total_force = lj_force(rel_pos, rel_dist)

        # Update velocities and positions
        pos = pos + vel*timestep
        vel = vel + total_force*timestep/MASS

        # Check periodic boundary conditions
        pos = pos - np.floor(pos/L)*L

        # plotting for checking the interactions
        if True:
            for i in range(N):
                # plot normal box and its eight neighbours
                plt.plot(pos[i,0], pos[i,1], "r.")
                plt.plot(pos[i,0]+L, pos[i,1], "r.")
                plt.plot(pos[i,0], pos[i,1]+L, "r.")
                plt.plot(pos[i,0]+L, pos[i,1]+L, "r.")
                plt.plot(pos[i,0]-L, pos[i,1], "r.")
                plt.plot(pos[i,0], pos[i,1]-L, "r.")
                plt.plot(pos[i,0]-L, pos[i,1]-L, "r.")
                plt.plot(pos[i,0]-L, pos[i,1]+L, "r.")
                plt.plot(pos[i,0]+L, pos[i,1]-L, "r.")

                for j in range(i):
                    plt.plot([pos[j,0], pos[j,0] + rel_pos_[i,j,0]], [pos[j,1], pos[j,1] + rel_pos_[i,j,1]], "y-") # no BC correction
                    plt.plot([pos[j,0], pos[j,0] + rel_pos[i,j,0]], [pos[j,1], pos[j,1] + rel_pos[i,j,1]], "b--") # correction for BC

            plt.plot([0,0,L,L,0],[0,L,L,0,0], "g-") # plot square for normal box
            plt.xlim(-L/2, 3*L/2)
            plt.ylim(-L/2, 3*L/2)
            plt.savefig("{}.png".format(k))
            plt.cla()

    return


def atomic_distances(pos, box_dim):
    """
    Calculates relative positions and distances between particles.

    parameters
    ----------
    pos : np.ndarray
        The positions of the particles in cartesian space
    box_dim : float
        The dimension of the simulation box

    returns
    -------
    rel_pos_ : np.ndarray
        Relative positions of particles without taking into account periodic boundary conditions (NxNx2 array)
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
    
    rel_pos : np.ndarray
        Relative positions of particles taking into account periodic boundary conditions (NxNx2 array)

    rel_dist : np.ndarray
        The distance between particles (NxN matrix)
        Each element is sqrt((xi-xj)**2+(yi-yj)**2), where i denotes row and j denotes column
    """

    N = len(pos(:,0)) # Number of particles
    L = box_dim #Shorten name of variable

    # Compute relative positions
    pos1 = np.repeat(pos[:, np.newaxis, :], N, axis=1) # NxNx2 matrix
    pos2 = np.repeat(pos[np.newaxis, :, :], N, axis=0) # NxNx2 matrix

    rel_pos = pos1 - pos2

    rel_pos_ = rel_pos.copy() # for plotting

    # minimum image convention (for the periodic boundary conditions)
    wrong_pairs = np.where(np.abs(rel_pos) > L/2)
    rel_pos[wrong_pairs] = rel_pos[wrong_pairs] - np.sign(rel_pos[wrong_pairs])*L

    # Compute relative distance
    rel_dist = np.linalg.norm(rel_pos, axis=2) # axis 2 contains the cartesian coordinates

    return  rel_pos_, rel_pos, rel_dist


def lj_force(rel_pos, rel_dist):
    """
    Calculates the net forces on each atom.

    Parameters
    ----------
    rel_pos : np.ndarray
        Relative particle positions as obtained from atomic_distances
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    np.ndarray
        The net force acting on particle i due to all other particles
    """

    rel_dist = rel_dist[:,:,np.newaxis] # add axis for LJ force calculation (so that it agrees with rel_pos dimensions)
    rel_dist[np.diag_indices(N)] = 1 # avoiding division by zero in the diagonal when calculating LJ force

    # Force calculation using the Lennard-Jones potential
    force = 24*EPSILON*rel_pos*(2*SIGMA**12/rel_dist**14 - SIGMA**6/rel_dist**8)

    # Total force exerted on each particle
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
    vel: np.ndarray
        Velocity of particle

    Returns
    -------
    float
        The total kinetic energy of the system.
    """

    return


def potential_energy(rel_dist):
    """
    Computes the potential energy of an atomic system.

    Parameters
    ----------
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    float
        The total potential energy of the system.
    """

    return


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
