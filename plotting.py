import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import imageio
import simulate as sim


def GIF_2D(gif_name, data_file, num_frames, box_dim):
    """
    Generates frames for the time evolution of particles in 2D
    and stores them in "tmp-plot" folder as "pair_int_2D{:05d}.png". 

    Parameters
    ----------
    gif_name : str
        Name of the GIF file to be generated
    data_file : str
        Name of the CSV file in which the data is stored
    num_frames : int
        Number of total frames generated (max is 99999)
    box_dim : float
        Dimensions of the simulation box

    Returns
    -------
    None
    """

    time, pos, vel = sim.load_data(data_file)
    num_tsteps = len(time) 
    save_frame = [int(i*(num_tsteps-1)/(num_frames-1)) for i in range(num_frames-1)] + [int(num_tsteps)-1] # timesteps in which to save frames

    if "tmp-plot" not in os.listdir():
        os.mkdir("tmp-plot")

    # Create figure and save initial position
    print("PLOTTING AND SAVING FRAMES... ({}/{})\r".format(1, num_frames), end="")
    fig = plt.figure(1)
    ax = fig.add_subplot(111)

    for f, t in enumerate(save_frame):
        print("PLOTTING AND SAVING FRAMES... ({}/{})\r".format(f+1, num_frames), end="")

        ax = plot_pos_2D(ax, pos[t], box_dim)
        ax.set_title("dimensionless t={:0.3f}".format(time[t]))
        fig.tight_layout()
        fig.savefig("tmp-plot/pair_int_2D{:05d}.png".format(f))
        plt.cla() # clear axis

    plt.clf()
    print("\n", end="")

    print("BUILDING GIF... ")
    with imageio.get_writer(gif_name, mode='I', duration=3/num_frames) as writer: # 30 fps
        for filename in ["tmp-plot/pair_int_2D{:05d}.png".format(f) for f in range(len(save_frame))]:
            image = imageio.imread(filename)
            writer.append_data(image)
    print("DONE")

    return

def GIF_3D(gif_name, data_file, num_frames, box_dim):
    """
    Generates frames for the time evolution of particles in 2D
    and stores them in "tmp-plot" folder as "pair_int_2D{:05d}.png". 

    Parameters
    ----------
    gif_name : str
        Name of the GIF file to be generated
    data_file : str
        Name of the CSV file in which the data is stored
    num_frames : int
        Number of total frames generated (max is 99999)
    box_dim : float
        Dimensions of the simulation box

    Returns
    -------
    None
    """

    time, pos, vel = sim.load_data(data_file)
    num_tsteps = len(time) 
    save_frame = [int(i*(num_tsteps-1)/(num_frames-1)) for i in range(num_frames-1)] + [int(num_tsteps)-1] # timesteps in which to save frames

    if "tmp-plot" not in os.listdir():
        os.mkdir("tmp-plot")

    # Create figure and save initial position
    print("PLOTTING AND SAVING FRAMES... ({}/{})\r".format(1, num_frames), end="")
    fig = plt.figure(1)
    ax = fig.add_subplot(111, projection='3d')

    for f, t in enumerate(save_frame):
        print("PLOTTING AND SAVING FRAMES... ({}/{})\r".format(f+1, num_frames), end="")

        ax = plot_pos_3D(ax, pos[t], box_dim)
        ax.set_title("dimensionless t={:0.3f}".format(time[t]))
        fig.tight_layout()
        fig.savefig("tmp-plot/pair_int_3D{:05d}.png".format(f))
        plt.cla() # clear axis

    plt.clf()
    print("\n", end="")

    print("BUILDING GIF... ")
    with imageio.get_writer(gif_name, mode='I', duration=3/num_frames) as writer: # 30 fps
        for filename in ["tmp-plot/pair_int_3D{:05d}.png".format(f) for f in range(len(save_frame))]:
            image = imageio.imread(filename)
            writer.append_data(image)
    print("DONE")

    return

def plot_pos_2D(ax, pos, L, central_box=True, relative_pos=False):
    """
    Plots positions of particles (and box) in 2D

    Parameters
    ----------
    ax : matplotlib axis
        Axis in which to plot the particles
    pos : np.ndarray(N,2)
        Positions of the atoms in Cartesian space
    L : float
        Dimensions of the simulation box
    central_box : bool
        If True, plots square box
    relative_pos : bool
        If True, plots line between closest pairs of all particles

    Returns
    -------
    ax : matplotlib axis
        Axis in which particles have been plotted
    """

    # plot central box and its eight neighbours
    for i in range(pos.shape[0]): # plot for all particles
        if central_box:
            ax.plot(pos[i,0]  , pos[i,1]  , ".", color="black") # central box
        else: 
            ax.plot(pos[i,0]  , pos[i,1]  , "r.") # central box
        ax.plot(pos[i,0]+L, pos[i,1]  , "r.")
        ax.plot(pos[i,0]  , pos[i,1]+L, "r.")
        ax.plot(pos[i,0]+L, pos[i,1]+L, "r.")
        ax.plot(pos[i,0]-L, pos[i,1]  , "r.")
        ax.plot(pos[i,0]  , pos[i,1]-L, "r.")
        ax.plot(pos[i,0]-L, pos[i,1]-L, "r.")
        ax.plot(pos[i,0]-L, pos[i,1]+L, "r.")
        ax.plot(pos[i,0]+L, pos[i,1]-L, "r.")

    if central_box: # plot square for central box
        ax.plot([0,0,L,L,0],[0,L,L,0,0], "g-") 

    if relative_pos:
        rel_pos, rel_dist = sim.atomic_distances(pos, L)
        for i in range(pos.shape[0]):
            for j in range(pos.shape[0]):
                if i == j: continue
                ax.plot([pos[i,0], pos[i,0]+rel_pos[j,i,0]], [pos[i,1], pos[i,1]+rel_pos[j,i,1]], "b--")

    ax.set_xlim(-L/2, 3*L/2)
    ax.set_ylim(-L/2, 3*L/2)

    ax.set_xlabel("dimensionless x coordinate")
    ax.set_ylabel("dimensionless y coordinate")

    # set axis' ticks inside figure
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')

    return ax

def plot_pos_3D(ax, pos, L, central_box=True, relative_pos=False, outer_boxes=False):
    """
    Plots positions of particles (and box) in 3D

    Parameters
    ----------
    ax : matplotlib axis
        Axis in which to plot the particles
    pos : np.ndarray(N,3)
        Positions of the atoms in Cartesian space
    L : float
        Dimensions of the simulation box
    central_box : bool
        If True, plots square box
    relative_pos : bool
        If True, plots line between closest pairs of all particles

    Returns
    -------
    ax : matplotlib axis
        Axis in which particles have been plotted
    """
    
    # plot central box and its eight neighbours
    for i in range(pos.shape[0]): # plot for all particles
        if central_box:
            ax.plot(pos[i,0]  , pos[i,1]  , pos[i,2]  , ".", color="black") # central box
        else: 
            ax.plot(pos[i,0]  , pos[i,1]  , pos[i,2]  , "r.") # central box
        if outer_boxes:
            ax.plot(pos[i,0]+L, pos[i,1]  , pos[i,2]  , "r.") # permutations + _ _
            ax.plot(pos[i,0]  , pos[i,1]+L, pos[i,2]  , "r.")
            ax.plot(pos[i,0]  , pos[i,1]  , pos[i,2]+L, "r.")

            ax.plot(pos[i,0]+L, pos[i,1]+L, pos[i,2]  , "r.") # permutations + + _
            ax.plot(pos[i,0]  , pos[i,1]+L, pos[i,2]+L, "r.")
            ax.plot(pos[i,0]+L, pos[i,1]  , pos[i,2]+L, "r.")

            ax.plot(pos[i,0]+L, pos[i,1]+L, pos[i,2]+L, "r.") # permutations + + + 

            ax.plot(pos[i,0]-L, pos[i,1]  , pos[i,2]  , "r.") # permutations - _ _
            ax.plot(pos[i,0]  , pos[i,1]-L, pos[i,2]  , "r.")
            ax.plot(pos[i,0]  , pos[i,1]  , pos[i,2]-L, "r.")

            ax.plot(pos[i,0]-L, pos[i,1]-L, pos[i,2]  , "r.") # permutations - - _
            ax.plot(pos[i,0]-L, pos[i,1]  , pos[i,2]-L, "r.")
            ax.plot(pos[i,0]  , pos[i,1]-L, pos[i,2]-L, "r.")

            ax.plot(pos[i,0]-L, pos[i,1]-L, pos[i,2]-L, "r.") # permutatinos - - - 
                
            ax.plot(pos[i,0]-L, pos[i,1]+L, pos[i,2]  , "r.") # permutations - + _
            ax.plot(pos[i,0]+L, pos[i,1]-L, pos[i,2]  , "r.")
            ax.plot(pos[i,0]-L, pos[i,1]  , pos[i,2]+L, "r.")
            ax.plot(pos[i,0]+L, pos[i,1]  , pos[i,2]-L, "r.")
            ax.plot(pos[i,0]  , pos[i,1]+L, pos[i,2]-L, "r.")
            ax.plot(pos[i,0]  , pos[i,1]-L, pos[i,2]+L, "r.")

            ax.plot(pos[i,0]+L, pos[i,1]-L, pos[i,2]-L, "r.") # permutations + - -
            ax.plot(pos[i,0]-L, pos[i,1]+L, pos[i,2]-L, "r.")
            ax.plot(pos[i,0]-L, pos[i,1]-L, pos[i,2]+L, "r.")

            ax.plot(pos[i,0]-L, pos[i,1]+L, pos[i,2]+L, "r.") # permutations + + -
            ax.plot(pos[i,0]+L, pos[i,1]-L, pos[i,2]+L, "r.")
            ax.plot(pos[i,0]+L, pos[i,1]+L, pos[i,2]-L, "r.")

            ax.set_xlim(-L, 2*L)
            ax.set_ylim(-L, 2*L)
            ax.set_zlim(-L, 2*L)

    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_zlim(0, L)
    if central_box: # plot square for central box
        ax.plot([0,L,L,0,0],[0,0,L,L,0],[0,0,0,0,0], "g-")
        ax.plot([0,0,L,L],[0,0,0,0],[0,L,L,0],"g-")
        ax.plot([0,0,0],[0,L,L],[L,L,0],"g-")
        ax.plot([0,L,L],[L,L,L],[L,L,0],"g-") 
        ax.plot([L,L],[0,L],[L,L],"g-")  

    if relative_pos:
        rel_pos, rel_dist = sim.atomic_distances(pos, L)
        for i in range(pos.shape[0]):
            for j in range(pos.shape[0]):
                    if i == j: continue
                    ax.plot([pos[i,0], pos[i,0]+rel_pos[j,i,0]], [pos[i,1], pos[i,1]+rel_pos[j,i,1]],[pos[i,2],pos[i,2]+rel_pos[j,i,2]], "b--")

    
    ax.set_xlabel("$x/\sigma$")
    ax.set_ylabel("$y/\sigma$")
    ax.set_zlabel("($z/\sigma$)")

    # set axis' ticks inside figure
    ax.tick_params(axis="x",direction="in")
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="z",direction="in")
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.zaxis.set_ticks_position('both')

    return ax


def E_vs_t(data_file, box_dim, kinetic_potential=False):
    """
    Plots energy as a function of time

    Parameters
    ----------
    data_file : str
        Name of the CSV file in which the data is stored
    box_dim : float
        Dimensions of the simulation box
    kinetic_potential : bool
        If True, plots also kinetic and potential energies

    Returns
    -------
    None
    """

    time, pos, vel = sim.load_data(data_file)

    E_total = []
    E_kinetic = []
    E_potential = []
    for k, t in enumerate(time):
        E_kinetic += [sim.kinetic_energy(vel[k])] 
        E_potential += [sim.potential_energy(pos[k], box_dim)]
        E_total += [E_kinetic[-1] + E_potential[-1]] 

    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.plot(time, E_total, "-", label="total E", color="black")
    if kinetic_potential:
        ax.plot(time, E_kinetic, "-", label="kinetic E", color="red")
        ax.plot(time, E_potential, "-", label="potential E", color="blue")
        ax.legend(loc="best")

    ax.set_xlim(0, np.max(time))

    ax.set_xlabel("dimensionless time")
    ax.set_ylabel("dimensionless energy")

    # set axis' ticks inside figure
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')

    plt.show()
    plt.clf()

    return


def E_conservation(data_file, box_dim):
    """
    Plots (E - average(E))/E as a function of time (where E = energy)

    Parameters
    ----------
    data_file : str
        Name of the CSV file in which the data is stored
    box_dim : float
        Dimensions of the simulation box

    Returns
    -------
    None
    """

    time, pos, vel = sim.load_data(data_file)

    E_total = []
    for k, t in enumerate(time):
        E_total += [sim.total_energy(pos[k], vel[k], box_dim)] 
    E_total = np.array(E_total)

    rel_AE = (E_total - np.average(E_total))/E_total

    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.plot(time, rel_AE, "-", color="black")

    ax.set_xlim(0, np.max(time))

    ax.set_xlabel("dimensionless time")
    ax.set_ylabel("relative energy difference")

    # set axis' ticks inside figure
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')

    plt.show()
    plt.clf()

    return


def reldist_vs_t(data_file, i, j, box_dim):
    """
    Plots relative distance between particles `i` and `j` as a function of time

    Parameters
    ----------
    data_file : str
        Name of the CSV file in which the data is stored
    i : int
        Number of the particle that constitutes the pair
    j : int
        Number of the particle that constitutes the pair
    box_dim : float
        Dimensions of the simulation box

    Returns
    -------
    None
    """

    time, pos, vel = sim.load_data(data_file)

    rel_dist = []
    for k, t in enumerate(time):
        rel_pos = pos[k, i, :] - pos[k, j, :]
        # check if using the minimum distance between pair due to BC
        wrong_pair = np.where(np.abs(rel_pos) > box_dim/2)
        rel_pos[wrong_pair] = rel_pos[wrong_pair] - np.sign(rel_pos[wrong_pair])*box_dim
        rel_dist += [np.linalg.norm(rel_pos)] 

    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.plot(time, rel_dist, "-", color="black")

    ax.set_xlabel("dimensionless time")
    ax.set_ylabel("relative_distance(i={},j={})/$\sigma$".format(i,j))

    # set axis' ticks inside figure
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')

    fig.tight_layout()

    plt.show()
    plt.clf()

    return

def GIF_potential_energy(gif_name, data_file, num_frames, i , j, box_dim):
    """
    Generates frames for the time evolution of the potential energy 
    of a pair of particles in a E_vs_t graph and a Lennard-Jones potential graph.
    The frames are saved in "LJ_gif/", while the GIF is saved in the main directory.

    Parameters
    ----------
    gif_name : str
        Name of the GIF which is created.
    data_file : str
        Name of the CSV file in which the data is stored
    num_frames : int
        Number of total frames generated (max is 99999)
    i : int
        Number of the particle that constitutes the pair
    j : int
        Number of the particle that constitutes the pair
    box_dim : float
        Dimensions of the simulation box

    Returns
    -------
    None
    """

    if "LJ-gif" not in os.listdir():
        os.mkdir("LJ-gif")

    time, pos, vel = sim.load_data(data_file)

    num_tsteps = len(time)
    save_frame = [int(i*(num_tsteps-1)/(num_frames-1)) for i in range(num_frames-1)] + [int(num_tsteps)-1] # timesteps in which to save frames
    print("PLOTTING AND SAVING FRAMES... ({}/{})\r".format(1, num_frames), end="")

    # Extract energies and relative distances at each timestep
    E_total = []
    E_kinetic = []
    E_potential = []
    rel_dist = []

    for k, t in enumerate(time):
        rel_pos = pos[k, i, :] - pos[k, j, :]
        wrong_pair = np.where(np.abs(rel_pos) > box_dim/2)
        rel_pos[wrong_pair] = rel_pos[wrong_pair] - np.sign(rel_pos[wrong_pair])*box_dim
        rel_dist += [np.linalg.norm(rel_pos)]

        E_kinetic += [sim.kinetic_energy(vel[k,i])+sim.kinetic_energy(vel[k,j])]
        E_potential += [4*(1/rel_dist[-1]**12-1/rel_dist[-1]**6)]
        E_total += [E_kinetic[-1] + E_potential[-1]]

    # Generate data to plat the Lennard-Jones potential from rmin to rmax
    rmin = np.min(rel_dist)
    rmax = np.max(rel_dist)
    N = 100
    rel_dist_LJ = np.linspace(rmin, rmax, N)
    LJ_potential = 4*(1/rel_dist_LJ**12-1/rel_dist_LJ**6)

    # Create figures at each timestep specified by save_frame
    time2 = [time[i] for i in save_frame]
    E_kinetic2 = [E_kinetic[i] for i in save_frame]
    E_potential2 = [E_potential[i] for i in save_frame]
    E_total2 = [E_total[i] for i in save_frame]

    for k, t in enumerate(save_frame):
        print("PLOTTING AND SAVING FRAMES... ({}/{})\r".format(k+1, num_frames), end="")

        fig, axis = plt.subplots(2)

        axis[0].scatter(time[t], E_kinetic[t], color = "red")
        axis[0].scatter(time[t], E_potential[t], color = "blue")
        axis[0].scatter(time[t], E_total[t], color = "black")
        axis[0].plot(time2,E_kinetic2, color = "red")
        axis[0].plot(time2,E_potential2, color = "blue")
        axis[0].plot(time2,E_total2, color = "black")
        axis[0].set_xlabel("dimensionless time")
        axis[0].set_ylabel("dimensionless energy")
        axis[1].scatter(rel_dist[t],E_potential[t], color = "blue")
        axis[1].plot(rel_dist_LJ,LJ_potential,"-",color="blue")
        axis[1].set_xlabel("dimensionless relative distance")
        axis[1].set_ylabel("dimensionless potential energy")

        axis[0].set_title("dimensionless t={:0.3f}".format(time[t]))
        fig.tight_layout()
        fig.savefig("LJ-gif/pair_pot_3D{:05d}.png".format(k))
        plt.cla() # clear axis
    
    plt.clf()
    print("\n", end="")

    print("BUILDING GIF... ")
    with imageio.get_writer(gif_name, mode='I', duration=3/num_frames) as writer: # 30 fps
        for filename in ["LJ-gif/pair_pot_3D{:05d}.png".format(f) for f in range(len(save_frame))]:
            image = imageio.imread(filename)
            writer.append_data(image)
    print("DONE")
    
    return

def merge_GIF_3D(gif_name, data_file1, data_file2, num_frames, box_dim):
    """
    Generates frames for the time evolution of particles in 2D, merging
    two different simulations and stores them in "tmp-plot" folder as 
    "pair_int_2D{:05d}.png". 

    Parameters
    ----------
    gif_name : str
        Name of the GIF file to be generated
    data_file1 : str
        Name of the CSV file in which the data is stored
    data_file1 : str
        Name of the CSV file in which the data is stored
    num_frames : int
        Number of total frames generated (max is 99999)
    box_dim : float
        Dimensions of the simulation box

    Returns
    -------
    None
    """

    time, pos1, vel1 = sim.load_data(data_file1)
    time, pos2, vel2 = sim.load_data(data_file2)
    num_tsteps = len(time) 
    save_frame = [int(i*(num_tsteps-1)/(num_frames-1)) for i in range(num_frames-1)] + [int(num_tsteps)-1] # timesteps in which to save frames

    if "tmp-plot" not in os.listdir():
        os.mkdir("tmp-plot")

    # Create figure and save initial position
    print("PLOTTING AND SAVING FRAMES... ({}/{})\r".format(1, num_frames), end="")
    fig = plt.figure(1)
    ax = fig.add_subplot(111, projection='3d')

    for f, t in enumerate(save_frame):
        print("PLOTTING AND SAVING FRAMES... ({}/{})\r".format(f+1, num_frames), end="")

        ax = plot_pos_3D(ax, pos1[t], box_dim, color="blue")
        ax = plot_pos_3D(ax, pos2[t], box_dim, color="red")
        ax.set_title("dimensionless t={:0.3f}".format(time[t]))
        fig.tight_layout()
        fig.savefig("tmp-plot/pair_int_3D{:05d}.png".format(f))
        plt.cla() # clear axis

    plt.clf()
    print("\n", end="")

    print("BUILDING GIF... ")
    with imageio.get_writer(gif_name, mode='I', duration=3/num_frames) as writer: # 30 fps
        for filename in ["tmp-plot/pair_int_3D{:05d}.png".format(f) for f in range(len(save_frame))]:
            image = imageio.imread(filename)
            writer.append_data(image)
    print("DONE")

    return

def plot_maxwell_distribution(init_vel, temp):
    """
    Generates a plot that shows a probability density of a a particle having a given velocity. 
    A gaussian distribution with standard deviation \sqrt(temperature) is shown on top of it.

    Parameters
    ----------
    init_vel : np.array(N,d)
        Initial distribution of velocities
    temp : float
        temperature of the system in units of KB/epsilon
    
    Returns
    -------
    None
    """

    plt.hist(init_vel[:,0],bins=20,density=True,label='Hist of velocities') 
    sigma = np.sqrt(temp)
    x = np.linspace(-3*sigma,3*sigma, 100)
    plt.plot(x, stats.norm.pdf(x, 0, sigma),label='Gauss$(\mu=0, \sigma=\sqrt{T})$')
    plt.xlabel("Velocity (dimensionless)")
    plt.ylabel("Density of probability")
    plt.legend(loc='upper left')
    plt.show()
    plt.clf()
    return 