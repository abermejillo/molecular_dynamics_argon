import numpy as np
import matplotlib.pyplot as plt
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
        fig.tight_layout()
        fig.savefig("tmp-plot/pair_int_2D{:05d}.png".format(f))
        plt.cla() # clear axis

    plt.clf()
    print("\n", end="")

    print("BUILDING GIF... ")
    with imageio.get_writer(gif_name, mode='I', duration=1/30) as writer: # 30 fps
        for filename in ["tmp-plot/pair_int_2D{:05d}.png".format(f) for f in range(len(save_frame))]:
            image = imageio.imread(filename)
            writer.append_data(image)
    print("DONE")

    return


def plot_pos_2D(ax, pos, L, central_box=True):
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

    ax.set_xlim(-L/2, 3*L/2)
    ax.set_ylim(-L/2, 3*L/2)

    ax.set_xlabel("x coordinate [m]")
    ax.set_ylabel("y coordinate [m]")

    # set axis' ticks inside figure
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')

    return ax


def E_vs_t(data_file, box_dim):
    """
    Plots energy as a function of time

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

    energy = []
    for k, t in enumerate(time):
        energy += [sim.total_energy(pos[k], vel[k], box_dim)] 

    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.plot(time, energy, "-")

    ax.set_xlim(0, np.max(time))

    ax.set_xlabel("time [s]")
    ax.set_ylabel("energy [J]")

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

    energy = []
    for k, t in enumerate(time):
        energy += [sim.total_energy(pos[k], vel[k], box_dim)] 
    energy = np.array(energy)

    rel_AE = (energy - np.average(energy))/energy

    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.plot(time, rel_AE, "-")

    ax.set_xlim(0, np.max(time))

    ax.set_xlabel("time [s]")
    ax.set_ylabel("relative energy difference")

    # set axis' ticks inside figure
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')

    plt.show()
    plt.clf()

    return