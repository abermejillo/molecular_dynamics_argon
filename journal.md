# Weekly progress journal

## Instructions

In this journal you will document your progress of the project, making use of the weekly milestones.

Every week you should 

1. write down **on the day of the lecture** a short plan (bullet list is sufficient) of how you want to 
   reach the weekly milestones. Think about how to distribute work in the group, 
   what pieces of code functionality need to be implemented.
2. write about your progress **until Monday, 23:59** before the next lecture with respect to the milestones.
   Substantiate your progress with links to code, pictures or test results. Reflect on the
   relation to your original plan.

We will give feedback on your progress on Tuesday before the following lecture. Consult the 
[grading scheme](https://computationalphysics.quantumtinkerer.tudelft.nl/proj1-moldyn-grading/) 
for details how the journal enters your grade.

Note that the file format of the journal is *markdown*. This is a flexible and easy method of 
converting text to HTML. 
Documentation of the syntax of markdown can be found 
[here](https://docs.gitlab.com/ee/user/markdown.html#gfm-extends-standard-markdown). 
You will find how to include [links](https://docs.gitlab.com/ee/user/markdown.html#links) and 
[images](https://docs.gitlab.com/ee/user/markdown.html#images) particularly
useful.

## Week 1

### Bullet List

(1) Write `main.py` with some basic working code. This includes:
   - Write the numerical constants for the simulation.
   - Initializing uniformly random velocity and position.
   - Calculate relative position and relative distance matrices for the LJ force.
   - Iterate over time, simulating the velocities and positions using Euler method.
   - Every iteration, check if particles are outside of domain to fulfill periodic boundary conditions (if so, displace them accordingly).
   - Calculate energy to check if it is conserved along time.

(2) Once code works, divide it in function blocks (and use `skeleton.py` as template for the functions)

### Progress

All of the milestones set for this week have been completed. First, we wrote all the code in `main.py`. It initialized the particles and run the time evolution of the system [link to particular commit](https://gitlab.kwant-project.org/computational_physics/projects/Project-1_abermejillo_dbedialaunetar_mserraperalta/-/commit/e25d039188a0e8e7fceb45092361b36e0a65c9bd). This completes most of block (1) of the bullet list.

Afterwards, we optimized the code and changed the relative positions taking into account the periodic boundary conditions [link to particular commit](https://gitlab.kwant-project.org/computational_physics/projects/Project-1_abermejillo_dbedialaunetar_mserraperalta/-/commit/f7deb3540ece2e4ad3cc08ce28c869d4e06e876b).

Then, we added the computation of the total energy of the system [link to particular commit](https://gitlab.kwant-project.org/computational_physics/projects/Project-1_abermejillo_dbedialaunetar_mserraperalta/-/commit/037d679ccd07f29414b13573f86c39a880d6f394). With this we fully completed block (1).

Finally, we structured the code by adding all functionalities into functions in a file called `simulate.py` and added plotting functionalities in a file called `plotting.py`. Now, in `main.py` we only need to specify the simulation parameters and call the appropriate functions. This fulfils (2) in the bullet list.

Some of the results obtained with the actual version of the code are shown hereafter.

1. GIF showing the dynamics of the system:

![Sample Video](results/W1_dynamics_1.gif)

2. Energy conservation:
<img src="results/W1_energy_1.png" width="500" >

We can see how the particles evolve smoothly and interact with eachother (noticeable when they get close to eachother). If we look at the energy, we notice an apparently big step that would suggest that the energy is not conserved. However, we can observe that the step is 5 orders lower that that of the actual energy. Nonetheless, when using a higher number of particles ($`N \approx 20`$) there are problems with energy conservation if the box is not very large. We expect those to be solved upon normalization and using the velocity-Verlet algorithm. 

3. Other checks

- Take minimum distance between particles when using periodic BC, see [this gif](results/W1_closest_relative_distance.gif)
- Attractive force when particles are further than $`\sigma`$, see [this gif](results/W1_attractive_force.gif)

In conclusion, all milestones have been more than fulfilled. We have a code that works, it is quite well structured and the few tests that have been done give good results given the stage of the project.  

(due 14 February 2022, 23:59)


## Week 2

### Bullet List

1. Convert to dimensionless units (including deriving the expression of the kinetic energy and changing the scripts) @abermejillo
2. Implement the minimal image convention (done by @mserraperalta in Week 1)
3. Change space dimension from 2D to 3D @dbedialaunetar
4. Update existing plotting functions from 2D to 3D (if necessary) @abermejillo @mserraperalta @dbedialaunetar
5. Write or update the corresponding documentation and `README.md` (from 2D to 3D, units) @abermejillo @mserraperalta @dbedialaunetar
6. Start working with the velocity-Verlet algorithm @abermejillo @mserraperalta @dbedialaunetar
7. Start implementing the initialization of positions onto an fcc lattice and velocities with Maxwell-Boltzmann distribution @abermejillo @mserraperalta @dbedialaunetar

### Progress
Next, we will describe the progress made along the week based on the bullet list above and the commits submitted.

- The minimal image convention was done in Week 1 by @mserraperalta.
- We derived the expression of the kinetic energy in dimensionless units, and then adapted the molecular dynamics simulation to these units (this was done by @abermejillo in [Transform to natural units](b0466ee64fe75954858a147e49161d169a9d36cf)).
- The functions to carry out the simulation of the system were valid also for 3D (this was checked by @dbedialaunetar). The input vector of positions and velocities were changed to be $`N\times 3 `$ matrices, the functions `load_data` and `save_data` and the documentation were adapted to account for arbitrary dimensions (this was done by @mserraperalta in [Update load/save function and documentation for arbitrary dimensionality](a1b15f0a59d773a4bc1b2b85d6196d03a3020759)).
- We wrote functions to save a 3D figure out of every timestep in the simulation, and to create a GIF out of these figures. We then created the GIF for the simulation of two particles near the boundaries (this was done by @abermejillo in [Add 3D plotting](56f6cd4f707702fefe151824180e2f17294fe41d) and in [Change plotting and add results](c440d78073bd902530ca568360f9f85ea4e2e35a)), as asked for in the milestones, and which is discussed below.
- We updated the `README.md` file to account for dimensionless parameters and arbitrary dimensions (this was done by @mserraperalta in [Update README.md for dimensionless parameters and d dimensions](e40c1c927bf4b8c03645223ce6b79d47aa1ad5ef)).
- We implemented the velocity-Verlet algorithm (this was done by @mserraperalta [Add velocity-Verlet algorithm](a828f280fdca68a17a092c51ff95f84f4ccc115b)). We can now simulate our system using both the Euler method and velocity-Verlet algorithm through an argument passed to the function `simulate`.
- We wrote a first version of the function for initializing atoms in an fcc lattice (this was done by @dbedialaunetar in [Create fcc_lattice function](38da9384172583994f0fd535c701af9e8671e9ff)), which was then shortened by @mserraperalta in [Improve fcc lattice function](5adc7fdf7be261142f8e101f45550661bc8c40c0).
- We modified our total energy plotting function so that it also plots the kinetic and potential energies with respect to time (this was done by @mserraperalta in [Add kinetic and potential energy plots](5c485ae416519ef51dcc3517f5433c52ac4d5315)).
- We wrote a script to time every function involved in the simulation, thanks to which we determined a two times speed-up in the calculation of `rel_dist` by using a different numpy function (this was done @mserraperalta in [Add timing calculation for each step of the simulation](b86746f923d0a89d9194d42302f0e740b512e2fd) and in [Optimize rel_dist calculation by 2x](e6246894e87b9f9e6a1d1d0d8f8e1d7265850a14)).
- We wrote a function to plot the relative distance with respect to time, and then obtained one of these plots for two particles, as asked for in the milestones (this was done by @mserraperalta in [Add relative distance plot between 2 particles](27c67dd91be2717924efd072b5484d091ff7bf39) and in [Solve rel_dist error in plotting and add rel_dist results](04db495398d2ef79e3532dd45877986d65b56432)).

**Results and comments**

We now show and discuss the results for the simulation of two particles in 3D.

First, we have a GIF that shows how the two particles move inside a box with periodic boundary conditions.

![Sample Video](results/W2_2_atoms_close_to_the_boundary.gif)

The particles were initialized close to each other. At times $`t=0.9`$ and $`t=3.8`$ we see how one particle crosses a side of the box and appears on the opposite side, which shows that the periodic boundary conditions work.

We know plot the relative distance between the two particles over time.

<img src="results/W2_rel_dist_2_particles.png" width="450" >

At times $`t\approx 1.8`$ and $`t\approx 3.5`$ we see that the relative distance takes a minimum, which represent two scattering events due to the repulsion in the Lennard-Jones potential.

Finally, we plot the total, kinetic and potential energies with respect to time.

![alt text](results/W2_KE_PE_TE_2_atoms_close_to_the_boundary.png)

The most important conclusion from this graph is that the total energy is, indeed, conserved. We also see that the kinetic and potential energies have opposite extremal points, i.e, when there is a minimum in kinetic energy, there is a maximum in potential energy, and vice versa. Relating the relative distance graph with this one, we see that the wide minima in kinetic energy coincide with the instants where the two particles are furthest apart.

In order to better understand all the peaks that we see in energy, we have created the following GIF that relates each point in time with a point in the Lennard-Jones potential.

![Sample Video](results/W2_GIF_Lennard-Jones.gif)

All in all, it is clear that we have achieved all milestones and obtained all results from Week 2. Due to the fact that our simulations show repulsion and attraction between particles in addition to scattering events, it seems like the Lennard-Jones potential is being correctly simulated; and since the total energy is conserved, we can infer no mistakes are being made in the coding.


(due 21 February 2022, 23:59)


## Week 3

### Bulletlist

1. Improve the _FCC lattice_ function to omit a 'for' loop and change it so that the whole box is filled @mserraperalta
2. _Energy conservation_ comparison of the Euler and Velociy-Verlet algorithms (plots) @dbedialaunetar
3. Implement _Maxwell velocity distribution_ for the initialization @abermejillo
4. Implement _rescaling_ of the velocities @mserraperalta  
5. Investigate if an _adaptive step verlet algorithm_ can be implemented (and do so if it is) @abermejillo
6. Implement a couple of _observables_ @dbedialaunetar

### Progress 

First, we will briefly comment how the tasks from the bulletlist have been performed and by whom. 

1. @mserraperalta improved the [FCC lattice](ba0bff1808cd77a76a5e20e0302df5c13dd83f8a) and changed it's funcitonality so that the box is completely filled.
2. @dbedialaunetar plotted the [comparison of energy conservation](7a260d04f6b65d94c4cce513cd282ae7f5161976) with the two algorithms.
3. @abermejillo implemented the [initialization of the particles velocities](f44ddcde776c63d3378df038a3600f2d4896ae9c) according to a gaussian distribution and analysed it was correct by computing a [density probability function](dbe60d9bb168aa7b3fb4418ac4043e082581a361).
4. @mserraperalta implemented the [rescaling of the velocities](9e211f71a1801dac2880b04aa339a47baf336905) until the system is thermalised to a given temperature
5. @abermejillo investigated on implementing an adaptive step, but this has not given fruitful results so far. We might still look into it during week 4.
6. @dbedialaunetar implemented a function to compute the [specific heat](066122def1e0b341e5458a4177a4dafa287eccc5) and @mserraperalta implemented a function to compute the [pair correlation](7618162af91439ab11f3ba26bc8d122eb30a0471).

Along the way we did more commits that correspond to solving mistakes and typos as well as including functionalities for plotting and saving data in csv files.

**Results and comments**

All of the milestones have been completed succesfully. Our code already worked for more than 2 particles (Milestone 1). The velocity-Verlet was also already implemented, and we had a function to analyse the energy of the system. Comparing the two algorithms we come across the following figure, which is a simulation of two particles with same initial conditions, same runtime, and number of steps.

Euler algorithm             |  Verlet algorithm
:-------------------------:|:-------------------------:
![alt text](results/W3_energy_Euler.png) |  ![alt text](results/W3_energy_Verlet.png)

We can see how for the Euler algorithm energy is not conserved when a particle goes 'inside' other, that is, when the potential increases and the gradient blows off. The velocity-Verlet algorithm, however, manages to keep the total energy of the system constant. We also simulated the trajectories of two particles with both algorithms so that we can see in real space how the two algorithms lead to different time evolutions. Here the blue particles correspond to the Euler algorithm and the red ones to the velocity-Verlet algorithm.

![Sample Video](results/W3_movie_Verlet_Euler.gif)

We can see how they start following the same trajectory, but in a very short time they behave completely different.

With this all of the milestones of the week are fulfilled. We have been keeping our code as structured as possible along the whole project; the structure it has now will most probably be kept unchanged. 

Apart from this we show some results regarding milestones of following weeks. We have included the function that provides as an initial position distribution an FCC lattice, which can be seen in the following simulation where we can (more or less) see the FCC structure at a temperature of 119.8 K.

![sample video](results/W3_movie_FCClattice.gif)

And also included the initialisation for the Maxwell distributed velocities. The following graph shows a density probability function computed from a histogram of the velocities with 20 bins. On top of it there is a gaussian with the same mean (zero) and standard deviation.

![alt text](results/W3_Maxwell_distribution.png)

Then, we have some results for the two observables we have already implemented. In both cases we still have to analyse whether they match with the experimental results for Argon. For example, in the case of the specific heat we are obtaining a value of 3.11 1/atom, which differs from the 1.5 we should get for an ideal gas. For the case of the pair correlation function, we obtain the following results for the following simulation parameters: Nparticles=500, lattice_const=1.54478, box_dimension=7.7239, temperature=0.2, run_time=1, num_tsteps=1000, algorithm_method=verlet.

<img src="results/W3_pair_correlation_function.png" width="500" >

which has similar peaks to [Wikipedia pair correlation function of Argon](https://upload.wikimedia.org/wikipedia/commons/thumb/b/b8/Simulated_Radial_Distribution_Functions_for_Solid%2C_Liquid%2C_and_Gaseous_Argon.svg/425px-Simulated_Radial_Distribution_Functions_for_Solid%2C_Liquid%2C_and_Gaseous_Argon.svg.png?20151030020319), but the peaks don't exactly match (i.e. the first peak is not at $`r/\sigma=1`$). This will be our main task for the upcoming week.

Finally, some comments on the adaptive step algorithm we wanted to implement. The idea is to be able to, at each step, update the timestep so that it is small when the gradient is big and viceversa. In this way, we can optimize the code so that it has high precision but is quicker that just setting the timestep very small. However, we found out that such an algorithm brings many problems for a many-body system like this. Due to the number of particles involved there are always pairs of particles strongly interacting. In general, this leads to always choose a very small timestep. 

Either way, what we tried was implementing a correction on the timestep following the idea explained in the following paper: [Numerical Solutions of Nonlinear Ordinary Differential Equations
by Using Adaptive Runge-Kutta Method](https://www.researchgate.net/publication/335841661_Numerical_Solutions_of_Nonlinear_Ordinary_Differential_Equations_by_Using_Adaptive_Runge-Kutta_Method). They compare two algorithms to compute the error at each step. In their case it is between fifth and forth order runge kutta, while for us it is between Euler and velocity-Verlet algorithm. With that we can update the timestep according to the following formula

- if error $`\geq`$ tolerance then $`h_{new}=\delta h (\frac{tolerance}{error})^{0.2}`$ 

- if error $`<`$ tolerance then $`h_{new}=\delta h (\frac{tolerance}{error})^{0.25}`$ 

where $`error=0.5h^2F(x(t))`$, the difference between Euler and Verlet positions. Here we can see that this error is given by a vector, that is, we have an error for each particle. If we decided to take the maximum of this errors as the error, then we always get a super high value (several orders of magnitude greater than the box_size). The aproach taken was then to take the median of this array. Following it, we are able to change the timestep without going to ridiculously small values. 

So far, this approach does not give better results than our previous code. We might discuss it in this week's meeting and decide whether we want to put more effort into investigating different approaches for adaptive step solving of ODEs.


(due 28 February 2022, 23:59)


## Week 4

# Bullet list

1. Implement diffusion (@mserraperalta)
2. Implement pressure (@dbedialaunetar)
3. Check that our observable implementations are consistent with literature (@mserraperalta, @abermejillo, @dbedialaunetar)
4. Start computing errors, for example, the error for the specific heat (@abermejillo)
5. Start structuring the report: theoretical background, methods... (@mserraperalta, @abermejillo, @dbedialaunetar)

### Progress 

First, we will briefly comment how the tasks from the bulletlist have been performed and by whom. 

1. @mserraperalta implemented the [mean-squared displacement](ef2fc08806559e941f7110b48855a7325fa81c7e) the [diffusion function, and the plotting of the mean-squared displacement](a2d174fd8599b09a13c827c196554f40a4610a89).
2. @dbedialaunetar implemented the [pressure function](c55009bf0d1be6af57d1a39256257770cffd7386).
3. @mserraperalta added [plots of the pair correlation function](bb5077a34a52480379adb656d4d76a53a1f669c7), @abermejillo implemented the code to [plot the specific heat with respect to system parameters](5de38fc52d080c23892e911dac59ef0e98291914) and obtained numerous plots to get insight in the possible phase transitions, and @dbedialaunetar compared the result of the pressure with that given by the ideal gas law. The results will be discussed below.
4. @mserraperalta added the [autocorrelation function](51c035f3a428b63e21af1af4902c0a2aab185784), [data-blocking function](93ccd45f0aa38ed576b1c4e35991b3204cbe8397) and the [pair-correlation function with both errors](d5f8946902a9f84c6129335859c64d7101154634). @abermejillo added the [correlation-time function](ad68699e09161bb3f8586845b41a5193fbfe0ab3) and the specific heat function with both errors ([autocorrelation](ad68699e09161bb3f8586845b41a5193fbfe0ab3) and [data-blocking](dcd8992c91321f1d758b8e28d8a914a2fbf2686a)). @dbedialaunetar added the [pressure function with both errors](569677e70f4cacdc2e6507a1d84d73b8e770aee9).
5. @mserraperalta, @abermejillo and @dbedialaunetar discussed some of the results that should go into the report.

Along the way we did numerous commits that correspond to solving mistakes and typos, including functionalities for plotting and improve structure, ease of use and readability of the code.

**Results and comments**

First, we analysed the rescaling of the velocities. Here, we plot the kinetic energy of the system while it is undergoing the rescaling process.

![alt text](results/W4_rescaling.png)

The initial and final velocity distributions are

Before rescaling process ($`T=10.020`$)             |  After rescaling process ($`T=10.016`$)
:-------------------------:|:-------------------------:
![alt text](results/W4_velocities_before_rescaling.png) |  ![alt text](results/W4_velocities_after_rescaling.png)

Therefore, we obtain a thermal equilibrium configuration with the desired temperature. Moreover, the equilibrium temperature can be tuned to arbitrary precision (the more precision, the longer the rescaling process).

Next, we plot the specific heat with respect to temperature (for a lattice constant of $`2\cdot 1.5471 \sigma `$)

![alt text](results/W4_specific_heat_time.png)

We see how for large temperatures, the specific heat per atom tends towards $`1.5`$, which corresponds to the $`\frac{3}{2}k_B T`$ that is found for a free particle gas, due to the equipartition theorem. On the other side of the graph, we see that the specific heat "hovers" around $`3`$, which corresponds to the $`3 k_B T`$ of the equipartition theorem, due to the additional vibrational degrees of freedom. 

We are able to consistently get a specific heat of $`1.5`$ with a small error for the situations where argon is a gas. For instance, for $`T=1`$, lattice_constant$`=8.25`$ and $`N=500`$ particles, we obtain the following specific heat

$`c_v = 1.500\pm 0.005,`$

where the error was calculated using the autocorrelation function. However, we still have not found a region for solid argon, where the specific heat gives $`3`$ consistently. For this same simulation we find a dimensionless pressure

$`\frac{\beta P}{\rho} = 1.000\pm 0.005,`$

which further supports that Argon is behaving as a gas in this configuration.

We also obtained three plots of the pair correlation function.

![alt text](results/W4_pair_correlation_function_states.png)

If we compare these functions to the ones provided in [Wikipedia](https://upload.wikimedia.org/wikipedia/commons/thumb/b/b8/Simulated_Radial_Distribution_Functions_for_Solid%2C_Liquid%2C_and_Gaseous_Argon.svg/425px-Simulated_Radial_Distribution_Functions_for_Solid%2C_Liquid%2C_and_Gaseous_Argon.svg.png?20151030020319), we can see how our plots correspond to a solid, liquid and gas, respectively.

These results suggest that our simulation is going in the right direction. However, we still need to solve some relatively big discrepancies that we have found between the errors computed with the autocorrelation function and with the data-blocking method. Furthermore, we have to carry out a more extensive comparison of our results with literature, and understand if our solid and liquid argon are behaving as they should.


(due 7 March 2022, 23:59)


## Week 5
(due 14 March 2022, 23:59)
