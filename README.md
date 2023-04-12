# Modelling a Carbon Nanotube Spin Qubit Experiencing Electric Dipole Spin Resonance 

This repository contains the code written in my MSci Physics project. My project simulates the reponse of a carbon nanotube spin qubit experiencing Electric Dipole Spin Resonance (EDSR). The first generation design of the carbon nanotube (CNT) quantum dot device used by the Quantum Devices was simulated. Hence, a realistic electric field and magnetic field were implemented into the simulation. 


### Magnetic Field Profile.txt
Mumax3 code to obtain the effective magnetic field in the $x$, $y$ and $z$ direction due to 3 micromagnets. 

### Magnetic Field Plot.py
A Python code for extracting and plotting the magnetic field profiles from the Mumax3 simulation. 

## `Tkwant CNT EDSR Simulator.py'
A Python code to define the CNT system in Kwant and evolve the ground state using Tkwant. The full EDSR simulation can be run using the `__init__` function. The system is defined by a set number of lattice points and the full simulation time specified in seconds. 

###### The `__init__` function:

**Arguements**:
* `hamiltonian`: the full EDSR Hamiltonian to describe the system
* `perturbation_type`: the type of time-dependent potential - specify as 'sin' or 'cos'
* `electric_field_type`: the type of electric field perturbation - specfiy as 'see-saw' or 'gaussian'
* `number_of_lattice`: the number of lattice points to define the CNT scattering region
* `magnetic_field_file`: file name of real magnetic field profile obtained from Mumax3
* `potential_type`: type of confinement potential - specify as 0: infinite square well or 1: parabolic potential

**Output**:
* Probability density function graph of the qubit states at $t = 0$.
* Eigenenergy level diagram of the qubit and excited states at $t = 0$.
* The PDF, electric field potential and expectation values of the Pauli spin matrices for each time step.
* the effective magnetic field in the x-direction calculated using the PDF values.
* Effective magnetic field in the x-direction against time graph.
* Rabi oscillation graph. 
* Animation of the PDF and electric field potential evovling over time.
* Animation of the PDF, the electric field potential, the effective magnetic field in the x-direction and the Rabi oscillations over time.
All the data and graphs are saved. 

**Advice**:
* For initial runs of the code I recommend defining the CNT system using 20 or 40 lattice points. In addition, evolve the state for a very small time range to check that the evolution occurs as expected, for example evolve up to $t = 1 \times 10^{-11}$ for 2 time steps. 
* A break can be implemented in the evolution loop to evaluate the system after a certain time frame.


## `Qiskit CNT EDSR Simulator.py'
A Python code to evolve an initial eigenvector using Qiskit Dynamics for a effective qubit Hamiltonian defined by parameters from Kwant. 
The parameters .. are defined by the Kwant system in the code .. and imported into the Qiskit simulation. 

The code outputs the Rabi Oscillation plot showing the evolution of the qubit state defined by $|0\rangle$. The code also outputs an animation of the Rabi oscillations over time and the corresponding oscillation about the Bloch Sphere.
All graphs and animations are saved.
