import kwant.continuum
import numpy as np
import matplotlib.pyplot as plt
import kwant
from numpy import vectorize
from matplotlib import animation
from tkwant import onebody
import tkwant
import os.path
import json
import time as timelib
from matplotlib import patches

from csv import writer

class System:
    def __init__(self, hamiltonian, pertubation_type="sin", magnetic_field_file="none", electric_field_type="see-saw", 
                 number_of_lattices=50, potential_type=0, file_name='none', folder_path='none'):

        # read arguments
        self.potential_type = potential_type
        self.hamiltonian = hamiltonian
        self.number_of_lattices = number_of_lattices
        self.magnetic_field_file = magnetic_field_file
        self.electric_field_type = electric_field_type
        self.pertubation_type = pertubation_type
        self.file_name = file_name
        self.folder_path = folder_path
        self.evolve_state = False

        # constants in S.I.
        self.g = 2
        self.hbar_SI = 1.054571817e-34 # in Js
        self.e_SI = 1.602176634e-19 # in J
        self.a_0_SI = 5.2917721090380e-11 # in m
        self.total_length_SI = 0.66e-6 # in m
        self.m_e_SI = 9.11e-31 # in kg
        self.mu_B_SI = 9.2740100783e-24 # in J/T
        self.lattice_size_SI = self.total_length_SI / self.number_of_lattices
        self.z_SI = np.arange(-self.total_length_SI / 2, self.total_length_SI / 2, self.lattice_size_SI,
                              dtype=np.double)  # in m
        
        # hamiltonian parameters in S.I.
        self.B_0_SI = 0.5 # in T - external field
        self.E_sl_eV = 1e-6 # in eV -  slanting field
        self.eV_0_eV = 1e-3 # in eV - electric field strength
        self.omega_0_eV = 1e-3  # in eV - parabolic potential 
        
        # constants in a.u.
        self.hbar_au = 1
        self.e_au = 1
        self.total_length_au = self.m_to_au(self.total_length_SI)  
        self.m_au = 0.01 # effective mass
        self.mu_B_au = .5
        self.lattice_size_au = self.total_length_au / self.number_of_lattices  # Distance in a.u. between lattice points
        self.z_au = np.arange(-self.total_length_au / 2, self.total_length_au / 2, self.lattice_size_au,
                              dtype=np.double)    
        
        # hamiltonian parameters in a.u.
        self.B_0_au = self.tesla_to_au(self.B_0_SI) 
        self.E_sl_au = self.ev_to_hartree(self.E_sl_eV)
        self.eV_0_au = self.ev_to_hartree(self.eV_0_eV)
        self.omega_0_au = self.ev_to_hartree(self.omega_0_eV)
        
        # slanted magnetic field value computed from E_sl in S.I. and a.u.
        self.b_sl_au = self.E_sl_au / (self.g * self.mu_B_au * self.total_length_au) 
        self.b_sl_SI = (self.E_sl_eV*self.e_SI)/ (self.g * self.mu_B_SI * self.total_length_SI) # in T/m
        
        # driving frequency in S.I. and a.u.        
        self.pulse_frequency_SI = (self.g * self.B_0_SI * self.mu_B_SI) / (2 * np.pi * self.hbar_SI)
        self.pulse_frequency_au = (self.g * self.B_0_au * self.mu_B_au) / (2 * np.pi * self.hbar_au)
    
    def cosine_v_ac(self, time, z, eV_0_au, pulse_frequency_au, total_length_au):
        """
        Function defining the electric field perturbation to the EDSR Hamiltonian for a cosine alternating voltage
        Inputs:
            time (a.u): time step in a.u
            z (a.u): position along quantum dot
            eV_0_au (a.u): stength of the oscillating electric field
            pulse_frequency_au (a.u): frequency of oscillation
            total_length_au (a.u): total length of quantum wire
        Output: the time-dependent pertubation due to electric field
        """
        # realistic Gaussian-shaped electric field
        if self.electric_field_type == "gaussian":
            # position gaussian pulse above the gate electrode
            self.z_shift = -self.m_to_au(170e-9)
            self.sigma = self.m_to_au(300e-10)
            return (((eV_0_au * np.cos(2 * np.pi *pulse_frequency_au * time)) * np.exp((-(z-self.z_shift)**2)/(2*self.sigma**2))))
        
        # simplified 'see-sawing' electric field
        elif self.electric_field_type == "see-saw":
            return ((eV_0_au * np.cos(2 * np.pi * pulse_frequency_au * time) * z) / total_length_au)
            
        else:
            print("Incorrect type of electric field perturbation specified - check spelling")
            
    def sine_v_ac(self, time, z, eV_0_au, pulse_frequency_au, total_length_au):
        """
        Function defining the electric field perturbation to the EDSR Hamiltonian for a sine alternating voltage
        Inputs:
            time (a.u): time step in a.u
            z (a.u): position along quantum dot
            eV_0_au (a.u): stength of the oscillating electric field
            pulse_frequency_au (a.u): frequency of oscillation
            total_length_au (a.u): total length of quantum wire
        Output: the time-dependent pertubation due to electric field
        """
        
        # realistic Gaussian-shaped electric field
        if self.electric_field_type == "gaussian":
            # position gaussian pulse above the gate electrode
            self.z_shift = -self.m_to_au(170e-9)
            self.sigma = self.m_to_au(300e-10)
            return (((eV_0_au * np.sin(2 * np.pi *pulse_frequency_au * time)) * np.exp((-(z-self.z_shift)**2)/(2*self.sigma**2))))
        
        # simplified 'see-sawing' electric field
        elif self.electric_field_type == "see-saw":
            return ((eV_0_au * np.sin(2 * np.pi * pulse_frequency_au * time) * z) / total_length_au)
            
        else:
            print("Incorrect type of electric field perturbation specified - check spelling")
   
    ## conversion functions
    def tesla_to_au(self, tesla):
        """
        Function to convert the magnetic flux density from SI to a.u. units
        tesla: magnetic flux density in teslas.
        return: magnetic flux density in a.u.
        """
        return tesla / 2.35e5

    def au_to_tesla(self, au):
        """
        Function to convert the magnetic flux density from a.u. to SI units
        au: magnetic flux density in a.u.
        return: magnetic flux density in teslas.
        """
        return au * 2.35e5

    def second_to_au(self, time):
        """
        Function to convert time from SI to a.u. units
        time: time in seconds
        return: time in a.u.
        """
        return time * 4.1341373336493e16

    def au_to_second(self, time):
        """
        Function to convert time from a.u. to SI units
        time: time in a.u.
        return: time in seconds
        """
        return time / 4.1341373336493e16

    def hartree_to_ev(self, hartree):
        """
        Function to convert energy from a.u. to SI units
        hartree: energy in hartree
        return: energy in electronvolts
        """
        return hartree * 2.72114e1

    def ev_to_hartree(self, ev):
        """
        Function to convert energy from SI to a.u. units
        ev: energy in electronvolts
        return: energy in hartree
        """
        return ev / 2.72114e1

    def au_to_m(self, au):
        return self.a_0_SI * au

    def m_to_au(self, m):
        return m / self.a_0_SI

    def hz_to_au(self, hz):

        return hz * 1.51983e-16

    def au_to_hz(self, au):
        # REF: https://onlinelibrary.wiley.com/doi/pdf/10.1002/3527605606.app9
        return au / 1.51983e-16 #* 4.13413732e16 #2.4e-17 #1.51983e-16 ###################### TODO: atomic unit conversion for frequency
     
    def potential(self, z, time):  # Infinite square well
        """
        Function to define the potential of the lead.
        z: the position in the CNT system
        time: time state evolved to
        return: the potential energy
        """

        if self.potential_type == 0:  # infinite-square well potential
            total_potential = 0  # zero potential inside the scattering region
                                 # outside the potential will be infinity
        elif self.potential_type == 1:  # parabolic potential
            total_potential = .5 * ((z * self.omega_0_au) ** 2)  # parabolic potential inside the scattering region
        
        if self.pertubation_type == "cos":
            self.pertubation = self.cosine_v_ac
        else:
            self.pertubation = self.sine_v_ac
        
        total_potential += self.pertubation(time, z, self.eV_0_au, self.pulse_frequency_au, self.total_length_au)
        
        return total_potential
    
    def import_mumax3_simulations(self):
        """
        Function to import magnetic fields from Mumax3 simulation file
        return: magnetic fields in the x,y and z directions
        """

        file_name = self.magnetic_field_file
        # load mumax3 file 
        f = np.load(r".\{}".format(file_name))

        # Effective magnetic field vectors across the wire
        B_x = f[0, 1, 0:100, 149]
        B_y = f[1, 1, 0:100, 149]
        B_z = f[2, 1, 0:100, 149]

        fig = plt.figure()
        plt.plot(self.z_SI, B_x, label="$B_x$")
        plt.plot(self.z_SI, B_y, label="$B_z$")
        plt.plot(self.z_SI, B_z, label="$B_y$")
        plt.ylabel("Effective Magnetic Field Strength (T)")
        plt.xlabel("$z$ (m)")
        plt.legend();        
        plt.savefig(r"{folder}\{file}_Mumax3_Bfield.png".format(folder = self.folder_path, file = self.file_name))
        print("Plot of Mumax3 magnetic field profile saved.")

        return B_x, B_z, B_y
    

    def kwant_shape(self, site):
        """
        Function to define the shape of the scattering region.
        syst: the system object.
        site: the current site.
        return: the a boolean saying whether the scattering site should be drawn.
        """
        (z,) = site.pos
        return (-self.total_length_au / 2 <= z < self.total_length_au / 2)

    def B_z_au_func(self, z):
        # generate an array of booleans where the current z matches one of the z values in array
        index = np.around(z, 3) == np.around(self.z_au, 3) 
        return self.B_z_au[index]
    def B_y_au_func(self, z):
        index = np.around(z, 3) == np.around(self.z_au, 3)  
        return self.B_y_au[index]
    def B_x_au_func(self, z):
        index = np.around(z, 3) == np.around(self.z_au, 3)  
        return self.B_x_au[index]
    
    def make_system(self):
        """
        Function to create the system
        syst: the system object.
        return: the system object
        """
        if self.potential_type == 1:  # parabolic potential
            self.potential_text = "parabolic"

        else:  # infinite square well
            self.potential_text = "infinite-well"

        # make template for lattice points based on the inputted hamiltonian/
        self.template = kwant.continuum.discretize(self.hamiltonian, grid=self.lattice_size_au) 
        # build the system
        self.syst = kwant.Builder()
        # add the nanotube to the system
        self.syst.fill(self.template, self.kwant_shape, (0,))
        # finalise the system
        self.syst = self.syst.finalized()
        # plot the system
        kwant.plot(self.syst)
        
        if self.magnetic_field_file != "none":  # if the user provides a magnetic field file
            self.B_x_SI, self.B_y_SI, self.B_z_SI = self.import_mumax3_simulations()  # import the magnetic fields
            self.B_x_au = self.tesla_to_au(self.B_x_SI)
            self.B_y_au = self.tesla_to_au(self.B_y_SI)
            self.B_z_au = self.tesla_to_au(self.B_z_SI)
            def B_constant(z):
                return -self.g * self.mu_B_au * self.B_z_au_func(z) * self.hbar_au / 2
            def C_constant(z):
                return -self.g * self.mu_B_au * self.B_x_au_func(z) * self.hbar_au / 2

            def D_constant(z):
                return -self.g * self.mu_B_au * self.B_y_au_func(z) * self.hbar_au / 2

        else:
            self.B_y_au = 0
            self.B_z_au = self.B_0_au
            self.B_x_au = self.b_sl_au
            
            def B_constant(z):
                return  -self.g * self.mu_B_au * self.B_z_au * self.hbar_au / 2

            def C_constant(z):
                return -self.g * self.mu_B_au * self.B_x_au * z * self.hbar_au / 2

            def D_constant(z):
                return -self.g * self.mu_B_au * self.B_y_au * self.hbar_au / 2
    
        # coefficient for the kinetic energy term
        A_constant = self.hbar_au ** 2 / (2 * self.m_au)  

        # import these function and coefficients for use in the full Hamiltonian used to define the system
        self.params = dict(A = A_constant, V = self.potential, B=B_constant, C=C_constant, D=D_constant)

        self.tparams = self.params.copy()  # copy the params array
        self.tparams['time'] = 0  # add another parameter, with the time in au 
        
        # compute the Hamiltonian matrix for this system using the above parameters.
        hamiltonian = self.syst.hamiltonian_submatrix(params=self.tparams)
        # From this Hamiltonian matrix compute the eigenvalues (energies) and eigenvectors (wavefunctions).
        eigenValues, eigenVectors = np.linalg.eig(hamiltonian)

        # sort the eigenvectors and eigenvalues according the ascending eigenvalues.
        idx = eigenValues.argsort()
        self.initial_eigenvalues = eigenValues[idx]
        eigenVectors = eigenVectors[:, idx]

        self.psi_1_init = eigenVectors[:, 0]
        self.psi_2_init = eigenVectors[:, 1]


        # tkwant object representing the spin-up state
        self.spin_up_state = onebody.WaveFunction.from_kwant(syst=self.syst,
                                                        psi_init=self.psi_1_init,
                                                        energy=eigenValues[0],
                                                        params=self.params,
                                                        perturbation_type=onebody.kernels.PerturbationExtractor)

        return self.syst
    
    def eigenstates(self):
        """
        Function to compute the eigenstates of the system.
        syst: the system object.
        return: the sorted eigenvalues and eigenvectors.
        """
        
        # compute the Hamiltonian matrix for this system using the above parameters.
        hamiltonian = self.syst.hamiltonian_submatrix(params=self.tparams)
        # From this Hamiltonian matrix compute the eigenvalues (energies) and eigenvectors (wavefunctions).
        eigenValues, eigenVectors = np.linalg.eig(hamiltonian)
        
        # Sort the eigenvectors and eigenvalues according the ascending eigenvalues.
        idx = eigenValues.argsort()
        eigenValues = eigenValues[idx]
        eigenVectors = eigenVectors[:, idx]

        return eigenValues, eigenVectors

    def initial_pdfs(self):

        """
        Function to plot the initial probability density functions of the spin-up and spin-down ground state.
        syst: the system object.
        return: PDFs of the spin-up and spin-down ground state.
        """
        eigenValues, eigenVectors = self.eigenstates()

        # https://kwant-project.org/doc/dev/tutorial/operators - this explains the output of the eigenvectors.
        psi1 = self.psi_1_init
        psi1_up, psi1_down = psi1[::2], psi1[1::2]
        # even indices give the spin up and odd indices give the spin down states
        density_1 = np.abs(psi1_up) ** 2 + np.abs(psi1_down) ** 2
        psi2 = self.psi_2_init
        psi2_up, psi2_down = psi2[::2], psi2[1::2]
        density_2 = np.abs(psi2_up) ** 2 + np.abs(psi2_down) ** 2

        fig = plt.figure()
        plt.plot(self.z_SI, density_1, label='$|G_-⟩$')
        plt.plot(self.z_SI, density_2, label='$|G_+⟩$')
        plt.ylabel("$|\psi(z)|^2$")
        plt.xlabel("z (m)")
        plt.legend(loc="upper right");

        # save PDF graph 
        plt.savefig(r"{folder}\{file}_InitialPDF.png".format(folder = self.folder_path, file = self.file_name))
        print("Plot of PDFs at t = 0 saved.")

        return density_1, density_2

    def initial_energies(self):
        """
        Function to plot energy levels of the system in meV
        syst: the system object.
        return: the energy eigenvalues of the system in meV
        """

        fig = plt.figure()

        y = self.hartree_to_ev(np.real(self.initial_eigenvalues))*1e3
        print("E_1 is", y[0], "eV.")
        print("E_2 is", y[1], "eV.")

        # Plot the energies of the different levels.
        plt.plot([0, 1], [y[0], y[0]], label=r"$G_-$")
        plt.plot([0, 1], [y[1], y[1]], label=r"$G_+$")
        plt.plot([0, 1], [y[2], y[2]], label=r"$E_-$")
        plt.plot([0, 1], [y[3], y[3]], label=r"$E_+$")
        plt.ylabel("$E$ (meV)")
        plt.legend(loc="best");

        # save initial eigenenergy graph with set file name
        plt.savefig(r"{folder}\{file}_InitialEnergies.png".format(folder = self.folder_path, file = self.file_name))
        print("Plot of eigenenergies at t = 0 saved.")
        return y

    def evolve(self, full_time, time_steps):
        """
        Function to evolve using the time-dependent EDSR Hamiltonian
        full_time (s): full simulation time to evolve the state
        time_steps: number of time steps
        return: PDF, electric field potential, expectation value of the Pauli spin matrices and time for each time step
        """
        self.evolve_state = True
        
        def x_onsite(site):
            """
            Function to compute the position operator matrix for a 'see-sawing' and Gaussian-shaped electric field
            site: position on Kwant grid
            """
            if self.electric_field_type == "gaussian": # for a Gaussian electric field
                return np.exp((-(site.pos[0]-self.z_shift)**2)/(2*self.sigma**2))* np.identity(2)
            
            else: # for a see-sawing electric field
                return [site.pos[0]] * np.identity(2) 

        # spin matrices
        sigma_x = np.array([[0, 1],
                            [1, 0]])
        sigma_y = np.array([[0, 1j],
                            [-1j, 0]])
        sigma_z = np.array([[1, 0],
                            [0, -1]])
        
        # extract lowest two energies - qubit states
        E_1 = np.real(self.initial_eigenvalues[0])
        E_2 = np.real(self.initial_eigenvalues[1])
        
        # compute the difference in these energies
        delta_E = np.abs(E_2 - E_1)
        # compute resonant frequency for the rabi oscillations
        omega_res = delta_E / self.hbar_au
        self.pulse_frequency_au = omega_res  / ( 2 * np.pi) # drive frequency of electric field equal to res. freq.
        
        # define the density operator of x
        rho_x = kwant.operator.Density(self.syst, x_onsite, sum=True)
        # compute this density operator on the ground states <2|x|1>:
        rho_x_2_1 = rho_x(self.psi_2_init, self.psi_1_init)

        # compute the energy E_x
        E_x = np.abs(np.real(2 * self.eV_0_au * rho_x_2_1))
        # compute t_pi
        self.t_pi = 2 * np.pi / E_x
        
        print("Resonant Frequency, f (Hz)", (self.au_to_hz(self.pulse_frequency_au)*2*np.pi))
        print('\u03B5_z (eV)', self.hartree_to_ev(delta_E))
        print('\u03B5_x (eV)', self.hartree_to_ev(E_x))
        print("t_pi (s)", self.au_to_second(self.t_pi))
    
    
        # save the parameters defining the system - can be used in the Qiskit Dynamics simulation
        # in a.u.
        if self.magnetic_field_file != "none":  # if the user provides a magnetic field file
    
            self.B_x_SI, self.B_y_SI, self.B_z_SI = self.import_mumax3_simulations()  # import the magnetic fields
            self.B_x_au = self.tesla_to_au(self.B_x_SI)
            self.B_y_au = self.tesla_to_au(self.B_y_SI)
            self.B_z_au = self.tesla_to_au(self.B_z_SI)          
            
            data = {
            'B_0': self.B_0_au,
            'lattice_points': self.number_of_lattices,
            'length': self.total_length_au,
            'eV_0': self.eV_0_au,
            'E_sl': self.E_sl_au,
            'b_sl': self.b_sl_au,
            'E_x': E_x,
            'E_z': delta_E,
            'pulse_freq': self.pulse_frequency_au,
            'perturbation': self.pertubation_type,
            'potential_type': self.potential_text,
            'effective_mass': self.m_au,
            't_pi': self.t_pi,
            'Real_Bfield': 'True',
            'B_z': list(np.float64(self.B_z_au)),
            'B_y': list(np.float64(self.B_y_au)),
            'B_x': list(np.float64(self.B_x_au))
        }

        else:
            
            data = {
            'B_0': self.B_0_au,
            'lattice_points': self.number_of_lattices,
            'length': self.total_length_au,
            'eV_0': self.eV_0_au,
            'E_sl': self.E_sl_au,
            'b_sl': self.b_sl_au,
            'E_x': E_x,
            'E_z': delta_E,
            'pulse_freq': self.pulse_frequency_au,
            'perturbation': self.pertubation_type,
            'potential_type': self.potential_text,
            'effective_mass': self.m_au,
            't_pi': self.t_pi,
            'Real_Bfield': 'False'
        }
           
        # save parameters
        json_string = json.dumps(data)

        with open(r'{folder}\{file}_QiskitParameters.json'.format(folder = self.folder_path, file = self.file_name), 'w') as outfile:
            outfile.write(json_string)
    
    
        # time array to evolve states in a.u.
        times = self.second_to_au(np.linspace(0, full_time, time_steps))
        
        # initial wavefunction - lowest qubit state
        psi = self.spin_up_state

        # kwant density operator
        density_operator = kwant.operator.Density(self.syst)
        # density operators for pauli spin matrices
        rho_sz = kwant.operator.Density(self.syst, sigma_z, sum=True)
        rho_sy = kwant.operator.Density(self.syst, sigma_y, sum=True)
        rho_sx = kwant.operator.Density(self.syst, sigma_x, sum=True)

        # empty array for data
        densities = []    
        potential = [] 
        exp_z = []
        exp_y = []
        exp_x = []
        
        # loop over each time step
        for time in times:

            print("Evolving state to time", self.au_to_second(time), "s")
            
            psi.evolve(time) # evolve the wavefunction according to the EDSR hamiltonian
            density = psi.evaluate(density_operator) # compute the PDF    
            perturb = self.potential(self.z_au, time) # determine the electric energy
            
            # compute the expectation values of the spin operators on the evolved state
            spin_z = np.real(psi.evaluate(rho_sz))
            spin_y = np.real(psi.evaluate(rho_sy))
            spin_x = np.real(psi.evaluate(rho_sx))
            
            # add values to list
            densities.append(density)
            potential.append(perturb)
            exp_z.append(spin_z)
            exp_y.append(spin_y)
            exp_x.append(spin_x)

        return densities, potential, exp_z, exp_y, exp_x, times
    
    
    def final_graphs(self, times, Bx_SI, spin_x, spin_y, spin_z):
        """Function to plot the effective magnetic field in the x-direction against time and the Rabi oscillations.
        Output:
            times (a.u): time array
            Bx_SI (T): effective magnetic field in the x-direction
            spin_x, spin_y, spin_z: expectation of the Pauli spin matrices for each time step
        return: saved plots"""
        
        plt.figure()
        plt.plot(self.au_to_second(times), Bx_SI*1e3, label  = '$\\langle B_x \\rangle$')
        plt.ylabel("Effective Magnetic Field Strength (mT)")
        plt.xlabel("Time (s)")
        plt.legend();
        
        plt.savefig(r"{folder}\{file}_Bxfield.png".format(folder = self.folder_path, file = self.file_name))
        print("Plot of effective magnetic field against time saved.")

        plt.figure()
        plt.plot(self.au_to_second(times), spin_x, label ='$\\langle X \\rangle$')
        plt.plot(self.au_to_second(times), spin_y, label ='$\\langle Y \\rangle$')
        plt.plot(self.au_to_second(times), spin_z, label ='$\\langle Z \\rangle$')
        plt.ylabel("Expectation Value")
        plt.xlabel("Time (s)")
        plt.legend();
        
        plt.savefig(r"{folder}\{file}_RabiOsc.png".format(folder = self.folder_path, file = self.file_name))
        print("Plot of Rabi oscillations saved.")
        
        return True
        
     
    def density_potential_animation(self, density, y3):
        """
        Function to visualise the dynamics of the PDF and electric field potential over time.
        Inputs: 
            density: PDF values
            y3 (a.u.): electric field pertubation 
            file_name: name of file to save animation
        return: animation of PDF and electric field potential
        """
        
        # create a figure with two subplots
        fig_animation, axs = plt.subplots(2,1)

        # 1920 x 1080
        w_in_inches = 10
        h_in_inches = 6
        dpi = 100
        fig_animation.set_size_inches(w_in_inches, h_in_inches, True)
        fig_animation.set_dpi(dpi)
        fig_animation.tight_layout(pad=5.0)

        # set variables for each of the axes
        ax1 = axs[0]
        ax2 = axs[1]

        # initialise two line objects (one in each axes)
        line1, = ax1.plot([], [], lw=2)
        line2, = ax2.plot([], [], lw=2, color='r', label='$H_1(t)$')
        line = [line1, line2]
        
        z_max = self.total_length_SI / 2
        y2_max = self.eV_0_eV*1e3
        y1 = density[0]
        y1_max = np.max(y1) * 2

        # PDF
        ax1.set_xlim(-z_max, z_max)
        ax1.set_ylim(0, np.max(density[-1])*3)
        ax1.set_ylabel("$|\psi(z,t)|^2$")
        ax1.set_xlabel("z ($m$)")
        ax1.set_title("PDF over time")

        # Potential
        ax2.set_xlim(-z_max, z_max)
        ax2.set_ylim(-y2_max, y2_max)
        ax2.set_ylabel("$E$ (meV)")
        ax2.set_xlabel("z ($m$)")
        ax2.set_title("Electric field potential over time")
        
        # draw gate on graph
        gate = patches.Rectangle((-220e-9,-y2_max), 100e-9, y2_max/2, fill = True, color = '#929591')
        ax2.add_patch(gate)
        plt.text(-185e-9, -y2_max/1.3, 'Gate')
        # draw magnet on graph
        magnet = patches.Rectangle((-100e-9,-y2_max), 300e-9, y2_max/2, fill = True,color = '#929591')
        ax2.add_patch(magnet)
        plt.text(30e-9, -y2_max/1.3, 'Magnet')

        # initialisation function: plot the background of each frame
        def init():
            line[0].set_data([], [])
            line[1].set_data([], [])
            return line

        def animate(i):
            # update line objects.
            line[0].set_data(self.z_SI, density[i])
            line[1].set_data(self.z_SI, self.hartree_to_ev(y3[i])*1e3)
            ax2.legend(loc="upper right")
            return line

        # call the animator.  blit=True means only re-draw the parts that have changed.
        anim = animation.FuncAnimation(fig_animation, animate, init_func=init,
                                       frames=len(density), interval=400, blit=True)

        anim.save(r"{folder}\{file}_density_potential_animation.mp4".format(folder = self.folder_path, file = self.file_name), writer='ffmpeg')
        plt.show()

        return True
    
    def full_system_animation(self,density, y3, Bfield, spinx, spiny, spinz, times):
        """
        Function to visualise the dynamics of the PDF, electric field potential, effective magnetic field and Rabi oscillations
        over time.
        Inputs: 
            density: PDF values
            y3 (a.u.): electric field pertubation 
            Bfield (T): effective magnetic field in the x-direction
            spinx, spiny, spinz: expectation of the Pauli spin matrices for each time step
            times (a.u.): time of state evolution
            file_name: name of file to save animation
        return: full animation of the EDSR system
        """
        # convert times into S.I.
        times = self.au_to_second(times)
        
         # create a figure with two subplots
        fig_animation, axs = plt.subplots(2, 2)

        # 1920 x 1080
        w_in_inches = 19.2
        h_in_inches = 10.8
        dpi = 100
        fig_animation.set_size_inches(w_in_inches, h_in_inches, True)
        fig_animation.set_dpi(dpi)
        fig_animation.tight_layout(pad=5.0)  # add padding to subplots.

        # set variables for each of the axes.
        ax1 = axs[0, 0]
        ax2 = axs[1, 0]
        ax3 = axs[0, 1]
        ax4 = axs[1, 1]

        # initialise two line objects (one in each axes)
        line1, = ax1.plot([], [], lw=2)
        line2, = ax2.plot([], [], lw=2, color='r', label='$H_1(t)$')
        line3, = ax3.plot([], [], lw=2, label=r'$\langle X \rangle$')
        line4, = ax3.plot([], [], lw=2, label=r'$\langle Y \rangle$')
        line5, = ax3.plot([], [], lw=2, label=r'$\langle Z \rangle$')
        line6, = ax4.plot([], [], lw=2, label=r'$\langle B_x \rangle$')
        line = [line1, line2, line3, line4, line5, line6]

        z_max = system.total_length_SI / 2
        y2_max = system.eV_0_eV*1e3
        y1 = density[0]
        y1_max = np.max(y1) * 2

        # PDF
        ax1.set_xlim(-z_max, z_max)
        ax1.set_ylim(0, np.max(density[-1])*3)
        ax1.set_ylabel("$|\psi(z,t)|^2$")
        ax1.set_xlabel("z ($m$)")
        ax1.set_title("PDF over time")

        # Potential
        ax2.set_xlim(-z_max, z_max)
        ax2.set_ylim(-y2_max, y2_max)
        ax2.set_ylabel("$E$ (meV)")
        ax2.set_xlabel("z ($m$)")
        ax2.set_title("Electric field potential over time")

        # Expectation values
        ax3.set_xlim(0, times[-1])
        ax3.set_ylim(-1.0, 1.0)
        ax3.set_ylabel("Expectation Value")
        ax3.set_xlabel("t (s)")
        ax3.set_title("Rabi Oscillation plot")

        # Magnetic Field
        ax4.set_xlim(0, times[-1])
        ax4.set_ylim(-1.3 * np.max(Bfield)*1e3, 1.1 * np.max(Bfield)*1e3)
        ax4.set_ylabel("Magnetic Field Strength (mT)")
        ax4.set_xlabel("t (s)")
        ax4.set_title("Effective magnetic field in the x-direction against time")

        # draw gate on graph
        gate = patches.Rectangle((-220e-9,-y2_max), 100e-9, y2_max/2, fill = True, color = '#929591')
        ax2.add_patch(gate)
        ax2.text(-185e-9, -y2_max/1.3, 'Gate')
        # draw magnet on graph
        magnet = patches.Rectangle((-100e-9,-y2_max), 300e-9, y2_max/2, fill = True,color = '#929591') 
        ax2.add_patch(magnet)
        ax2.text(30e-9, -y2_max/1.3, 'Magnet')

        # initialisation function: plot the background of each frame
        def init():
            line[0].set_data([], [])
            line[1].set_data([], [])
            line[2].set_data([], [])
            line[3].set_data([], [])
            line[4].set_data([], [])
            line[5].set_data([], [])
            return line

        def animate(i):
            # update line objects.
            line[0].set_data(self.z_SI, density[i])
            line[1].set_data(self.z_SI, system.hartree_to_ev(y3[i])*1e3)
            line[2].set_data(times[0:i + 1], spinx[0:i + 1])
            line[3].set_data(times[0:i + 1], spiny[0:i + 1])
            line[4].set_data(times[0:i + 1], spinz[0:i + 1])
            line[5].set_data(times[0:i + 1], Bfield[0:i + 1]*1e3)

            ax2.legend(loc="upper right")
            ax3.legend(loc="upper right")
            ax4.legend(loc="upper right")
            return line

        # call the animator.  blit=True means only re-draw the parts that have changed.
        anim = animation.FuncAnimation(fig_animation, animate, init_func=init,
                                       frames=len(density), interval=400, blit=True)

        anim.save(r"{folder}\{file}_full_system_animation.mp4".format(folder = self.folder_path, file = self.file_name), writer='ffmpeg')
        plt.show()

        return True
        
def main():
    file_name = 'EDSR_Simulation'
    folder_path = r"."

    # set up the system
    lattices = 20 # number of lattice points
    full_time = 1e-12 # full simulation time in s
    steps = 2 # number of time steps

    potential = 0  # infinite square-well potential
    # input into the system to include a realistic magentic field
    magnetic_field_file = "B_eff000000.npy"
    # type of electric field
    electric_field_type = "gaussian" #"see-saw"

    system = System("((A * k_z**2) + V(z, time)) * identity(2) + B(z) * sigma_z + C(z) * sigma_x + D(z) * sigma_y",
                    pertubation_type="sin", electric_field_type=electric_field_type, number_of_lattices=lattices,
                    potential_type=potential, file_name=file_name, folder_path=folder_path)#, magnetic_field_file = magnetic_field_file)

    system.make_system();
    system.initial_pdfs();
    system.initial_energies();

    # evolve the state - note that potential and times are given in a.u. units
    density, potential, spin_z, spin_y, spin_x, times = system.evolve(full_time, steps);
    
    # calculate effective magnetic field using density
    Bx_SI = []
    for i in range(0,len(density)):
        Bx_SI = np.append(Bx_SI, np.sum(density[i]*system.z_SI*system.b_sl_SI))
    
    # save data
    data = {'density': np.array(density).tolist(), 'potential': np.array(potential).tolist(), 'Bx_SI': np.array(Bx_SI).tolist(), 'spin_z': np.array(spin_z).tolist(),
            'spin_y': np.array(spin_y).tolist(), 'spin_x': np.array(spin_x).tolist(), 'times': np.array(times).tolist()}
    json_string = json.dumps(data)

    with open(r'{folder}\{file}.json'.format(folder = folder_path, file = file_name), 'w') as outfile:
        outfile.write(json_string)

    # display the effective magnetic field graphs and rabi oscillation graphs
    system.final_graphs(times, Bx_SI, spin_x, spin_y, spin_z);
    
    # display the animations of the state
    system.density_potential_animation(density, np.array(potential));
    system.full_system_animation(density, potential, Bx_SI, spin_x, spin_y, spin_z, times);

if __name__ == '__main__':
    main()