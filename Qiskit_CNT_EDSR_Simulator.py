import matplotlib.pyplot as plt
import numpy as np
from qiskit.quantum_info import Operator
from qiskit_dynamics import Solver, Signal
from qiskit_dynamics.signals import SignalSum
from qiskit.quantum_info.states import Statevector
from qiskit.quantum_info import DensityMatrix
import qutip
from matplotlib import pyplot, animation
from mpl_toolkits.mplot3d import Axes3D
import json

def load_parameters(folder_path ,file_name):
    """Function to obtain parameters for the effective qubit Hamiltonian from Kwant simulation
    Input:
        folder_path: location of saved Qiskit parameters
        file_name: file name of saved Qiskit parameters"""
    with open(r'{folder}\{file}.json'.format(folder = folder_path, file = file_name)) as json_file:
            parameters = json.load(json_file)
    return parameters

folder_path = r'.'
file_name = 'QiskitParameters'
parameters = load_parameters(folder_path, file_name)

## conversion functions
def ev_to_hartree(ev):
    return ev / 2.72114e1

def second_to_au(time):
    return time * 4.1341373336493e16

def au_to_second(time):
    return time / 4.1341373336493e16

def hartree_to_ev(hartree):
    return hartree * 2.72114e1

def au_to_tesla(au):
    return au * 2.35e5

def au_to_m(au):
    a_0_SI = 5.2917721090380e-11
    return a_0_SI * au

def m_to_au(m):
    a_0_SI = 5.2917721090380e-11
    return m / a_0_SI

def hz_to_au(hz):
    return hz / 4.13413732e16

def tesla_to_au(tesla):
    return tesla / 2.35e5

def au_to_hz(au):
    return au * 4.13413732e16 

# constants in a.u.
hbar_au = 1
e_au = 1
mu_B_au = .5
e_SI = 1.602176634e-19 # in J
mu_B_SI = 9.2740100783e-24 # in J/T

# import parameters from Kwant data
ex = parameters['E_x']
ez = parameters['E_z']
pulse_freq = parameters['pulse_freq']
t_pi = parameters['t_pi']
perturb_type = parameters['perturbation']
m_au = parameters['effective_mass']
number_of_lattices = parameters['lattice_points']
total_length_au = parameters['length']

# system parameters
B_0 = parameters['B_0']
E_sl = parameters['E_sl']
b_sl = parameters['b_sl']
eV = parameters['eV_0']
n_steps = 50

## output parameters for hamiltonian
print("frequency (Hz):", au_to_hz(pulse_freq))
print("ex (eV):", hartree_to_ev(ex))
print("ez (eV):", hartree_to_ev(ez))
print("\t")
print("B_0 (T):", au_to_tesla(B_0))
print("E_sl (eV):", hartree_to_ev(E_sl))
print("eV0 (eV):", hartree_to_ev(eV))
print("effective mass (au):", m_au)
print("\t")
print("t_pi (s):", au_to_second(t_pi))
print("number of lattices:", number_of_lattices)

if parameters['Real_Bfield'] == 'True':
    B_x = au_to_tesla(np.array(parameters['B_x']))
    B_y = au_to_tesla(np.array(parameters['B_y']))
    B_z = au_to_tesla(np.array(parameters['B_z']))


    plt.figure()
    plt.plot(z_SI, B_x, label="$B_x$")
    plt.plot(z_SI, B_z, label="$B_z$")
    plt.plot(z_SI, B_y, label="$B_y$")
    plt.ylabel("Effective Magnetic Field Strength (T)")
    plt.xlabel("$z$ (m)")
    plt.legend(loc = 'right');

## setup solver with Hamiltonian model
# pauli spin matrices
X = Operator.from_label('X')
Y = Operator.from_label('Y')
Z = Operator.from_label('Z')

# solve with the static hamiltonian and time-dependent part
solver = Solver(static_hamiltonian = 0.5 * ez * Z,
                hamiltonian_operators = [0.5 * ex * X]) 

## solve the system
t_final = t_pi # evolve up to t_pi
y0 = Statevector([1,0]) # initial state vector is 0
t_eval = np.linspace(0., t_final, n_steps) # array of time values to evolve over

# type of sinusodial perturbation
if perturb_type == 'cos':
    enve = 1
else:
    enve = -1j
# define the perturbation to the state  
define_signal = Signal(envelope= -1j, carrier_freq = pulse_freq)
signals = [define_signal]
#solve the Hamiltonian
sol = solver.solve(t_span=[0., t_final], y0=y0, signals=signals, t_eval=t_eval, atol=1e-8, rtol=1e-8)

def plot_qubit_dynamics(sol, t_final, X, Y, Z, n_steps, folder_path, file_name):
    """Function to plot the Rabi oscillations of the evolved state.
    Inputs:
        t_final (s): final time state evolved up to
        X, Y, Z: Pauli spin matrix operators
        n_steps: number of time steps
        folder_path: location of saved graph
        file_name: name of saved graph
    return: expectation values for the Pauli spin matrices for each time step
    """
    fontsize = 16
    t_eval = np.linspace(0., t_final, n_steps)
    n_times = len(sol.y)
    x_data = np.zeros((n_times,))
    y_data = np.zeros((n_times,))
    z_data = np.zeros((n_times,))

    for t_i, sol_t in enumerate(sol.y):
        x_data[t_i] = sol_t.expectation_value(X).real
        y_data[t_i] = sol_t.expectation_value(Y).real
        z_data[t_i] = sol_t.expectation_value(Z).real
 
    _, ax = plt.subplots(figsize = (10, 6))
    plt.rcParams.update({'font.size': fontsize})
    plt.plot(t_eval, x_data, label = '$\\langle X \\rangle$')
    plt.plot(t_eval, y_data, label = '$\\langle Y \\rangle$')
    plt.plot(t_eval, z_data, label = '$\\langle Z \\rangle$')
    plt.legend(fontsize = fontsize)
    ax.set_ylabel("Expectation Value")
    ax.set_xlabel('Time', fontsize = fontsize)
    plt.show()
    
    plt.savefig(r"{folder}\{file}_QiskitRabiOsc.png".format(folder = folder_path, file = file_name))
    
    return x_data, y_data, z_data

def RabiOsc_animation(x_data, y_data, z_data, t_final, folder_path):
    """Function to animate the Rabi oscillations over time
    Inputs:
        x_data, y_data, z_data: the expectation values of the Pauli spin matrices for each time step
        t_final (s): final time state evolved up to
        folder_path: location of saved animation
    """
    t_eval_SI = au_to_second(np.linspace(0., t_final, n_steps))

    # create a figure with two subplots
    fig_animation, axs = plt.subplots(1,1)

    # 1920 x 1080
    w_in_inches = 10 #19.2
    h_in_inches = 6.4#10.8
    dpi = 100
    fig_animation.set_size_inches(w_in_inches, h_in_inches, True)
    fig_animation.set_dpi(dpi)

    ax1 = axs
    # initialise line objects for axees
    line1, = ax1.plot([], [], lw=2, label=r'$\langle X \rangle$')
    line2, = ax1.plot([], [], lw=2, label=r'$\langle Y \rangle$')
    line3, = ax1.plot([], [], lw=2, label=r'$\langle Z \rangle$')
    line = [line1, line2, line3]

    # Expectations
    ax1.set_xlim(0, t_eval_SI[-1])
    ax1.set_ylim(-1.0, 1.0)
    ax1.set_ylabel("Expectation Value")
    ax1.set_xlabel("t (s)")
    ax1.legend(loc="upper right")

    # initialisation function to plot the background of each frame
    def init():
        line[0].set_data([], [])
        line[1].set_data([], [])
        line[2].set_data([], [])
        return line

    def animate(i):
        line[0].set_data(t_eval_SI[0:i + 1], x_data[0:i + 1])
        line[1].set_data(t_eval_SI[0:i + 1], y_data[0:i + 1])
        line[2].set_data(t_eval_SI[0:i + 1], z_data[0:i + 1])
        return line
       
    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig_animation, animate, init_func=init,
                                   frames=len(t_eval), interval=80, blit=False)

    anim.save('{}/QiskitRabiOsc_animation.mp4'.format(folder_path), writer='ffmpeg')
    plt.show()

def BlochSphere_animation(x_data, y_data, z_data,folder_path):
    """Function to animate the Rabi oscillations over time on a Bloch sphere
    Inputs:
        x_data, y_data, z_data: the expectation values of the Pauli spin matrices for each time step
        folder_path: location of saved animation
    """
    fig = plt.figure();
    ax = Axes3D(fig, azim=-40, elev=30, auto_add_to_figure=False)
    sphere = qutip.Bloch(axes=ax)
    sphere.make_sphere()
    fig.add_axes(ax) 

    def animate(i):
        # add vector and point for the X, Y and Z expectation values at a given time step
        sphere.clear() # clear previous data points
        sphere.add_vectors([x_data[i], y_data[i], z_data[i]])
        sphere.add_points([x_data[i], y_data[i], z_data[i]])
        sphere.make_sphere()
        return ax

    def init():
        # initalise sphere
        sphere.vector_color = ['r']
        return ax

    ani = animation.FuncAnimation(fig, animate, init_func=init, frames=len(t_eval), interval=80, blit=True, repeat=False);
    ani.save('{}/QiskitBlochSphere_animation.mp4'.format(folder_path), writer='ffmpeg')
    plt.show()

# plot the rabi oscillations 
x_data, y_data, z_data = plot_qubit_dynamics(sol, au_to_second(t_final), X, Y, Z, n_steps, folder_path, file_name)
# obtain an animation of the Rabi oscillation graph
RabiOsc_animation(x_data, y_data, z_data, t_final, folder_path)
# obtain an animation of the Rabi oscillations on a Bloch sphere
BlochSphere_animation(x_data, y_data, z_data,folder_path)