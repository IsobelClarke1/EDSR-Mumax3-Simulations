import numpy as np
import matplotlib.pyplot as plt
from numpy import vectorize
from matplotlib import animation
import os.path
import json
import time as timelib
from IPython.display import Image
from csv import writer
from matplotlib import patches

def import_mumax_file(file_loc, file_name):
    """
    Function to import magnetic fields values from Mumax3 simulation file
    Inputs:
        file_loc: location of saved Mumax3 file
        file_name: name of saved Mumax3 file
    """
    f = np.load(file_loc.format(file_name))
    return f

def plot_Bfield(B_x, B_y, B_z, number_of_cells, twofields, folder_path):
    """
    Function to plot the effective magnetic field in the x, y, z direction. 
    Input:
        B_x, B_y, B_z: x, y and z-direction components
        number_of_cells: number of effective magnetic field values 
        twofields:  Boolean (True if 2 external fields applied or False if only one applied)
        folder_path: location to save the final plot of the magnetic field profile
    """
    # define dot region 
    total_length_SI = 0.66e-6 # length of dot
    lattice_size_SI = total_length_SI / number_of_cells
    y_SI = np.arange(-total_length_SI / 2, total_length_SI / 2, lattice_size_SI, dtype=np.double) 
    
    fig, ax = plt.subplots()
    plt.plot(y_SI, B_x, label="$B_x$")
    plt.plot(y_SI, B_y, label="$B_z$")
    plt.plot(y_SI, B_z, label="$B_y$")
    plt.ylabel("Effective Magnetic Field Strength (T)")
    plt.xlabel("$z$ (m)")
    plt.legend(loc="upper left", bbox_to_anchor=(0.8,0.9))
    
    # add magnet to diagram - position varies if a second, smaller magnetic field is applied
    if twofields == True: 
        magnet = patches.Rectangle((-100e-9,-0.005), 300e-9, 0.002, fill = True,color = '#929591')
        ax.add_patch(magnet)
        plt.text(10e-9,-0.004, 'magnet')
    else: 
        magnet = patches.Rectangle((-100e-9,-0.05), 300e-9, 0.03, fill = True,color = '#929591')
        ax.add_patch(magnet)
        plt.text(10e-9,-0.040, 'magnet')

    # add source and drain labels
    plt.text(- total_length_SI / 2 - 25e-9 , 0.105, 'S')
    plt.text(total_length_SI / 2 + 2e-9 , 0.105, 'D')
    
    # vertical lines at source and drain
    plt.axvline(x = -330e-9, color = '#808080')
    plt.axvline(x = 325e-9, color = '#808080')  

    plt.savefig(r"{folder}\Mumax3_MagneticFieldProfile.png".format(folder = folder_path))


def magnetic_field(x_disp, y_disp, angle, number_of_cells, file_loc):
    """
    Function to obtain the effective magnetic field in the x, y, z direction for a central, rotated or translated CNT.
    Input:
        x_disp (m): displacement in x direction from the centre position
        y_disp (m): displacement in y direction in diagrams corresponds to the z direction in Mumax3 
        angle (degrees): angle of rotation
        file_loc: location of saved Mumax3 file
    """
    
    global tube_height, cx, cy, cz, xGrid

    # convert displacements to number of cells
    x_disp = int(x_disp/cx)
    z_disp = int(y_disp/cz) # in mumax3 the y-direction corresponds to the z-direction project diagrams

    # import data
    data = import_mumax_file(file_loc, "B_eff000000.npy")
    
    # determine center of mumax grid
    gridCenter = (xGrid/2) - 1
    
    # non-rotated CNT
    if angle == 0:
        x_value = int(gridCenter + x_disp)
        B_x = data[0, z_disp, 0:number_of_cells, x_value] 
        B_y = data[1, z_disp, 0:number_of_cells, x_value]
        B_z = data[2, z_disp, 0:number_of_cells, x_value]
     
    # rotated CNT       
    else:
        
        # distance along x-direction
        opp = (tube_height/2)* np.tan(angle*(np.pi/180)) # np.sin is in radians
        # first x cell
        x_cells_min = opp/cx
        y_cells = 100 # distance between electrodes in y-direction

        # start point of x array
        x_start = (gridCenter + 1 + x_disp) - x_cells_min
       # new length along x axis
        new_opp = (100*cy)* np.tan(angle*(np.pi/180))
        x_cells_diff = new_opp/cx
        # end point of x array
        x_end = (x_start + x_cells_diff)

        y_start = 0
        y_end = y_cells
        x_start = int(x_start)
        x_end = int(x_end)

        # array of x and y cells
        y_array = np.arange(y_start, y_end, 1)
        x_array = np.arange(x_start, x_end+2, 1)
        
        # number of y values at each x point
        factor = int(y_cells/x_cells_diff)
        array = np.arange(0, y_end, factor)
        
        # empty array for magnetic field values
        B_x = np.zeros((y_end,))
        B_y = np.zeros((y_end,))
        B_z = np.zeros((y_end,))

        factor_array = np.arange(0,factor,1)

        for i, j in zip(array, y_array):
            
            for x in factor_array:
                print(x)          

                B_x[i + x,] = data[0, z_disp, y_array[i+x],x_array[j]]
                B_y[i + x,] = data[1, z_disp, y_array[i+x],x_array[j]]
                B_z[i + x,] = data[2, z_disp, y_array[i+x],x_array[j]]
                
                
                print((0,z_disp,y_array[i+x],x_array[j]))
        
            # if x_array finishes before y_array - not full diagonal
            if x_array[j] == x_end+1 and y_array[i+x] != y_end-1:
                print("Diagonal terminated early. Not equivalent number of cells in x and y direction.")
                break
            
    return B_x, B_y, B_z

## define parameters
# grid size in number of cells
xGrid = 300
yGrid = 276
zGrid = 64

# cell size in meters
cx = 5e-9
cy = 66e-10
cz = 5e-9
#tube height in meters
tube_height = 1820e-9

file_loc = r'C:\Users\bella\OneDrive\Documents\Year 4\MSci Project\Coding\My Code\Magnetic_Field_Profile.out\{}'
folder_save = r'C:\Users\bella\OneDrive\Documents\Year 4\MSci Project\Tkwant Data WORKING'

# one cell in the z direction
# note: the the z-displacement in Mumax3 is the y-direction in diagrams in this project
B_x, B_y, B_z = magnetic_field(0, cz, 10, 100, file_loc)

plot_Bfield(B_x, B_y, B_z, len(B_x), False, folder_save)
