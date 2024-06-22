import math
import time
import numpy as np
import pandas as pd
import scipy.optimize
from scipy import integrate
import xlsxwriter
import matplotlib.pyplot as plt
from plotly.subplots import make_subplots
import plotly.graph_objs as go


def Inverse_Kinetostatic_traj(p_ee, angles, init_guess, RPY=False):

    '''input output function
    input arguments as follows: 
    pos_EE = position of the EE--> (3,) shape,
    RPY = vector containing roll, pitch, yaw values--> (3,) shape
    init_guess = intial guess vector of size 36--> [nx1,ny1,nz1, mx1, my1, mz1,...,mxL, myL, mzL, L1,...L6] 
    plot_bool = for activating the plotting function
    saveData_bool = for activating the saving data as excel file


    output vector contains following values: 
    n = force vector at the cross-section s=L (last cross-section at EE)-->3D vector 
    m = moment vector at the cross-section s=L (last cross-section at EE)-->3D vector 
    L = Optimised length vector (Linear actuator)-->6D vector '''
    
    def rotation_matrix_3d(theta_x, theta_y, theta_z): #angles are in rad
        '''create a rotation matrix given angles about x,y,z axes'''
        cos_x = np.cos(theta_x)
        sin_x = np.sin(theta_x)
        Rx = np.array([[1, 0, 0],
                    [0, cos_x, -sin_x],
                    [0, sin_x, cos_x]])
        cos_y = np.cos(theta_y)
        sin_y = np.sin(theta_y)
        Ry = np.array([[cos_y, 0, sin_y],
                    [0, 1, 0],
                    [-sin_y, 0, cos_y]])
        cos_z = np.cos(theta_z)
        sin_z = np.sin(theta_z)
        Rz = np.array([[cos_z, -sin_z, 0],
                    [sin_z, cos_z, 0],
                    [0, 0, 1]])
        R = Rz.dot(Ry.dot(Rx))
        return R

    
    def RPY_to_matrix(roll, pitch, yaw):
        '''Map roll, pitch, yaw to matrix'''
        matrix = np.array([[np.cos(yaw)*np.cos(pitch),-np.sin(yaw)*np.cos(roll)+np.cos(yaw)*np.sin(pitch)*np.sin(roll),np.cos(yaw)*np.sin(pitch)*np.cos(roll)+np.sin(yaw)*np.sin(roll)],
                [np.sin(yaw)*np.cos(pitch), np.cos(yaw)*np.cos(roll)+np.sin(yaw)*np.sin(pitch)*np.sin(roll),np.sin(yaw)*np.sin(pitch)*np.cos(roll)-np.cos(yaw)*np.sin(roll)],
                [-np.sin(pitch),np.cos(pitch)*np.sin(roll),np.cos(pitch)*np.cos(roll)]])
        return matrix     

 #####################################################################################################################
    '''Defining the rigid kinematics for the 6-RUS parallel continuum robot'''  
    #initialisation of the parameters of the geometry for the end-effector and base-platform
    t = 0.29711 		
    b = 0.21841     	     
    d = 0.05374     	
    D = 0.0675     
    l1 = 0.23164	#Length of the rigid crank	
    L = 0.53		#Total length of the flexible link at a stress free height of 0.4m
    Dist = np.array([D,0,0])
    dist = np.array([[d],[0],[0]])
    #120 degree separation between the motor pairs at the base and spherical joint pairs at the end-effector
    Rz_120 = rotation_matrix_3d(np.deg2rad(0),np.deg2rad(0),np.deg2rad(120))
  
    '''Motor position and orientation in world frame'''
    B = np.array([0,(2/3)*b*np.cos(30*np.pi/180),0])

	#Motor positions M1, M2,...,M6 in world coordinates
    M1 = B + Dist
    M2 = B - Dist
    M3 = Rz_120 @ M1 #B2 + Rz120.dot(Dist)
    M4 = Rz_120 @ M2#B2 - Rz120.dot(Dist)
    M5 = Rz_120.T @ M1
    M6 = Rz_120.T @ M2

    #Concatenate motor position into a single vector
    Mm = np.array([M1,M2,M3,M4,M5,M6])

    #Concatenate corresponding motor orientation into a single vector
    Rm = np.array([np.eye(3),np.eye(3),Rz_120,Rz_120,Rz_120.T,Rz_120.T])

    '''spherical joint position and orientation in world frame'''
    T = np.array([[0],[-(2/3)*t*np.cos(30*np.pi/180)],[0]]) 
    #Spherical joints positions C1, C2,...,C6 in EE coordinates
    C5 = T + dist
    C4 = T - dist
    C6 = Rz_120 @ C4	#T2 - Rz120 .dot(dist)
    C1 = Rz_120 @ C5	#T2 + Rz120 .dot(dist)
    C2 = Rz_120.T @ C4	#T3 - Rz120_ .dot(dist)
    C3 = Rz_120.T @ C5	#T3 + Rz120_.dot(dist)

    Ree = RPY_to_matrix(angles[0], angles[1], angles[2]) #orientation of the EE platform
    #Homogeneous transformation matrix--> EE to base frame coordinate
    O_Tee_EE = np.hstack((np.vstack((Ree, np.zeros((1,3)))), np.vstack(((np.reshape(p_ee,(3,1))), 1))))

    #each column of the matrix below is the spherical joint in world coordinates
    World_position_SphericalJoints =  O_Tee_EE @ (np.hstack((np.vstack((C1, 1)),
                                    np.vstack((C2, 1)),np.vstack((C3, 1)),np.vstack((C4, 1)),
                                    np.vstack((C5, 1)),np.vstack((C6, 1)))))

    #extracting each column from the above matrix which contains spherical joint positions in world
    C1 = World_position_SphericalJoints[:,0][:-1]
    C2 = World_position_SphericalJoints[:,1][:-1]
    C3 = World_position_SphericalJoints[:,2][:-1]
    C4 = World_position_SphericalJoints[:,3][:-1]
    C5 = World_position_SphericalJoints[:,4][:-1]
    C6 = World_position_SphericalJoints[:,5][:-1]

    #orientation of the spherical joint wrt world
    RC = [Ree @ Rz_120, Ree @ Rz_120.T, Ree @ Rz_120.T, Ree @ np.eye(3), Ree @ np.eye(3), Ree @ Rz_120]

    #Concatenate spherical joint positions into a single vector
    r = np.array([C1, C2, C3, C4, C5, C6])

#############################################################################################################
    '''Defining the flexible rod material properties for the 6-RUS parallel continuum robot''' 
    # Independent Parameters for the flexible rod

    E = 110.3e9                 # Young's modulus (N/m2)
    nu = 0.31                   #Poisson's ratio
    G = E/(2*(1+nu))            # Shear modulus (N/m2)
    rho = 4428.8                # Density (kg/m3)
    rad = 0.002                 #radius of rod cross-section
    # h = 0.0012 #thickness of the rectangular cross-section (m)
    # w = 0.077  #width of the cross-section (m)
    ee_mass = 0.5             #mass of the end-effector platform (Kg)
    g = np.array([0,0,-9.81])   #gravity in -z direction

    # Dependent Parameters for the flexible rod
    A = math.pi * rad ** 2      #Cross-sectional area
    I = math.pi * rad ** 4 / 4  # Area moment of inertia
    J = 2 * I                   # Polar moment of inertia

    Kse = np.diag([G * A, G * A, E * A])    #stiffness matrix for shear and extension
    KseInv = np.linalg.inv(Kse)
    Kbt = np.diag([E * I, E * I, G * J])    #stiffness matrix for bending and torsion
    KbtInv = np.linalg.inv(Kbt)

    F = np.array([0, 0, -9.81*ee_mass]) #External force acting at the center of the end-effector
    M = np.zeros(3)                     #External moment acting at the center of the end-effector

####################################################################################################################
    '''Formulation of the shooting method'''
    counter = 0     #counter variable for the optimisation loop

    def residual_function(guess):
        '''Residual vector function which includes the constraints that need to be zero for a solution'''
        # Update guessed initial conditions
        nonlocal counter #count the function calls
        counter += 1
        # print(counter)
        
        residual = np.empty(42) #vector to concatenate the constraint equations
        EF = F
        EM = M

        #loop over each flexible rod to compute the residual vector
        for i in range(6):
            #unpack the unknown variable vector for each rod
            n0 = guess[4*i:4*i+3]       #n_x(0), n_y(0), n_z(0)
            m0 = guess[4*i+3]           #m_z(0)
            q = guess[3*i+24:3*i+27]    #actuator variables + universal joint angles

            #position of the base of the flexible rod in motor coordinate frame
            p0_local = np.array([0, l1*np.cos(q[0]), l1*np.sin(q[0])])

            #position vector of the base of the flexible rod in world coordinate frame
            p0 = Mm[i] + (Rm[i] @ p0_local)

            #orientation at the base of the flexible rod in world coordinate frame
            #orientation of the universal joint in local joint frame
            Rx = rotation_matrix_3d(q[2], 0, 0) #theta_x, theta_y, theta_z
            Ry = rotation_matrix_3d(0, q[1], 0) #theta_x, theta_y, theta_z
            
            #orientation at the base of the flexible rod in world coordinate frame 
            R0 = Rm[i] @ Ry @ Rx
            
            '''Solve the IVP problem for each rod'''
            #concatenate the flexible rod states at the base
            y0 = np.concatenate([p0, np.reshape(R0, 9), n0,[0,0], [m0]])
            nonlocal Y
            num_steps = 100  # Total number of sampled integration steps for uniformity
            s_eval = np.linspace(0, L, num_steps) #samples along the rod length

            #integrate to solve the IVP problem till the rod tip
            #RK45' (default): Explicit Runge-Kutta method of order 5(4)
            #solution of the IVP:
            Y = integrate.solve_ivp(rod_ode, (0.0, L), y0, t_eval = s_eval, atol=1e-5, rtol=1e-5).y
            
        
            #residual constraints at the distal end including the force balance
            '''extract the solution of the IVP for each rod'''
            pL_shot = Y[0:3, -1]        #extract the position vector of the tip of the rod 
            nL = Y[12:15, -1]           #extract the force vector at the tip of the rod 
            mL = Y[15:18, -1]           #extract the moment at the tip of the rod 
            # print(f"pL_shot: {R_ee.shape}")
            residual[6*i:6*i+3] = pL_shot - r[i]    #residual for the position error
            residual[6*i+3:6*i+6] = Ree.T @ mL      #residual for the moment error for spherical joint
            
            #F and M respectively = 0
            EF = EF - nL
            EM = EM - (mL + np.cross(r[i],nL))
        
        residual[36:39] = EF    #residual for the external force balance
        residual[39:42] = EM    #residual for the external moment balance
        # print(f"residual vector: {residual}")
        # print(f"residual magnitude: {np.linalg.norm(residual)}")
        return residual #residual vector for the system
    
    def rod_ode(s, y):
        '''System of ODEs for the states of the rod'''  
        del s  # Integration variable unused in autonomous ODE
        
        #Unpack the state vector
        R = np.reshape(y[3:12], (3, 3)) #orientation at the base of the flexible rod
        n = y[12:15]                    #internal force vector
        m = y[15:18]                    #internal moment vector

        '''Linear Constitutive relation between the deformation and the internal forces'''
        v = KseInv @ (R.T @ n) #linear strain in local frame at the base of the flexible rod
        # print(f"v: {v}")
        v[2] = v[2] + 1        #For a straight rod
        u = KbtInv @ (R.T @ m) + np.array([1/0.3,0,0]) #angular strain in local frame-->for a curved rod about x-axis

        y_out = np.empty(y.shape)
        
        #ODEs for the rod states: 
        y_out[0:3] = R @ v                      #change in the position p0
        # np.matmul(R, v, y_out[0:3]) 
        y_out[3:12] = hat_postmultiply(R, u)    #change in the orientation R0    
        y_out[12:15] = -rho * A * g             #internal force per unit length of the rod
        y_out[15:18] = -np.cross(y_out[0:3], n) #internal moment per unit length of the rod

        # Pack state vector derivative
        return y_out

    def hat_postmultiply(m, v):
        '''multiplication operation of rotation matrix with skew symmetric matrix
        ouput = (9,) vector '''

        return np.array([m[0,1]*v[2] - m[0,2]*v[1], m[0,2]*v[0] - m[0,0]*v[2], m[0,0]*v[1] - m[0,1]*v[0],
        m[1,1]*v[2] - m[1,2]*v[1], m[1,2]*v[0] - m[1,0]*v[2], m[1,0]*v[1] - m[1,1]*v[0],
        m[2,1]*v[2] - m[2,2]*v[1], m[2,2]*v[0] - m[2,0]*v[2], m[2,0]*v[1] - m[2,1]*v[0]])
        
####################################################################################################################################### 
    '''Non-linear optimization loop to update the unknown variables'''
    Y = None
    options = {'ftol': 1e-5, 'xtol': 1e-5} #tolerance for the solver solution
    # Record the start time
    start_time = time.time()

    ''' levenberg marquardt algorithm 'lm' '''
    sol = scipy.optimize.root(residual_function, init_guess, method='lm', jac = False, tol=1e-5, options=options).x
    # print(f"sol: {sol}")
    # Record the end time of code
    end_time = time.time()
    # Calculate the total time taken
    total_time = end_time - start_time
    # Print the total time taken
    # print(f"sol:  {sol}")
    print(f"Total time taken: {total_time} seconds")

       
    #return variables for analysis
    return sol, total_time



def save_excel(path, valid_array):
    '''function to save into excel file (mainly the motor angles as an input to the FK model)'''
    with xlsxwriter.Workbook(path) as workbook:
        worksheet = workbook.add_worksheet()

        red_format = workbook.add_format({'bg_color': '#FF0000', 'font_color': '#FF0000'})
        green_format = workbook.add_format({'bg_color': '#00FF00', 'font_color': '#00FF00'})
        
        row_color = []  # To keep track of row colors

        for row_num, data in enumerate(valid_array):
            worksheet.write_row(row_num, 0, data)
            last_values = data[-14:]


def helicalTraj():
    '''function to sample 3D points from a helical trajectory'''
    # Parameters
    radius = 0.25  # Radius of the helix
    pitch = 1.0  # Pitch of the helix
    z_start = 0.4  # Starting z height
    z_end = 0.55  # Ending z height

    # Create time values
    t = np.linspace(0, 2 * np.pi, 500)

    # Parametric equations for the helix
    x = radius * np.cos(t)
    y = radius * np.sin(t)
    z = np.linspace(z_start, z_end, len(t))

    # Plot the helix in 3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z, label='Helical Trajectory', linewidth=2)

    # Set axis labels
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_zlabel('Z-axis')

    # Show the plot
    plt.legend()
    # plt.show()
    return np.column_stack((x, y, z)) #samples (500,3)


saveData_bool=True #saving data to an excel file 
IK_vec = [] #saves the results for plotting#Initialization for the inverse kinematics

#actuator and universal joint angles initialization
qi = np.array([0,0,0,
               0,0,0,
               0,0,0,
               0,0,0,
               0,0,0,
               0,0,0])

#initializing the guess vector for the IK model
init_guess = np.concatenate([np.zeros(24),qi]) #42 variables
#orientation of the end-effector platform about x,y,z axis in world coordinate
angles = np.array([np.deg2rad(0),np.deg2rad(0),np.deg2rad(0)])

#extracting 10 samples from the helical trajectory 
rand_sample_index = np.random.choice(helicalTraj().shape[0], size=10, replace=False)
helical_samples = helicalTraj()
#end-effector position as an input for the IK model
p_ee = np.array(helical_samples[rand_sample_index])
i =0

##############################################################################################
'''compute joint angles for each end-effector position '''
while i<len(p_ee):
    print(i) #print the sample number

    #excuting the IK function
    optimised_states, total_time = Inverse_Kinetostatic_traj(p_ee[i], angles, init_guess)
    # print(f"optimised_states: {restrack.shape}")
    IK_vec.append(np.concatenate(([total_time, optimised_states[24]],[optimised_states[27]],[optimised_states[30]],[optimised_states[33]],[optimised_states[36]],[optimised_states[39]], p_ee[i])))
    # print(f"linAct_error: {linAct_error}")
    # init_guess = optimised_states
    i+=1

save_excel('path-to-IK_PCR_Trajectory.xlsx',
                   np.array(IK_vec))
