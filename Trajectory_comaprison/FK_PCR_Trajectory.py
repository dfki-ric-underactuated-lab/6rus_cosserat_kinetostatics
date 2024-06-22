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


def Forward_Kinetostatic_traj(qm, init_guess):

    ''' '''
    
    #calculate rotation matrices @ x,y,z
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
  
    #
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
   
    
    def findEElinkpoints(p_ee, Ree):

        '''for a given pose of the end-effector pose it computes the 
        spherical joint position and orientation in world frame'''

        T = np.array([[0],[-(2/3)*t*np.cos(30*np.pi/180)],[0]]) 
        #Spherical joints positions C1, C2,...,C6 in EE coordinates
        C5 = T + dist
        C4 = T - dist
        C6 = Rz_120 @ C4	#T2 - Rz120 .dot(dist)
        C1 = Rz_120 @ C5	#T2 + Rz120 .dot(dist)
        C2 = Rz_120.T @ C4	#T3 - Rz120_ .dot(dist)
        C3 = Rz_120.T @ C5	#T3 + Rz120_.dot(dist)

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
        return np.array([C1, C2, C3, C4, C5, C6])
    
    
    p0 = np.empty((6,3))
    for i in range(6):

        '''compute the position of the base of the flexible link from 
            the motor angles qm'''
        
        N_local = np.array([0, l1*np.cos(qm[i]), l1*np.sin(qm[i])])
        #universal in world coo
        p0[i] = Mm[i] + (Rm[i] @ N_local)
    

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


    def residual_function(guess):
        '''objective function aka cost function'''
        # Update guessed initial conditions
        residual = np.empty(42) #42 constraint equations
        EF = F
        EM = M
        for i in range(6):
            n0 = guess[6*i:6*i+3]   #n_x(0), n_y(0), n_z(0)
            m0 = guess[6*i+3]       #m_z(0)
            q = guess[6*i+4:6*i+6]  #universal joint angles
            p_ee = guess[36:39]     #position of the center of the end-effector
            R_ee = guess[39:42]     #orientation of the center of the end-effector
            R_ee = RPY_to_matrix(R_ee[0], R_ee[1], R_ee[2])#mapping RPY to matrix

            #compute the spherical joint position for the pose of the end-effector
            r = findEElinkpoints(p_ee, R_ee)

            #orientation at the base of the flexible rod in world coordinate frame
            #orientation of the universal joint in local joint frame    
            Rx = rotation_matrix_3d(q[1], 0, 0) #theta_x, theta_y, theta_z
            Ry = rotation_matrix_3d(0, q[0], 0) #theta_x, theta_y, theta_z

            #orientation at the base of the flexible rod in world coordinate frame 
            R0 = Rm[i] @ Ry @ Rx
            
            '''Solve the IVP problem for each rod'''
            #concatenate the flexible rod states at the base
            y0 = np.concatenate([p0[i], np.reshape(R0, 9), n0,[0,0], [m0]])

            # Numerically solve the IVP
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
            residual[6*i+3:6*i+6] = R_ee.T @ mL      #residual for the moment error for spherical joint
            
            #F and M respectively = 0
            EF = EF - nL
            EM = EM - (mL + np.cross(r[i],nL))
        
        residual[36:39] = EF    #residual for the external force balance
        residual[39:42] = EM    #residual for the external moment balance
        # print(f"residual vector: {residual}")
        # print(f"residual magnitude: {np.linalg.norm(residual)}")
        return residual #Error vector
    
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

    '''sol.x = contains the optimised guess vector which minimizes the residual_function
    ptrack = stores the states of each leg for every iteration of the optimization
    '''
    
    return sol, total_time

######################################################################################################################


def save_excel(path, valid_array):
    '''
    function to save into excel file:
    valid_array = type/shape of the solution array should be valid for the below excel library
    path  = '/path/of/the/excelfile.xlsx'
    
    '''
    #save the date of the states in an excel file
    with xlsxwriter.Workbook(path) as workbook:
        worksheet = workbook.add_worksheet()

        for row_num, data in enumerate(valid_array):
            worksheet.write_row(row_num, 0, data)
    return None

def helicalTraj():
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
    return x, y, z
def FK_plots():

    file_paths = ['/home/dfki.uni-bremen.de/vrodrigues/6_rus_Stewart_platform/Simulation/Flexible_links/Cosserat_rod_theory/excel_files/Trajectories/PACOMA_IK_trajec_helical_5N_ROD_test.xlsx',
        '/home/dfki.uni-bremen.de/vrodrigues/6_rus_Stewart_platform/Simulation/Flexible_links/Cosserat_rod_theory/excel_files/Trajectories/PACOMA_FK_trajec_helical_5N_ROD_test.xlsx']
 
    
    dfIK5C  =pd.read_excel(file_paths[0], header=None)
    dfFK5C  =pd.read_excel(file_paths[1], header=None)
    
    fontsize = 36
    
    distances = np.linalg.norm(dfIK5C[dfIK5C.columns[7:10]].values - dfFK5C[dfFK5C.columns[1:4]].values, axis=1)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    x,y,z = helicalTraj()
    ax.set_xlabel('x (m)', fontdict = {'size': fontsize},labelpad=10)
    ax.set_ylabel('y (m)', fontdict = {'size': fontsize},labelpad=20)
    ax.set_zlabel('z (m)', fontdict = {'size': fontsize}, labelpad=20)
    ax.scatter(x,y,z,label='Reference trajectory',c='g',linewidth=1.0)
    # ax.scatter(dfIK5C[7], dfIK5C[8], dfIK5C[9], label='Reference EE position',marker="o", c='r', alpha=0.5)
    ax.scatter(dfFK5C[1], dfFK5C[2], dfFK5C[3], label='FK solution', s = 80,facecolors='none', edgecolors='k')


    ax.tick_params(axis='both', which='both', labelsize=fontsize-2)

    # Set the number of ticks or specify a locator
    ax.xaxis.set_major_locator(plt.MaxNLocator(6))  # Adjust the number (6 in this example)
    ax.yaxis.set_major_locator(plt.MaxNLocator(6))  # Adjust the number (6 in this example)
    ax.zaxis.set_major_locator(plt.MaxNLocator(6))  # Adjust the number (6 in this example)
    ax.set_zticks(np.arange(0.3, 0.6, 0.05))

    ax.legend(loc='best',fontsize=fontsize-8)
    
    # Plot the Euclidean distance error
    fig, ax_error = plt.subplots()

    ax_error.xaxis.set_major_locator(plt.MaxNLocator(10))
    ax_error.yaxis.set_major_locator(plt.MaxNLocator(10))

    #Plot the distances
    ax_error.scatter( dfIK5C.index, distances, marker='o',s =40)

    ax_error.set_xlabel('Point Index', fontdict={'size': fontsize})
    ax_error.set_ylabel('Distance Error (m)', fontdict={'size': fontsize})
    # from matplotlib.ticker import FormatStrFormatter
    plt.tick_params(axis='both', which='both', labelsize=fontsize-2)
    
    # # Move the legend to the upper left corner
    ax_error.legend(loc='upper left')

    plt.show()



saveData_bool=True #saving data to an excel file 

ls_value = []
FK_vec = []#pee calculated from FK



#intial guess for the pose of the end-effector
p_ee = np.array([0.25,0,0.4])#pee for FK
R_ee = np.array([np.deg2rad(0), np.deg2rad(0), np.deg2rad(0)])

#universal joint angle initialization
qi = np.array([0,0,
               0,0,
               0,0,
               0,0,
               0,0,
               0,0])

'''For trajectory comparison'''
# Read the joint angles from the saved IK Excel file
df = pd.read_excel('path-to-IK_PCR_Trajectory.xlsx', header=None)
df_m = df.iloc[:, 1:7]  #unpack the motor angles
# Convert each row into a NumPy array
motor_angle = [np.array(row) for _, row in df_m.iterrows()]
motor_angle = np.array(motor_angle) #motor angles for FK

#initializing the guess vector for the FK model
init_guess = np.concatenate([np.zeros(24),qi,p_ee,R_ee]) #42 variables
i =0

##############################################################################################
'''compute pose of the end-effector for given motor angles '''
while i < len(motor_angle):
    # L = np.array([0.43218561, 6.1670081 , 2.17345386, 0.43806599, 1.40285005 ,0.61917025]) #known actuator vaues Li = qi + li
    print(i)
    # print(motor_angle[i])
    optimised_states, total_time = Forward_Kinetostatic_traj(motor_angle[i], init_guess)
    FK_vec.append(np.concatenate(([total_time], optimised_states[36:39])))
    # init_guess = optimised_states
    i+=1

FK_vec = np.array(FK_vec)


save_excel('path-to-FK_PCR_Trajectory.xlsx',FK_vec)

FK_plots()