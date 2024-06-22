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


def Workspace(p_ee, R_ee, init_guess):

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

    Ree = RPY_to_matrix(R_ee[0], R_ee[1], R_ee[2]) #orientation of the end-effector plate
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

    restrack = [] #tracks the residual magnitude for each initial update
    counter = 0     #counter variable for the optimisation loop

    def residual_function(guess):
        '''Residual vector function which includes the constraints that need to be zero for a solution'''
        # Update guessed initial conditions
        nonlocal counter #count the function calls
        counter += 1
        # print(counter)
        residual = np.empty(42) #42 constraint equations
        res_mag = np.empty(14) #store the residual vector magnitude
        
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
            res_mag[2*i] = np.linalg.norm(residual[6*i:6*i+3])
            residual[6*i+3:6*i+6] = Ree.T @ mL      #residual for the moment error for spherical joint
            res_mag[2*i+1] = np.linalg.norm(residual[6*i+3:6*i+6])
            res_mag[13] = np.linalg.norm(residual[39:42])
            
            #F and M respectively = 0
            EF = EF - nL
            EM = EM - (mL + np.cross(r[i],nL))
        
        residual[36:39] = EF    #residual for the external force balance
        res_mag[12] = np.linalg.norm(residual[36:39])
        residual[39:42] = EM    #residual for the external moment balance
        # print(f"residual vector: {residual}")
        # print(f"residual magnitude: {np.linalg.norm(residual)}")
        restrack.append(res_mag)
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
    return sol, restrack[-1], total_time


def save_excel(path, valid_array):
    '''
    valid_array = type/shape of the solution array should be valid for the below excel library
    path  = '/path/of/the/excelfile.xlsx'
    
    '''
    #save the date of the states in an excel file
    with xlsxwriter.Workbook(path) as workbook:
        worksheet = workbook.add_worksheet()

        for row_num, data in enumerate(valid_array):
            worksheet.write_row(row_num, 0, data)

def workspace_analysis(path):

    df = pd.read_excel(path, header=None)

    merged_condition = (
    (df[df.columns[22]] <= 5e-10) &     #position errors
    (df[df.columns[23]] <= 5e-10) &     #        .         
    (df[df.columns[24]] <= 5e-10) &     #        . 
    (df[df.columns[25]] <= 5e-10) &     #        . 
    (df[df.columns[26]] <= 5e-10) &     #        . 
    (df[df.columns[27]] <= 5e-10) &     #        .
    (np.rad2deg(df[df.columns[1]]) <=90)& (np.rad2deg(df[df.columns[1]]) > -20)&    #motor angle limitations
    (np.rad2deg(df[df.columns[2]]) <=90)& (np.rad2deg(df[df.columns[2]]) > -20)&    #       .            
    (np.rad2deg(df[df.columns[3]]) <=90)& (np.rad2deg(df[df.columns[3]]) > -20)&    #       .            
    (np.rad2deg(df[df.columns[4]]) <=90)& (np.rad2deg(df[df.columns[4]]) > -20)&    #       .            
    (np.rad2deg(df[df.columns[5]]) <=90)& (np.rad2deg(df[df.columns[5]]) > -20)&    #       .             
    (np.rad2deg(df[df.columns[6]]) <=90)& (np.rad2deg(df[df.columns[6]]) > -20)&    #       .            
    (df[df.columns[28]] <= 5e-10) &     #internal moment limitations
    (df[df.columns[29]] <= 5e-10) &     #           .
    (df[df.columns[30]] <= 5e-10) &     #           .
    (df[df.columns[31]] <= 5e-10) &     #           .
    (df[df.columns[32]] <= 5e-10) &     #           .
    (df[df.columns[33]] <= 5e-10) &     #           .
    (df[df.columns[34]] <=5e-10)& #external force limitations
    (df[df.columns[35]] <=5e-10)) #external moment limitations
    
    df[36] = merged_condition.astype(int)
    print(len(df))
    df = df[merged_condition].copy()
    print(len(df))

    fontsize = 32

    # Scatter plot for Motor 1
    plt.scatter(df.index, np.rad2deg(df[df.columns[1]]), label='M1', s=3, marker='o')
    plt.scatter(df.index, np.rad2deg(df[df.columns[2]]), label='M2', s=3, marker='s')
    plt.scatter(df.index, np.rad2deg(df[df.columns[3]]), label='M3', s=3, marker='^')
    plt.scatter(df.index, np.rad2deg(df[df.columns[4]]), label='M4', s=3, marker='D')
    plt.scatter(df.index, np.rad2deg(df[df.columns[5]]), label='M5', s=3, marker='p')
    plt.scatter(df.index, np.rad2deg(df[df.columns[6]]), label='M6', s=3, marker='h')

    
    plt.ylabel('Motor angles (deg)', fontdict = {'size': fontsize})
    plt.xlabel('samples at each 2d cross section along the z-height', fontdict = {'size': fontsize})
    plt.legend(loc='best',fontsize=fontsize-6)    # axs[i].legend()


    plt.tick_params(axis='both', which='both', labelsize=fontsize-6)

    plt.xticks(np.arange(0, 4000, step=500))
    plt.yticks(np.arange(0, 89, step=5))
    # plt.yticks(np.concatenate([np.arange(0, 100, step=25), np.arange(100, 800, step=50)]))


    #################################################################################################
    #plotting the filtered EE pose
    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Scatter plot for points where merged_condition is True (1)
    ax.scatter(df.loc[df[36] == 1, 19], df.loc[df[36] == 1, 20], df.loc[df[36] == 1, 21], label='Reachable workspace', marker='o', c='darkgreen', alpha=1)
    
    plt.tick_params(axis='both', which='both', labelsize=fontsize-2)

    plt.legend(loc='best',fontsize=fontsize-6)

    # Show the plot
    plt.show()



plot_bool=False #plotting function
saveData_bool=True #saving data to an excel file 
NormalPlot = True
InteractivePlot = True
ls_value = []
#Initialization for the inverse kinematics
# p_ee = np.array([0,0,0.4])
#joint angle initialization
qi = np.array([0,0,0,
               0,0,0,
               0,0,0,
               0,0,0,
               0,0,0,
               0,0,0])

init_guess = np.concatenate([np.zeros(24),qi]) #42 variables
# print(f"len of guess vec: {init_guess.shape}")
R_ee = np.array([np.deg2rad(0),np.deg2rad(0),np.deg2rad(0)])
#excuting the main function
linAct_error = [] #saves the error magnitude and joint angle to evaluate a valid solution

def generate_random_points_inside_circle(radius, num_points, height):
    # Generate random points inside a circle
    theta = np.random.uniform(0, 2 * np.pi, num_points)
    r = np.sqrt(np.random.uniform(0, radius**2, num_points))
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    z = np.full_like(x, fill_value=height)
    return x, y, z

def generate_points_arrays(radius_vector, height_vector, num_points=250):
    all_points = []

    for radius, height in zip(radius_vector, height_vector):
        x, y, z = generate_random_points_inside_circle(radius, num_points, height)
        points = np.column_stack((x, y, z))
        all_points.append(points)

    all_points = np.concatenate(all_points)

    return all_points

# Example usage
circle_radius = np.array([0.4,0.4,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.3,0.3,0.2,0.2,0.2,0.1,0.05])

height_vector = np.arange(0.4, 0.72, 0.02)

points = generate_points_arrays(circle_radius, height_vector, num_points=250)

p_ee = points
# print(p_ee.shape)

# print(len(p_ee))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(p_ee[:,0],p_ee[:,1],p_ee[:,2],label='Reference trajectory')

q1_vec = np.empty(6)
q2_vec = np.empty(6)
q3_vec = np.empty(6)
# p_ee = np.array([0, 0, 0.4])
i =0
while i<5:
    print(i)
    print(p_ee[i])

    #excuting the main function
    optimised_states, restrack, total_time = Workspace(p_ee[i], R_ee, init_guess)
    # print(f"optimised_states: {restrack.shape}")
    for j in range(6):
        q1_vec[j] = optimised_states[3*j+24]
        q2_vec[j] = optimised_states[3*j+25]
        q3_vec[j] = optimised_states[3*j+26]
    linAct_error.append(np.concatenate(([total_time],q1_vec, q2_vec, q3_vec, p_ee[i],restrack.flatten())))
    # print(f"linAct_error: {linAct_error}")
    # init_guess = optimised_states
    i+=1

save_excel('/home/dfki.uni-bremen.de/vrodrigues/6_rus_Stewart_platform/Simulation/Flexible_links/Cosserat_rod_theory/excel_files/PACOMA_Rotation_EE_workspace_ROD_test.xlsx',
                   np.array(linAct_error))

workspace_analysis('/home/dfki.uni-bremen.de/vrodrigues/6_rus_Stewart_platform/Simulation/Flexible_links/Cosserat_rod_theory/excel_files/PACOMA_Rotation_EE_workspace_ROD.xlsx')
