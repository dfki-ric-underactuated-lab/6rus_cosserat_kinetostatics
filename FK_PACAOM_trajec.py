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


def main(qm, init_guess):

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
    
    
    #initialisation of the geometry for the rigid part
    t = 0.29711 		#length of edge of the end-effector triangle 
    b = 0.21841     	#length of edge of the base triangle       
    d = 0.05374     	#distance between the spherical connection points at the end-effector
    D = 0.0675     	#distance between the motors center point at the base
    l1 = 0.23164		#length of the link connecting motor and universal joint
    L = 0.53		#curved length of the flexible link at a stress free height of 0.4m
    Dist = np.array([D,0,0])
    dist = np.array([[d],[0],[0]])
    Rz_120 = rotation_matrix_3d(np.deg2rad(0),np.deg2rad(0),np.deg2rad(120))
  
    #base triangle vertices in world coordinates
    B = np.array([0,(2/3)*b*np.cos(30*np.pi/180),0])

	#motor position in world coordinates using the base traingle vertices
	#where M1, M2,...,M6 are the motor position in world coordinates
    M1 = B + Dist
    M2 = B - Dist
    M3 = Rz_120 @ M1 #B2 + Rz120.dot(Dist)
    M4 = Rz_120 @ M2#B2 - Rz120.dot(Dist)
    M5 = Rz_120.T @ M1
    M6 = Rz_120.T @ M2
    #orientation of the motor joint wrt world
    Mm = np.array([M1,M2,M3,M4,M5,M6])
    #rotation of motor coo to world coo
    Rm = np.array([np.eye(3),np.eye(3),Rz_120,Rz_120,Rz_120.T,Rz_120.T])
    #End-effector triangle vertices in the EE coordinate
   
    #Transformation matrix--> EE to base frame 
    
    def findEElinkpoints(p_ee, Ree):
        # '''position of the rod connections at the EE'''
        # theta_L = (-alpha1 + ((i+1) - (i+1)%2)*alpha1 + (i-i%2)*alpha2) / 2
        # #r is the position of the each attachment at EE in the world frame
        # r = scrib_R * np.array([np.cos(theta_L), np.sin(theta_L), 0])
        # r = p_ee + (R_ee @ r)
        # return r
        T = np.array([[0],[-(2/3)*t*np.cos(30*np.pi/180)],[0]]) #length of vector from EE to the moving platform edge
    
        #EE position in world coordinates using the EE traingle vertices
        #where C1, C2,...,C6 are the spherical joint position in EE coordinates
        C5 = T + dist
        C4 = T - dist
        C6 = Rz_120 @ C4	#T2 - Rz120 .dot(dist)
        C1 = Rz_120 @ C5	#T2 + Rz120 .dot(dist)
        C2 = Rz_120.T @ C4	#T3 - Rz120_ .dot(dist)
        C3 = Rz_120.T @ C5	#T3 + Rz120_.dot(dist)
        O_Tee_EE = np.hstack((np.vstack((Ree, np.zeros((1,3)))), np.vstack(((np.reshape(p_ee,(3,1))), 1))))
        World_position_SphericalJoints =  O_Tee_EE @ (np.hstack((np.vstack((C1, 1)),
                                        np.vstack((C2, 1)),np.vstack((C3, 1)),np.vstack((C4, 1)),
                                        np.vstack((C5, 1)),np.vstack((C6, 1))))) #each column of the matrix is the new spherical joint 
        #orientation of the spherical joint wrt world
        RC = [Ree @ Rz_120, Ree @ Rz_120.T, Ree @ Rz_120.T, Ree @ np.eye(3), Ree @ np.eye(3), Ree @ Rz_120]
        #extracting each column from the above matrix which contains transformed spherical joint positions
        C1 = World_position_SphericalJoints[:,0][:-1]
        C2 = World_position_SphericalJoints[:,1][:-1]
        C3 = World_position_SphericalJoints[:,2][:-1]
        C4 = World_position_SphericalJoints[:,3][:-1]
        C5 = World_position_SphericalJoints[:,4][:-1]
        C6 = World_position_SphericalJoints[:,5][:-1]

        return np.array([C1, C2, C3, C4, C5, C6])
    
    p0 = np.empty((6,3))
    
    for i in range(6):
        N_local = np.array([0, l1*np.cos(qm[i]), l1*np.sin(qm[i])])
        #universal in world coo
        p0[i] = Mm[i] + (Rm[i] @ N_local)
    
    # print(f"shape of p0 is {p0[0].shape}")

    # Independent Parametersfor the flexible link
    E = 110.3e9  # Young's modulus (N/m2)
    nu = 0.31   #tranversal elongation/axial compression
    G = E/(2*(1+nu))  # Shear modulus (N/m2)
    rho = 4428.8  # Density (kg/m3)
    rad = 0.002
    h = 0.0012 #thickness of the rectangular cross-section (m)
    w = 0.077  #width of the cross-section (m)
    ee_mass = 0.5#0.000001 #Kg #neglible mass 0.000001
    g = np.array([0,0,-9.81]) 

    # Dependent Parameters for the flexible link
    # A = h * w  #Cross-sectional area  # Cross-sectional area
    # Ixx = (w * h**3)/12  #Area moment of inertia
    # Iyy = (h * w**3)/12  # Area moment of inertia
    # J = Ixx + Iyy
    A = math.pi * rad ** 2  # Cross-sectional area
    I = math.pi * rad ** 4 / 4  # Area moment of inertia
    J = 2 * I  # Polar moment of inertia
    Kse = np.diag([G * A, G * A, E * A])  # shear and extension
    KseInv = np.linalg.inv(Kse)
    Kbt = np.diag([E * I, E * I, G * J]) #bending and torsion
    KbtInv = np.linalg.inv(Kbt)
    # print(f"Kse: {Kse}")
    # print(f"Kbt: {Kbt}")

    F = np.array([0, 0, -9.81*ee_mass]) #+ np.array([10,0,0])
    M = np.zeros(3)
    restrack = [] #tracks the residual magnitude for each initial update

 
    #Sub-functions
    def obj_fun(guess):
        '''objective function aka cost function'''
        # Update guessed initial conditions
        residual = np.empty(42) #42 constraint equations
        res_mag = np.empty(14) #store the residual vector magnitude
        EF = F
        EM = M
        for i in range(6):
            n0 = guess[6*i:6*i+3]
            m0 = guess[6*i+3]
            q = guess[6*i+4:6*i+6] #only universal joint angles (FK)
            p_ee = guess[36:39]
            R_ee = guess[39:42]
            R_ee = RPY_to_matrix(R_ee[0], R_ee[1], R_ee[2])#mapping RPY to matrix
            r = findEElinkpoints(p_ee, R_ee)
                      
            Rx = rotation_matrix_3d(q[1], 0, 0) #theta_x, theta_y, theta_z
            Ry = rotation_matrix_3d(0, q[0], 0) #theta_x, theta_y, theta_z
            #rotation matrix for the universal joints about y and x axis given by Ryx = Ry @ Rx
            R0 = Rm[i] @ Ry @ Rx
            # circradius = R0 @ np.array([0,0.2222,0])
            # print("circular radius",np.linalg.norm(circradius))
            y0 = np.concatenate([p0[i], np.reshape(R0, 9), n0,[0,0], [m0]])

            # Numerically solve the IVP
            nonlocal Y
            num_steps = 100  # Total number of sampled integration steps
            t_eval = np.linspace(0, L, num_steps)
            #‘RK45’ (default): Explicit Runge-Kutta method of order 5(4)
            Y = integrate.solve_ivp(rod_ode, (0.0, L), y0, t_eval = t_eval, atol=1e-5, rtol=1e-5).y
        
            # Calculate distal constraint violation
            pL_shot = Y[0:3, -1]
            nL = Y[12:15, -1]
            mL = Y[15:18, -1]
            # print(f"pL_shot: {R_ee.shape}")
            residual[6*i:6*i+3] = pL_shot - r[i]
            res_mag[2*i] = np.linalg.norm(residual[6*i:6*i+3])
            residual[6*i+3:6*i+6] = R_ee.T @ mL #spherical joint
            res_mag[2*i+1] = np.linalg.norm(residual[6*i+3:6*i+6])
            
            #F and M respectively = 0
            EF = EF - nL
            EM = EM - (mL + np.cross(r[i],nL))
        
        residual[36:39] = EF
        res_mag[12] = np.linalg.norm(residual[36:39])
        residual[39:42] = EM
        res_mag[13] = np.linalg.norm(residual[39:42])
        # print(f"residual vector: {residual}")
        # print(f"residual magnitude: {np.linalg.norm(residual)}")
        restrack.append(res_mag)
        return residual #Error vector
    
    def rod_ode(s, y):
        '''State vector derivative function'''  
        del s  # Integration variable unused in autonomous ODE
        # Unpack state vector
        R = np.reshape(y[3:12], (3, 3))
        n = y[12:15]
        m = y[15:18]

        # Constitutive equation
        # Constitutive equation
        v = KseInv @ (R.T @ n) #+ np.array([0,-np.cos(np.deg2rad(15)), np.sin(np.deg2rad(15))])
        # print(f"v: {v}")
        v[2] = v[2] + 1 #radius of p(0) initial configuration
        u = KbtInv @ (R.T @ m) + np.array([1/0.3005,0,0]) #0.2222

        y_out = np.empty(y.shape)
        
        y_out[0:3] = R @ v
        # np.matmul(R, v, y_out[0:3]) 
        y_out[3:12] = hat_postmultiply(R, u)
        y_out[12:15] = -rho * A * g
        y_out[15:18] = -np.cross(y_out[0:3], n)

        # Pack state vector derivative
        return y_out

    def hat_postmultiply(m, v):
        '''multiplication operation of rotation matrix with skew symmetric matrix
        ouput = (9,) vector '''

        return np.array([m[0,1]*v[2] - m[0,2]*v[1], m[0,2]*v[0] - m[0,0]*v[2], m[0,0]*v[1] - m[0,1]*v[0],
        m[1,1]*v[2] - m[1,2]*v[1], m[1,2]*v[0] - m[1,0]*v[2], m[1,0]*v[1] - m[1,1]*v[0],
        m[2,1]*v[2] - m[2,2]*v[1], m[2,2]*v[0] - m[2,0]*v[2], m[2,0]*v[1] - m[2,1]*v[0]])

        
    
    # init_guess = np.concatenate([np.zeros(30),p_ee[2]*np.ones(6)])
    Y = None
    options = {'ftol': 1e-5, 'xtol': 1e-5}
    # Record the start time
    start_time = time.time()
    sol = scipy.optimize.root(obj_fun, init_guess, method='lm', jac = False, tol=1e-5, options=options).x
    # print(f"sol: {sol}")
    # sol = 0
    # Record the end time of code
    end_time = time.time()
    # Calculate the total time taken
    total_time = end_time - start_time
    # Print the total time taken
    # print(f"sol:  {sol}")
    print(f"Total time taken: {total_time} seconds")

    '''sol.x = contains the optimised guess vector which minimizes the obj_fun
    ptrack = stores the states of each leg for every iteration of the optimization
    '''
    
    return sol, total_time, restrack[-1],qm, F




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
    return None



plot_bool=False #plotting function
saveData_bool=True #saving data to an excel file 
NormalPlot = True
InteractivePlot = True
ls_value = []
p_ee_FK = []#pee calculated from FK

# Read the Excel file into a DataFrame
df = pd.read_excel('/home/dfki.uni-bremen.de/vrodrigues/6_rus_Stewart_platform/Simulation/Flexible_links/Cosserat_rod_theory/excel_files/Trajectories/PACOMA_IK_trajec_helical_5N_ROD.xlsx', header=None)
df_m = df.iloc[:, 1:7]#motor angles
# print(df_m.head())
# df_pee = df.iloc[:, 6:9]#pee

# Convert each row into a NumPy array
motor_angle = [np.array(row) for _, row in df_m.iterrows()]
# p_ee = [np.array(row) for _, row in df_pee.iterrows()]

# p_ee = np.array(p_ee)#pee for FK
motor_angle = np.array(motor_angle) #motor angles for FK
p_ee = np.array([0.25,0,0.4])#pee for FK
R_ee = np.array([np.deg2rad(0), np.deg2rad(0), np.deg2rad(0)])

#joint angle initialization
qi = np.array([0,0,
               0,0,
               0,0,
               0,0,
               0,0,
               0,0])

init_guess = np.concatenate([np.zeros(24),qi,p_ee,R_ee]) #42 variables

i =0
while i < len(motor_angle):
    # L = np.array([0.43218561, 6.1670081 , 2.17345386, 0.43806599, 1.40285005 ,0.61917025]) #known actuator vaues Li = qi + li
    print(i)
    # print(motor_angle[i])
    optimised_states, total_time, restrack, motorangle, F = main(motor_angle[i], init_guess)
    p_ee_FK.append(np.concatenate(([total_time], F, optimised_states[36:39],motorangle, restrack.flatten())))
    # init_guess = optimised_states
    i+=1

p_ee_FK = np.array(p_ee_FK)
save_excel('/home/dfki.uni-bremen.de/vrodrigues/6_rus_Stewart_platform/Simulation/Flexible_links/Cosserat_rod_theory/excel_files/Trajectories/PACOMA_FK_trajec_helical_5N_ROD_winit.xlsx',p_ee_FK)