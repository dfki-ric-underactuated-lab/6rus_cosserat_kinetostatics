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


def Forward_Kinetostatic(init_guess, qm):

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
    ptrack = [] #tracks the position of the flexible link in space
    nLtrack = [] #track the force acting at the last rod segment
    pLtrack = [] #track the position of the last cross-section of the rod 
    rtrack =[]   #track the position of the spherical joint position

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
            rtrack.append(r)

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
            ptrack.append(Y[0:3,:])
        
             #residual constraints at the distal end including the force balance
            '''extract the solution of the IVP for each rod'''
            pL_shot = Y[0:3, -1]        #extract the position vector of the tip of the rod
            pLtrack.append(Y[0:3, -1]) 
            nL = Y[12:15, -1]           #extract the force vector at the tip of the rod
            nLtrack.append(Y[12:15, -1]) 
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
    ptrack = np.array(ptrack)
    nLtrack = np.array(nLtrack)
    ptrack = ptrack.reshape(ptrack.shape[0] * ptrack.shape[1], ptrack.shape[2])

    return sol, np.array(rtrack), Mm, ptrack, nLtrack,F, pLtrack[-6:],p0

######################################################################################################################
def plotting(path_excel, path_img_interactive, path_img):

    '''
    path_excel1 = '/path/of/the/saved/excel1.xlsx'
    path_img = '/path/of/the/output.html'
    '''
    all_sheets = pd.read_excel(path_excel, sheet_name=None)
    # df1 = pd.read_excel(path_excel1, header=None)
    # df2 = pd.read_excel(path_excel2, header=None)
    df1 = all_sheets['Sheet1']
    df2 = all_sheets['Sheet2']
    nL_norm_sum = np.zeros(3)
    nL_norm = all_sheets['Sheet3'].tail(6).values
    for i in range(len(nL_norm)):
        nL_norm_sum = nL_norm_sum + nL_norm[i]
    #printing the normal forces at the last cross-section of the rod.
    print(f"nL magnitude: {np.linalg.norm(nL_norm_sum)} and the individual forces norm: {np.linalg.norm(nL_norm[0])},{np.linalg.norm(nL_norm[1])},{np.linalg.norm(nL_norm[2])},{np.linalg.norm(nL_norm[3])},{np.linalg.norm(nL_norm[4])},{np.linalg.norm(nL_norm[5])}")


    r = df2.tail(6).values
    last_18_rows = df1.tail(18)
    N = 100

    # Reshape the DataFrame to get 6 sets of (3, N) coordinates
    coordinates = last_18_rows.values.reshape(6, 3, N)

    if NormalPlot:
         # # Create a 3D plot
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        label_font = {'family': 'sans-serif', 'weight': 'normal', 'size': 14}


        # # Plot each set of (3, 3) coordinates
        labels = ['leg 1','leg 2', 'leg 3', 'leg 4', 'leg 5', 'leg 6']
        for coord_set, label in zip(coordinates, labels):
            ax.scatter(coord_set.T[:, 0], coord_set.T[:, 1], coord_set.T[:, 2], s = 10, label= label)

        ax.scatter(p0[:, 0], p0[:, 1], p0[:, 2],c='black', marker='o', label='Universal joints')
        for i, txt in enumerate(range(1, len(p0) + 1)):
            # x,y,z = [p0[i][0]],[p0[i][0]],[p0[i][0]]

            x,y,z = [r[i-1][0], r[i][0]], [r[i-1][1], r[i][1]], [r[i-1][2], r[i][2]]
            ax.plot(x,y,z, color='black',linewidth=5.0)
            ax.text(r[i][0], r[i][1], r[i][2], str(txt), ha='right', va='bottom')
            ax.text(p0[i][0], p0[i][1], p0[i][2], str(txt), ha='right', va='bottom')

        # ax.legend(p0)
        #plot motor positions
        # Plot points from Mm as red markers
        ax.scatter(Mm[:, 0], Mm[:, 1], Mm[:, 2], c='red', marker='o', label='Motors')

        # Connect each pair of points with a line
        for i in range(len(p0)):
            ax.plot([p0[i, 0], Mm[i, 0]],
                    [p0[i, 1], Mm[i, 1]],
                    [p0[i, 2], Mm[i, 2]], c='green')
        # Define the length of the unit vectors
        axis_length = 0.05
        # Plot the unit vectors at the EE
        ex = R_ee @ np.array([[axis_length], [0], [0]])
        ey = R_ee @ np.array([[0], [axis_length], [0]])
        ez = R_ee @ np.array([[0], [0], [axis_length]])

        ax.plot([p_ee[0], p_ee[0] + ex[0,0]], [p_ee[1], p_ee[1]+ex[1,0]], [p_ee[2], p_ee[2]+ex[2,0]],'r-')
        ax.plot([p_ee[0], p_ee[0]+ ey[0,0]], [p_ee[1], p_ee[1]+ ey[1,0]], [p_ee[2], p_ee[2]+ey[2,0]],'g-')
        ax.plot([p_ee[0], p_ee[0]+ ez[0,0]], [p_ee[1], p_ee[1]+ ez[1,0]], [p_ee[2], p_ee[2] + ez[2,0]],'b-')
        
        # Plot the axes at the base/world unit vectors
        ax.plot([0, 0 + axis_length], [0,0], [0,0], 'r-')
        ax.plot([0,0],[0, 0 + axis_length], [0,0], 'g-')
        ax.plot([0,0], [0,0], [0, 0 + axis_length], 'b-')

       # # Set axis labels
        ax.set_xlabel('x-axis (m)',fontdict=label_font,labelpad=10)
        ax.set_zlabel('z-axis (m)',fontdict=label_font,labelpad=10)
        ax.set_ylabel('y-axis (m)',fontdict=label_font,labelpad=10)

        # Set the number of ticks or specify a locator
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))  # Adjust the number (5 in this example)
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))  # Adjust the number (5 in this example)
        ax.zaxis.set_major_locator(plt.MaxNLocator(5))  # Adjust the number (5 in this example)

        # Increase the font size of tick labels
        ax.tick_params(axis='x', labelsize=14)#,pad=1)
        ax.tick_params(axis='y', labelsize=14)#,pad=1)
        ax.tick_params(axis='z', labelsize=14)#,pad=1)
        # ax.legend()

        plt.savefig(path_img)

        # # Show the plot
        plt.show()

        plt.close()

    if InteractivePlot:
        # Create an empty 3D scatter plot
        fig = go.Figure()
        # Add each set of (3, N) points as separate scatter plots
        for i, coord_set in enumerate(coordinates):
            x, y, z = coord_set[0, :], coord_set[1, :], coord_set[2, :]
            fig.add_trace(go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode='markers',
                marker=dict(size=2),
                name=f'Set {i + 1}'
            ))
        # Add points from p0 as blue markers
        fig.add_trace(go.Scatter3d(
            x=p0[:, 0],
            y=p0[:, 1],
            z=p0[:, 2],
            mode='markers',
            marker=dict(size=5, color='blue'),
            name='p0'
        ))

        # Add points from Mm as red markers
        fig.add_trace(go.Scatter3d(
            x=Mm[:, 0],
            y=Mm[:, 1],
            z=Mm[:, 2],
            mode='markers',
            marker=dict(size=5, color='red'),
            name='Mm'
        ))

        # Connect each pair of points with a line
        for i in range(len(p0)):
            fig.add_trace(go.Scatter3d(
                x=[p0[i, 0], Mm[i, 0]],
                y=[p0[i, 1], Mm[i, 1]],
                z=[p0[i, 2], Mm[i, 2]],
                mode='lines',
                line=dict(width=2, color='green'),
                showlegend=False
            ))

        # Customize the layout
        fig.update_layout(
            scene=dict(
                xaxis_title='X',
                yaxis_title='Y',
                zaxis_title='Z',
            ),
            title='Interactive 3D Scatter Plot'
        )
        # Add points as blue markers
        fig.add_trace(go.Scatter3d(
            x=r[:, 0],
            y=r[:, 1],
            z=r[:, 2],
            mode='markers',
            marker=dict(size=5, color='blue'),
            name='Points'
        ))

        # Connect each point to the next one with a line
        for i in range(len(r)):
            next_i = (i + 1) % len(r)  # Wrap around to connect the last point to the first one
            fig.add_trace(go.Scatter3d(
                x=[r[i, 0], r[next_i, 0]],
                y=[r[i, 1], r[next_i, 1]],
                z=[r[i, 2], r[next_i, 2]],
                mode='lines',
                line=dict(width=2, color='green'),
                showlegend=False
            ))

        # Save the interactive plot as an HTML file
        fig.write_html(path_img_interactive)

    return None

def save_excel(new_file_path, data_list):
    '''to save multiple arrays as a new sheet in the same generated excel file'''
    # Create a new Excel writer with a specified file path
    with pd.ExcelWriter(new_file_path, engine='openpyxl') as writer:
        # Iterate over the list of arrays and write each array as a new sheet
        for i, data in enumerate(data_list):
            sheet_name = f'Sheet{i + 1}'  # Generate a new sheet name
            df = pd.DataFrame(data)  # Convert NumPy array to a DataFrame
            df.to_excel(writer, sheet_name=sheet_name, index=False)


def RPY_to_matrix(roll, pitch, yaw):
    '''Map roll, pitch, yaw to matrix'''
    matrix = np.array([[np.cos(yaw)*np.cos(pitch),-np.sin(yaw)*np.cos(roll)+np.cos(yaw)*np.sin(pitch)*np.sin(roll),np.cos(yaw)*np.sin(pitch)*np.cos(roll)+np.sin(yaw)*np.sin(roll)],
            [np.sin(yaw)*np.cos(pitch), np.cos(yaw)*np.cos(roll)+np.sin(yaw)*np.sin(pitch)*np.sin(roll),np.sin(yaw)*np.sin(pitch)*np.cos(roll)-np.cos(yaw)*np.sin(roll)],
            [-np.sin(pitch),np.cos(pitch)*np.sin(roll),np.cos(pitch)*np.cos(roll)]])
    return matrix 



saveData_bool=True #saving data to an excel file 
plot_bool=True #plotting function
NormalPlot = True
InteractivePlot = False


#intial guess for the pose of the end-effector
p_ee = np.array([0,0,0.4])#pee for FK
R_ee = np.array([np.deg2rad(0), np.deg2rad(0), np.deg2rad(0)])

#intial guess for the Motor angle in radians
qm = np.array([0.39497946,0.39497946,0.39497946,0.39497946,0.39497946,0.39497946])  #in degrees

#universal joint angle initialization
qi = np.array([0,0,
               0,0,
               0,0,
               0,0,
               0,0,
               0,0])

#initializing the guess vector for the FK model
init_guess = np.concatenate([np.zeros(24),qi,p_ee,R_ee]) #42 variables
optimised_states,rtrack, Mm, ptrack, nLtrack, F, pLtrack,p0 = Forward_Kinetostatic(init_guess, qm)

rtrack = rtrack.reshape(rtrack.shape[0] * rtrack.shape[1], rtrack.shape[2])
print(f"shape of rtrack {rtrack.shape}")
p_ee = optimised_states[36:39]
R_ee = optimised_states[39:42]

R_ee = RPY_to_matrix(R_ee[0], R_ee[1], R_ee[2])

print(f"external force magnitude: {np.linalg.norm(F)}")
if saveData_bool: #if True then save as excel
        save_excel('/home/dfki.uni-bremen.de/vrodrigues/6_rus_Stewart_platform/Simulation/Flexible_links/Cosserat_rod_theory/excel_files/PACOMA_FK.xlsx',
                   [ptrack,rtrack, nLtrack])
if plot_bool: #if True then plot
    plotting('/home/dfki.uni-bremen.de/vrodrigues/6_rus_Stewart_platform/Simulation/Flexible_links/Cosserat_rod_theory/excel_files/PACOMA_FK.xlsx', 
            '/home/dfki.uni-bremen.de/vrodrigues/6_rus_Stewart_platform/Simulation/Flexible_links/Cosserat_rod_theory/Images/PACOMA_test_FK.html',
            '/home/dfki.uni-bremen.de/vrodrigues/6_rus_Stewart_platform/Simulation/Flexible_links/Cosserat_rod_theory/Images/PACOMA_test_FK.png')
