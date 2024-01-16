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


def main(init_guess, qm):

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
    
    print(f"shape of p0 is {p0[0].shape}")

    # Independent Parametersfor the flexible link
    E = 110.3e9  # Young's modulus (N/m2)
    nu = 0.31   #tranversal elongation/axial compression
    G = E/(2*(1+nu))  # Shear modulus (N/m2)
    rho = 4428.8  # Density (kg/m3)
    h = 0.0012 #thickness of the rectangular cross-section (m)
    w = 0.077  #width of the cross-section (m)
    ee_mass = 0.000001 #Kg #neglible mass 0.000001
    g = np.array([0,0,-9.81]) 

    # Dependent Parameters for the flexible link
    A = h * w  #Cross-sectional area  # Cross-sectional area
    Ixx = (w * h**3)/12  #Area moment of inertia
    Iyy = (h * w**3)/12  # Area moment of inertia
    J = Ixx + Iyy

    Kse = np.diag([G * A, G * A, E * A])  # shear and extension
    KseInv = np.linalg.inv(Kse)
    Kbt = np.diag([E * Ixx, E * Iyy, G * J]) #bending and torsion
    KbtInv = np.linalg.inv(Kbt)
    # print(f"Kse: {Kse}")
    # print(f"Kbt: {Kbt}")

    F = np.array([0, 0, -9.81*ee_mass]) #+ np.array([10,0,0])
    M = np.zeros(3)
    # p_ee = np.array([0,0,0.4])
    #whether to use RPY or rotation matrix--> boolean argument in the main function

    ptrack = [] #tracks the position of the flexible link in space
    nLtrack = [] #track the force acting at the last rod segment
    pLtrack = [] #track the position of the last cross-section of the rod 
    rtrack =[]   #track the position of the spherical joint position
    counter = 0
    #Sub-functions
    def obj_fun(guess):
        '''objective function aka cost function'''
        # Update guessed initial conditions
        nonlocal counter #count the function calls
        counter += 1
        print(counter)
        residual = np.empty(42) #42 constraint equations
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
            rtrack.append(r)
                      
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
            # print(f"shape of r: {Y[0:3,:].shape}")
            ptrack.append(Y[0:3,:])
            # Rtrack.append(Y[3:12,:])
            
        
            # Calculate distal constraint violation
            pL_shot = Y[0:3, -1]
            pLtrack.append(Y[0:3, -1])
            nL = Y[12:15, -1]
            nLtrack.append(Y[12:15, -1])
            mL = Y[15:18, -1]
            # print(f"pL_shot: {R_ee.shape}")
            residual[6*i:6*i+3] = pL_shot - r[i]
            residual[6*i+3:6*i+6] = mL #spherical joint
            
            #F and M respectively = 0
            EF = EF - nL
            EM = EM - (mL + np.cross(r[i],nL))
        
        residual[36:39] = EF
        residual[39:42] = EM
        # print(f"residual vector: {residual}")
        print(f"residual magnitude: {np.linalg.norm(residual)}")
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
    
    def hat(y):
        '''Mapping the vector into skew-symmetric matrix'''
        return np.array([[0, -y[2], y[1]],
                         [y[2], 0, -y[0]],
                         [-y[1], y[0], 0]])

    def inv_hat(skew):
        '''mapping the skew-symmentric matrix to vector '''
        return np.array([skew[2, 1], skew[0, 2], skew[1, 0]])
        
    
    # init_guess = np.concatenate([np.zeros(30),p_ee[2]*np.ones(6)])
    Y = None
    options = {'ftol': 1e-5, 'xtol': 1e-5}
    # Record the start time
    start_time = time.time()
    sol = scipy.optimize.root(obj_fun, init_guess, method='lm', jac = False, tol=1e-5, options=options).x
    print(f"sol: {sol}")
    # sol = 0
    # Record the end time of code
    end_time = time.time()
    # Calculate the total time taken
    total_time = end_time - start_time
    # Print the total time taken
    # print(f"sol:  {sol}")
    print(f"Total time taken: {total_time} seconds")

    ptrack = np.array(ptrack)
    # p0track = np.array(p0track)
    nLtrack = np.array(nLtrack)

    # # ptrack, Rtrack = np.array(ptrack), np.array(Rtrack)
    # # print(Rtrack.shape)
    ptrack = ptrack.reshape(ptrack.shape[0] * ptrack.shape[1], ptrack.shape[2])
    # Rtrack = Rtrack.reshape(Rtrack.shape[0] , Rtrack.shape[1]*Rtrack.shape[2])
    # print(ptrack.shape)

    '''sol.x = contains the optimised guess vector which minimizes the obj_fun
    ptrack = stores the states of each leg for every iteration of the optimization
    '''
    
    return sol, np.array(rtrack), Mm, ptrack, nLtrack,F, pLtrack[-6:],p0


# def plotting(path_excel1, path_excel2, path_img):
#         '''
#         path_excel  = '/path/of/the/saved/excelfile.xlsx'
#         path_img = '/path/of/the/image.png'
#         '''
  
#         df1 = pd.read_excel(path_excel1,header=None)
#         df2 = pd.read_excel(path_excel2,header=None)
        
#         last_18_rows = df1.tail(18)
#         p0 = np.array(df2.tail(6))


#         N = 100
#         # Reshape the DataFrame to get 6 sets of (3, N) coordinates
#         coordinates = last_18_rows.values.reshape(6, 3, N)

#         # # Create a 3D plot
#         fig = plt.figure()
#         ax = fig.add_subplot(111, projection='3d')

#         # # Plot each set of (3, 3) coordinates
#         for coord_set in coordinates:
#             ax.scatter(coord_set.T[:, 0], coord_set.T[:, 1], coord_set.T[:, 2], s = 10)


#         for n, m in zip(p0, Mm):
#             plt.plot([n[0], m[0]], [n[1], m[1]], [n[2], m[2]], marker='o')
        
#         # Plot lines connecting the points in sequence
#         ax.scatter(r[:, 0], r[:, 1], r[:, 2], c='b', marker='o')

#         # Plot lines connecting the points in sequence
#         for i in range(len(r)):
#             ax.plot([r[i][0], r[(i + 1) % len(r)][0]],
#                     [r[i][1], r[(i + 1) % len(r)][1]],
#                     [r[i][2], r[(i + 1) % len(r)][2]], c='r')


#         # plt.savefig(path_img)

#         # # Show the plot
#         plt.show()

#         plt.close()

#         return None

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

    for i in range(6):
        ls_value.append(np.linalg.norm(p0[i]-pLtrack[i]))

    print(f"dist of pL from world frame = {pLtrack}")
    print(f"dist bet pL and p0 along the tension cable = {ls_value}")
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

# def save_excel(valid_array, path):
#     '''
#     valid_array = type/shape of the solution array should be valid for the below excel library
#     path  = '/path/of/the/excelfile.xlsx'
    
#     '''
#     #save the date of the states in an excel file
#     with xlsxwriter.Workbook(path) as workbook:
#         worksheet = workbook.add_worksheet()

#         for row_num, data in enumerate(valid_array):
#             worksheet.write_row(row_num, 0, data)


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

plot_bool=True #plotting function
saveData_bool=True #saving data to an excel file 
NormalPlot = True
InteractivePlot = False
ls_value = []
#Initialization for the inverse kinematics
p_ee = np.array([0.1,0,0.6])

R_ee = [np.deg2rad(0), np.deg2rad(0), np.deg2rad(0)]

#Motor angle in radians
qm = np.array([np.deg2rad(52.26244404),np.deg2rad(47.12311329), 
               np.deg2rad(54.68361852),np.deg2rad(57.2157811),
               np.deg2rad(40.48599139),np.deg2rad(43.03864202)])  #in degrees

#joint angle initialization
qi = np.array([0,0,
               0,0,
               0,0,
               0,0,
               0,0,
               0,0])

init_guess = np.concatenate([np.zeros(24),qi,p_ee,R_ee]) #42 variables
# print(f"len of guess vec: {init_guess.shape}")
# angles = np.array([np.deg2rad(0),np.deg2rad(0),np.deg2rad(0)])
#excuting the main function
optimised_states,rtrack, Mm, ptrack, nLtrack, F, pLtrack,p0 = main(init_guess, qm)

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

#tuning paramters: J, ee_mass, v, u