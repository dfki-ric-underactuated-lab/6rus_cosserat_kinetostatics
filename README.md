# 6rus_cosserat_kinetostatics
## Abstract:
Parallel Continuum Robots (PCR) are closed-loop mechanisms but use
elastic kinematic links connected in parallel between the end-effector (EE) and
the base platform. PCRs are actuated primarily through large deflections of the
interconnected elastic links unlike by rigid joints in rigid parallel mechanisms.
In this work, Cosserat rod theory-based forward and inverse kinetostatic models
of 6-RUS PCR are proposed. A set of simulations are performed to analyze the
proposed PCR structure which includes maneuverability in 3-dimensional space
through trajectory following, deformation effects due to the planar rotation of the
EE platform, and axial stiffness evaluation at the EE.

**Maintainers:**
- Vinay Rodrigues [rodriguesvinay10@gmail.com](mailto:rodriguesvinay10@gmail.com)
- Bingbin Yu [bingbin.yu@dfki.de](mailto:bingbin.yu@dfki.de)

## Introduction
![test](./Images/PCR_schematic.png?raw=false "Conceptual design of $'6-\overbar{R}US'$ Parallel Continuum Robot")

In this work, boundary conditions for both IK and FK are formulated for a 6-RUS PCR using Cosserat rod theory. A shooting method is used to iteratively solve the IVP by updating the guessed variables till the boundary value constraints are within the desired tolerance. The kinetostatic model has been analysed on a different aspect in simulation. Trajectory simulation shows the FK was able to find a solution with an error of the order $1\times10^{-7}$ under constant load condition of 5 N for a helical trajectory. Maximum load capacity and axial stiffness is estimated for the PCR by applying compressing forces at the EE. The solution for different EE rotation is studied to evalute the range of motion for the PCR. A reachable workspace is estimated for the proposed PCR using the IK model. Motor angles range for each rod are also visualised for the reachable workspace. The future work includes experimental validation of this model on the physical prototype. 


## Kinetostatic model 
### Inverse Kinetostatic (IK) model:

```py
p_ee = np.array([0,0,0.5]) #position of center of the end-effector position in world coordinate
R_ee = np.array([np.deg2rad(10),np.deg2rad(0),np.deg2rad(0)]) #orientation of the end-effector platform about x,y,z axis in world coordinate

#initializing the actuator variables + universal joint values for each rod--> q=[q1, q2, q3] 
qi = np.array([0,0,0,
               0,0,0,
               0,0,0,
               0,0,0,
               0,0,0,
               0,0,0])

##initializing the guess vector for the IK model
init_guess = np.concatenate([np.zeros(24),qi]) #42 variables

q1 = Inverse_Kinetostatic(p_ee, R_ee, init_guess)
6-element Vector{Real}:
0.48785652
0.48785652
0.42037453
0.28121425
0.28121425
0.42037453
```

### Forward Kinetostatic (FK) model:
```py

#intial guess for the pose of the end-effector
p_ee = np.array([0,0,0.4])
R_ee = np.array([np.deg2rad(0), np.deg2rad(0), np.deg2rad(0)])

#intial guess for the Motor angle in radians (values taken from IK model for p_ee=[0,0,0.5], R_ee=[0, 0, 0])
qm = np.array([0.48785652,
                0.48785652,
                0.42037453,
                0.28121425,
                0.28121425,
                0.42037453])  #in degrees

#universal joint angle initialization
qi = np.array([0,0,
               0,0,
               0,0,
               0,0,
               0,0,
               0,0])

#initializing the guess vector for the FK model
init_guess = np.concatenate([np.zeros(24),qi,p_ee,R_ee]) #42 variables
p_ee, R_ee = Forward_Kinetostatic(init_guess, qm)

Optimized pose of the EE: p_ee=[9.21906358e-09 7.30108121e-04 4.97398602e-01] and R_ee=[ 1.77036209e-01 -3.45502705e-08  2.06590012e-08]
```

### Trajectory comparison
<!-- ![GitHub Octocat GIF](./Images/helical2.gif) -->
<!-- <p align="center">
  <img src="./Images/helical2.gif" alt="Helical Trajectory following GIF">
</p> -->
<div>
  <img src="./Images/PCR_schematic.png" alt="Left Image" width="800">
  <img src="/Images/helical2.gif" alt="Helical Trajectory following GIF" width="800">
</div>



In this simulation, the FK model is validated by comparing the obtained solution of the EE position with samples from a reference helical trajectory under a constant load of 5 N at the EE. Euclidean distance is calculated to measure the error between the FK model and the reference trajectory. The error is of the order $1\times10^{-7}$ for the samples which shows the validity of the boundary conditions for the FK model for the PCR.
```py
ee_mass = 0.5             #mass of the end-effector platform (Kg)

#compute motor joint angles for samples representing the position of the EE from a
#helical trajectory using IK model
IK_vec = Inverse_Kinetostatic_traj(p_ee[i], angles, init_guess)

#IK_vec=[optimised_states, total_time] where optimised_states


#Provide these motor angles as input to FK model:
ee_mass = 0.5             #mass of the end-effector platform (Kg)
#initializing the guess vector for the FK model
init_guess = np.concatenate([np.zeros(24),qi,p_ee,R_ee]) #42 variables

FK_vec = Forward_Kinetostatic_traj(motor_angle[i], init_guess)
```
##### Acknowledgements
The work presented in this paper is supported by the PACOMA project (Grant No. ESA-TECMSM-SOW-022836) subcontracted to us by Airbus Defence \& Space GmbH (Grant No. D.4283.01.02.01) with funds from the European Space Agency. The authors also want to acknowledge John Till's GitHub [tutorial](https://github.com/JohnDTill/ContinuumRobotExamples) on PCR and his guidance on deriving the boundary condition equations for the proposed PCR.
