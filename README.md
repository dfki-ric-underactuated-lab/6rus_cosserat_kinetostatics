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

## Installation
```jl
pkg> add NovelWrist
```

## Working of the code
### Kinetostatic model 
## Inverse Kinetostatic (IK) model:
IK model can be described as $'\textbf{q}_1=IK(\textbf{p}_e,\textbf{R}_e,\\\textbf{F},\textbf{M},\textbf{B}_{IK},\textbf{H}_0)'$ where vector $'\textbf{q}_1'$ consists of actuator variables and $'\textbf{B}_{IK}'$ contains the unknown variables for the system for the IK model.
In IK, additionally $'q_{1i}$, $q_{2i}'$, and $'q_{3i}'$ is included as unknown variables for each rod such that $'\textbf{B}_{IK} \in \mathbb{R}^{42}'$, expressed as 
\(\textbf{B}_{IK}=[n_{x,1}(0), n_{y,1}(0), n_{z,1}(0),m_{z,1}(0),q_{11}, q_{21}, q_{31},\dots]^T\).

## Forward Kinetostatic (FK) model:
Similarly, for the FK model can be described as function $\textbf{p}_e, \textbf{R}_e =FK(\textbf{q}_1,\textbf{F},\textbf{M},\textbf{B}_{FK},\textbf{H}_0)$) where $\textbf{B}_{FK}$ contains the unknown variables for the system for the FK model. In this case, $\textbf{p}_e$, $\textbf{R}_e$, $q_{2i}$, and $q_{3i}$ is considered as unknown variables such that $\textbf{B}_{FK} \in \mathbb{R}^{42}$ for the coupled system is given by:
\(\textbf{B}_{FK}=[n_{x,1}(0), n_{y,1}(0), n_{z,1}(0),m_{z,1}(0),q_{21}, q_{31},\dots,\textbf{p}_e, \textbf{R}_e]^T\).

### Kinematics
#### Inverse Kinematics 
Computation of the actuator length from a given pose defined by *inclination* ($`\alpha`$) and *tilt* ($`\gamma`$). Note that 'solution' defines which intersection points to pick from both sides of the circle-sphere intersectio. All functions consider *intrinsic* rotation of the end-effector but it can be changed via `intrinsic = false`.

```jl
julia> x = [0, 0]; # angles in rad  

julia> q = inverse_kinematics(x, RH5_wrist; solution = [1,2])
2-element Vector{Real}:
 0.13347357815533836
 0.13347357815533836
```

#### Forward Kinematics
Computes the end-effector orientation `α` and `γ`, given the solution for the actuator lengths `q`:

```jl
julia> α, γ = forward_kinematics(q, RH5_wrist, solution = [2,1,1]) 
(-2.2204460492503136e-16, 0.0)
```



#### Constrained Jacobian
To get the Jacobian $`\mathbf{J}`$ as product of the inverted work space Jacobian $`\mathbf{J}_x`$ and the joint space Jacobian $`\mathbf{J}_q`$:

```jl
julia> J = Jacobian(x, RH5_wrist; specsol = [1,2] split = false)
2×2 Matrix{Real}:
 15.9633   15.9633
 17.737   -17.737
```
When `split = true`, $`\mathbf{J}_x`$ and $`\mathbf{J}_q`$ are returned componentwise. 

### Performance Analysis
#### Conditioning
The condition index of the novel mechanism can be plotted over `α` and `γ`:

```jl
julia> plot_conditioning(RH5_wrist, α = (-π, π), γ = (-π, π), solution = [1,2], resol = 500) # increasing resol will give a higher resolution
```
![test](./assets/condition_index.png?raw=true "Conditioning")
The dashed lines indicate the workspace limits imposed by `actuator_limits`.

#### Configuration Space
The actuator lengths for plotting the the configuration space are computed for end-effector orientations between -π and π: 
```jl
julia> plot_configuration_space(RH5_wrist; solution = [1,2], intrinsic = true, resol = 100)
```
![test](./assets/c_space.png?raw=true "Configuration space")
Here, for better visibility, the `actuator_limits` are visualized using a red rectangle. 

#### Comparison to Conventional Wrist Designs

Computes and plots the **difference of the condition index** between $`2S\underbar{P}U+2RSU+1U`$ and $`2S\underbar{P}U+1U`$ mechanism (positive values indicate increased dexterity of the novel design): 

```jl
julia> plot_comparative_conditioning(RH5_wrist, α = (-π, π), γ = (-π, π), solution = [1,2], resol = 400)
```
![test](./assets/conditioning_comparison.png?raw=true "Comparison of conditioning")


The **singularity curves** of novel design and comparative design are obtained by sampling through the work space. Note, that in order to get closed contures, a high value for `resol` has to be set. This however increases the computing time considerably.        

```jl
julia> plot_comparative_singularities(RH5_wrist, α = (-π, π), γ = (-π, π), solution = [1,2], intrinsic = true, resol = 5000)
```
![test](./assets/singularities_C.png?raw=true "Comparison of singularity curves")
The theoretically feasible work space for the novel design is denoted by the blue coloured "shadow".

Plots of **Torque** and **Speed** at pure inclination and pure tilt movements can be computed. Additionally, characteristic values are printed to the console:

```jl
julia> plot_torque_C(RH5_wrist, α = (-π, π), γ = (-π, π), solution = [1,2], resol=600)
    Pure inclination/tilt characteristics - new wrist:
    Inclination range: -0.74/1.83 rad, 
    Maximum inclination torque: 62.94 Nm, correspondent inclination velocity: 6.36 rad/s, 
    Tilt range: -0.97/0.98 rad, 
    Maximum tilt torque: 56.38 Nm, correspondent tilt velocity: 7.09 rad/s
s
    Pure inclination/tilt characteristics - comparative design:
    Inclination range: -0.74/1.76 rad, 
    Maximum inclination torque: 59.86 Nm, correspondent inclination velocity: 6.68 rad/s, 
    Tilt range: -0.97/0.98 rad, 
    Maximum tilt torque: 53.86 Nm, correspondent tilt velocity: 7.43 rad/s
```
![test](./assets/torque_and_speed.png?raw=true "Comparison of torque and speed at pure inclination/ tilt")

##### Acknowledgements
The work presented in this paper is supported by the PACOMA project (Grant No. ESA-TECMSM-SOW-022836) subcontracted to us by Airbus Defence \& Space GmbH (Grant No. D.4283.01.02.01) with funds from the European Space Agency. The authors also want to acknowledge John Till's GitHub [tutorial](https://github.com/JohnDTill/ContinuumRobotExamples) on PCR and his guidance on deriving the boundary condition equations for the proposed PCR.
