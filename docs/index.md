---
title: "Kinetostatic Analysis for 6RUS Parallel Continuum
Robot using Cosserat Rod Theory"

Github: https://github.com/dfki-ric-underactuated-lab/6rus_cosserat_kinetostatics
# pdf: https://www.researchgate.net/publication/376186142_Fast_Dynamic_Walking_with_RH5_Humanoid_Robot
# youtube: https://www.youtube.com/watch?v=39GL2vPedGY&ab_channel=GermanResearchCenterforArtificialIntelligence

authors:
  - {name: Vinayvivian Rodrigues, affiliation_key: [1,2] link: https://de.linkedin.com/in/vinayrodrigues.html}
  - {name: Bingbin Yu, affiliation_key: 2, link: https://robotik.dfki-bremen.de/de/ueber-uns/mitarbeiter/person/biyu01.html}
  - {name: Christoph Stoeffler, affiliation_key: 2, link: https://robotik.dfki-bremen.de/de/ueber-uns/mitarbeiter/person/chst04.html}
  - {name: Shivesh Kumar, affiliation_key: [2,3], link: https://robotik.dfki-bremen.de/de/ueber-uns/mitarbeiter/shku02.html}
  
affiliations:
  - {name: Rheinisch-Westfalische Technische Hochschule Aachen, 52062 Aachen, Germany, affiliation_key: 1, link: https://www.rwth-aachen.de/go/id/a/?lidx=1}
  - {name: Robotics Innovation Center, German Research Center for Artificial Intelligence (DFKI GmbH), 28359 Bremen, Germany, affiliation_key: 2, link: https://robotik.dfki-bremen.de/de/startseite}
  - {name: Dynamics Division, Department of Mechanics and Maritime Sciences, Chalmers University of Technology, 41296 Goteborg, Sweden, affiliation_key: 3, link: https://www.chalmers.se/en/departments/m2/}
  
---

## Introduction
Humanoid robots have the potential of becoming general purpose robots augmenting the human work-force in industries. However, they must match the agility and versatility of humans. It is particularly challenging for humanoids actuated with electric drives to achieve that as one must strive for the right balance between mass-inertial distribution in the robot as well as velocity and force transmissions in its actuation concept. In addition to optimal design of the robot, the control system must be designed to exploit the full potential of the robot. In this paper, we perform experimental investigations on the dynamic walking capabilities of a series-parallel hybrid humanoid named RH5. We demonstrate that it is possible to walk up to speeds of 0.43 m/s with a position controlled robot without full state feedback which makes it one of the fastest walking humanoids with similar size and actuation modalities.

## Experimental Design
The experimental design focuses on the RH5 humanoid capabilities for dynamic walking. RH5 \cite{rh5_ebetaer2021design} is a 2 m tall series-parallel hybrid humanoid robot with 32 degrees of freedom (DoFs) and a weight of 62.5 kg. The robot has 6 DoFs per leg, 7 DoFs on each arm, 3 DoFs for the head and 3 DoFs for the body. Head and arms, accounting for $23\%$ of the total mass, are not used in the experiments.

The robot was tested in three different experiments:
<ul>  
  <li>Fast dynamic walking using several combinations of step time and step length,</li>
  <li>Fast stepping in place using several combination of step time and step height,</li>
  <li>Long stride walking, up to 0.9m step stride.</li>
</ul>

<figure>
  <!-- 
  <img src="static/figures/rh5_fast_walk_ver2.png" width="40%" alt="">
  <img src="static/figures/rh5_long_stride_contrast.png" width="40%" alt="">
  <img src="static/figures/rh5_step_in_place.png" width="60%" alt=""> 
  <img src="static/figures/fast_walk.gif" width="50%" alt="">
  <img src="static/figures/long_stride.gif" width="20%" alt="">
  <img src="static/figures/step_in_place.gif" width="60%" alt="">
  -->
  <img src="static/figures/paper31a.png" width="50%" alt="">
  <img src="static/figures/paper31a.png" width="44%" alt="">
  <br>
  <img src="static/figures/paper31a.png" width="60%" alt="" style="margin-top: 5px;">
  <figcaption>Fig.1 - Shots of different walk experiments</figcaption>
</figure>

Each motion is generated starting from a velocity and direction input given through a joystick. These define the next footsteps to follow. Given a Zero-Moment Point (ZMP) trajectory, a corresponding Center of Mass (CoM) trajectory is obtained by solving a LQR problem modeling the dynamics as a linear inverted pendulum model.
The CoM trajectory is then stabilized online by means of admittance strategies that leverage the force/torque feedback at the ankle joints of while the overall behavior of the robot is optimized using a Task Space Inverse Dynamics (TSID) based Whole Body Controller (WBC).

<figure>
  <img src="static/figures/paper31a.png" width="100%" alt="">
  <figcaption>Fig.2 - RH5 Control Architecture</figcaption>
</figure>

## Results and Discussion

<table>
  <tbody>
    <tr align="center" valign="center">
      <td>
        <figure>
          <img src="static/figures/paper31a.png" width="95%" alt="">
          <figcaption>Fig.3 - Measured footsteps for fast walking</figcaption>
        <figure>
      </td>
      <td>
        </figure>
          <img src="static/figures/paper31a.png" width="95%" alt="">
          <figcaption>Fig.4 - Measured vertical GRF while stepping in place</figcaption>
        </figure>
      </td>
    </tr>
    <tr align="center" valign="center">
      <td>
        <figure>
          <img src="static/figures/paper31a.png" width="95%" alt="">
          <figcaption>Fig.5 - Measured CoP on fastest walk</figcaption>
        </figure>
      </td>
      <td>
        <figure>
          <img src="static/figures/paper31a.png" width="95%" alt="">
          <figcaption>Fig.6 - Measured footsteps for long stride walk</figcaption>
        </figure>
      </td>
    </tr>
  </tbody>
</table>

<p>
To summarise the work,
<ul>  
  <li> RH5 is able to walk up to 0.43 m/s which is among the fastest robots with similar size and actuation modalities.</li>
  <li>  Saturation of velocity limit for knee and hip joints while using a fraction of the available torque. Possible improvements by changing screw pitch for selected linear actuators.</li>
  <li> While using only the 77% of the total mass of the robot (no arms, no head), the used effort is below 70% of the limit. This suggests that the robot is able to mantain the same performance with a complete upper body.</li>
  <li> Future introduction of an updated upper body for manipulation capabilities and momentum compensation.</li>
</ul>
</p>