VTVL Rocket 6-DoF Simulator

A 6-DoF simulator and hybrid guidance, navigation and control (GNC) system for trajectory optimisation and online replanning of a small-scale vertical take-off, vertical landing (VTVL) launch vehicle.

This repository contains the MATLAB/Simulink code for my BEng Aerospace Engineering dissertation at the University of Sheffield (supervisor: Dr Yuanbo Nie, May 2026): Trajectory Optimisation and Online Replanning for a Small-Scale Vertical Take-Off Vertical Landing Launch Vehicle

The work is a pre-flight simulation study supporting the Project Sunride entry to the Sunride Lander Challenge. The vehicle in simulation is the 16.88 kg "Desert Winds" parameter set, calibrated against the Hopper demonstrator's conceptual design.


Overview

The GNC stack combines an offline optimal reference with an online convex replanner:

- Offline guidance: a four-phase fuel-optimal trajectory (ascent, hover, descent), transcribed by Hermite-Simpson collocation and solved with ICLOCS2 / IPOPT.
- Online guidance (Mode B): a sequential convex programming (SCP) model predictive controller with lossless convexification (LCvx) of the terminal thrust constraint, built with YALMIP and solved with ECOS.
- Dispatcher: a five-flag validity-gated, two-mode dispatcher with a three-strike fallback latch that arbitrates between the offline reference (Mode A) and the online replanner (Mode B).
- Navigation: a multiplicative extended Kalman filter (MEKF) with 17 nominal states and 16 error states.
- Control: a mass-gain-scheduled linear-quadratic-integral (LQI) attitude regulator.
- Plant: a 14-state quaternion 6-DoF vehicle model, verified against seven analytical benchmarks.


Repository structure
Folder	Contents
offline/	ICLOCS2 four-phase optimal control problem: driver script, per-phase problem, dynamics and settings files
guidance/	Online guidance: SCP-MPC with LCvx (SOCP build, SCP iteration, MPC solve)
nav/	State estimation (MEKF)
control/	Attitude control and thrust allocation
actuators/	TVC servo and throttle valve models
plant/	14-state quaternion 6-DoF vehicle model
sensors/	Sensor models
scheduler/	Mode dispatcher and scheduling logic
sim/	Closed-loop simulation
tests/	Verification and unit test scripts (Mode B steps, closed-loop descent, and so on)



Requirements
MATLAB (version: R2026a) with Simulink
ICLOCS2 (offline trajectory optimisation). It is not included here, so download it separately and add it to the MATLAB path.
IPOPT (used by ICLOCS2)
YALMIP
ECOS (SOCP solver)


Getting started
Clone the repository: git clone https://github.com/Lucasxwhite04/vtvl-rocket-6dof-simulator.git
Install ICLOCS2, YALMIP and ECOS, and add them to the MATLAB path.
Add this repository to the MATLAB path, including subfolders: addpath(genpath(pwd));

Thesis

The full dissertation (PDF):

Acknowledgements
Supervisor: Dr Yuanbo Nie, University of Sheffield
ICLOCS2 by Yuanbo Nie, Omar Faqir and Eric Kerrigan (Imperial College London), MIT licence
Project Sunride, University of Sheffield


Author

Lucas White, BEng Aerospace Engineering, University of Sheffield





