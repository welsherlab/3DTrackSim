# 3DTrackSim
Numerical simulations of fully 3D active-feedback single-particle tracking for comparisons between combined online Bayesian and windowed estimation of background and signal localization (COBWEBS) and Kalman mapped gaussian assumed density filter localization. 
For additional conceptual information see 2D algorithm development, https://doi.org/10.1063/5.0118317 and 3D communication conference proceeding, 

MATLAB code written in MATLAB R2022a.
Requires Signal Processing Toolbox and and Statistics and Machine Learning Toolboxes

To run simple even background simulations, download galvoimpulseresponse.mat, Impulse_Response.txt, and desired tracking simulation function. 

For Kalman XY, Z download track_Kalman_3D.m For COBWEBS XY + Kalman Z download track_XYBayesZKalm.m  For COBWEBS XY + COBWEBS Z download track_XYBayesZBayes.m

All track functions share inputs:

D: Desired diffusion coefficient of particle in meters squared per second. Suggested value: 6e-12

s: max particle count rate in counts per second (cps). Suggested value range: (2.8e5, 7.0e7) 

bofr: background intensity in cps on each side of any defined radii. Suggested value range: (0,1e6)

r: radial position of background intensity changes in meters relative to initial particle position. Suggested value: Even bg-[] Uneven bg: 2e-6

ki: integral feedback control parameter along X/Y. Suggested value: 0.05

kiz: integral feedback control parameter along Z. Suggested value: 0.04

N: number of steps desired to run given trajectories. Suggested value: 5e4 (1 second trajectory)

ogtau: Bin time in seconds prior to coercion. Suggested value: 20e-6

To pass uniform background set bofr to desired value and pass null value [] for r. 

For additional documentation on outputs, originof toolbox dependencies, and suggested parameter values see the readme pdf. 
