# 3DTrackSim
Numerical simulations of fully 3D active-feedback single-particle tracking for comparisons between combined online Bayesian and windowed estimation of background and signal localization (COBWEBS) and Kalman mapped gaussian assumed density filter localization. 
For additional conceptual information see 2D algorithm development, https://doi.org/10.1063/5.0118317 and 3D communication conference proceeding, 

MATLAB code written in MATLAB R2022a.
Requires Signal Processing Toolbox and and Statistics and Machine Learning Toolboxes

To run simple even background simulations, download galvoimpulseresponse.mat, Impulse_Response.txt, and desired tracking simulation function. 

For Kalman XY, Z download track_Kalman_3D.m For COBWEBS XY + Kalman Z download track_XYBayesZKalm.m  For COBWEBS XY + COBWEBS Z download track_XYBayesZBayes.m

All track functions share inputs:

D: Desired diffusion coefficient of particle in meters squared per second.

s: max particle count rate in counts per second (cps). 

bofr: background intensity in cps on each side of any defined radii.

r: radial position of background intensity changes in meters relative to initial particle position.

ki: integral feedback control parameter along X/Y.

kiz: integral feedback control parameter along Z.

N: number of steps desired to run given trajectories. 

ogtau: Bin time in seconds prior to coercion.

To pass uniform background set bofr to desired value and pass null value [] for r. 

For documentation on outputs, origin toolbox dependencies, and additional functionality see the wiki. 
