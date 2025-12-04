% Assignment 6
% Rob Włodarski
% December 2025

% Clear the space
clear; 
close all;

%% 1. Prepare the space and load parameters

% Add paths
addpath _functions\
addpath _scripts\

% Load parameters
LoadParameters;

%% 2. Calibrate eta & steady state
Parameters.pEta         = fnCalibrateEta(Parameters);
ResultsSS               = fnSolveSSEquilibrium(Parameters);

%% 3. Uncertainty equilibrium
Results                 = fnEquilibriumUncertainty(Parameters, Grids,ResultsSS);

% Run version w/ no constraint
Parameters.pPsi         = 0;
ResultsNC               = fnEquilibriumUncertainty(Parameters, Grids,ResultsSS);